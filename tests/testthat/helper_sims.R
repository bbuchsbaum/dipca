# Simulation and diagnostic helpers for DiPCA/DiCCA tests
# Adapted from Dong & Qin (2018) test suite

set.seed(123)

# ---- Diagnostic metrics -----------------------------------------------------

#' ACF Energy Metric
#'
#' Sum of squared autocorrelations at lags 1:lag.max, averaged over columns.
#' Lower values indicate whiter (less autocorrelated) residuals.
#'
#' @param x Matrix or vector of residuals
#' @param lag.max Maximum lag to include
#' @return Scalar energy metric
acf_energy <- function(x, lag.max = 10) {
  x <- as.matrix(x)
  mean(apply(x, 2, function(col) {
    ac <- acf(col, plot = FALSE, lag.max = lag.max)$acf
    sum(ac[-1]^2)  # Exclude lag 0
  }))
}

#' Mean Squared Prediction Error
#'
#' @param y_true True values (matrix)
#' @param y_pred Predicted values (matrix)
#' @return Scalar MSPE
mspe <- function(y_true, y_pred) {
  mean(rowSums((y_true - y_pred)^2))
}

#' One-Step AR Prediction for Scores
#'
#' Predict scores using AR(s) model with fitted betas.
#' Adapted to work with matrix betas (l x s) instead of list.
#'
#' @param scores N x l matrix of latent scores
#' @param betas_matrix l x s matrix of AR coefficients (each row = one component)
#' @param s Lag order
#' @return N x l matrix of predicted scores (first s rows are 0)
one_step_predict_scores_AR <- function(scores, betas_matrix, s) {
  N <- nrow(scores)
  l <- ncol(scores)
  t_hat <- matrix(0, N, l)

  for (j in 1:l) {
    beta_j <- betas_matrix[j, ]  # s-vector for component j
    for (k in (s + 1):N) {
      # AR(s): t_hat[k] = sum_{i=1}^s beta_i * t[k-i]
      t_hat[k, j] <- sum(beta_j * scores[(k - 1):(k - s), j])
    }
  }
  t_hat
}

# ---- Simulation functions ---------------------------------------------------

#' Simulate VAR(1) Latent Process with Static Loadings
#'
#' Generates dynamic, stationary data following:
#' t_k = A * t_{k-1} + v_k, v_k ~ N(0, sigma_v^2 I)
#' x_k = P * t_k + e_k, e_k ~ N(0, sigma_e^2 I)
#'
#' @param N Number of time points
#' @param m Number of observed variables
#' @param h Number of latent components
#' @param A_diag Diagonal elements of VAR(1) transition matrix (length h)
#' @param sigma_v Standard deviation of latent noise
#' @param sigma_e Standard deviation of observation noise
#' @return List with X (data), T (latent scores), P (loadings)
sim_VAR <- function(N = 300, m = 8, h = 3,
                   A_diag = c(0.8, 0.6, 0.3),
                   sigma_v = 0.3, sigma_e = 0.2) {
  stopifnot(length(A_diag) == h)

  # Orthonormal loadings P (m x h)
  Q <- qr.Q(qr(matrix(rnorm(m * h), m, h)))
  P <- Q

  # Generate latent VAR(1) process with diagonal transition
  Tlat <- matrix(0, N, h)
  for (k in 2:N) {
    Tlat[k, ] <- A_diag * Tlat[k - 1, ] + rnorm(h, 0, sigma_v)
  }

  # Observed data X = T * P' + noise
  X <- Tlat %*% t(P) + matrix(rnorm(N * m, 0, sigma_e), N, m)

  list(X = X, T = Tlat, P = P)
}

#' Simulate ARIMA Latent Process
#'
#' First factor is I(1), others are stationary ARMA processes.
#' Useful for testing ARIMA-DiCCA vs AR-DiCCA.
#'
#' @param N Number of time points
#' @param m Number of observed variables
#' @param h Number of latent components
#' @param d Integration orders for each component (length h)
#' @param ar List of AR coefficients for each component (length h)
#' @param ma List of MA coefficients for each component (length h)
#' @param sigma_v Standard deviation of latent innovations
#' @param sigma_e Standard deviation of observation noise
#' @return List with X (data), T (latent scores), P (loadings)
sim_ARIMA_latent <- function(N = 300, m = 8, h = 3,
                             d = c(1, 0, 0),  # First is I(1), others I(0)
                             ar = list(0.7, 0.6, 0.3),
                             ma = list(numeric(0), 0.4, 0.0),
                             sigma_v = 0.3, sigma_e = 0.2) {

  # Orthonormal loadings
  Q <- qr.Q(qr(matrix(rnorm(m * h), m, h)))
  P <- Q

  Tlat <- matrix(0, N, h)

  for (j in 1:h) {
    # Simulate ARIMA(p, d, q) by differencing d then ARMA
    e <- rnorm(N, 0, sigma_v)
    z <- numeric(N)

    # AR part
    p_j <- length(ar[[j]])
    q_j <- length(ma[[j]])

    for (t in 2:N) {
      ar_term <- if (p_j > 0) {
        sum(ar[[j]] * rev(z[max(1, t - p_j):(t - 1)]))
      } else 0

      ma_term <- if (q_j > 0) {
        sum(ma[[j]] * rev(e[max(1, t - q_j):(t - 1)]))
      } else 0

      z[t] <- ar_term + ma_term + e[t]
    }

    # Integrate (cumsum) d[j] times
    tj <- z
    if (d[j] > 0) {
      for (ii in 1:d[j]) {
        tj <- cumsum(tj)
      }
    }

    Tlat[, j] <- tj
  }

  # Observed data
  X <- Tlat %*% t(P) + matrix(rnorm(N * m, 0, sigma_e), N, m)

  list(X = X, T = Tlat, P = P)
}
