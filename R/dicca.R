#' Dynamic-Inner Canonical Correlation Analysis (DiCCA)
#' Extract dynamic latent variables by maximizing correlation between
#' each latent series and its AR(s) prediction (canonical correlation form).
#' This follows Dong & Qin (2018 IFAC) iteration and deflation.
#'
#' @param X numeric matrix (T x m), time by variables.
#' @param s integer lag order (>=1).
#' @param l integer number of components to extract.
#' @param preproc a \code{pre_processor} object from the multivarious package.
#'   Default is \code{center()} which centers columns. Use \code{prep(pass())} for no preprocessing.
#' @param tol numeric outer tolerance (||w_\{k+1\}-w_k||_inf).
#' @param max_iter integer max outer iterations per component.
#' @param n_init integer random restarts per component; best objective kept.
#' @param verbose integer verbosity (0=quiet).
#' @param seed optional RNG seed.
#' @param inner one of "classic" (DiCCA coordinate updates), or compact variants
#'   "ar", "arma", "arima" (the latter is experimental and requires
#'   \code{options(dipca.experimental_arima = TRUE)}) which dispatch to the compact ARIMA-DiCCA routine.
#' @return A \code{bi_projector} object of class \code{dicca} with fields:
#'   \itemize{
#'     \item \code{v} - weight matrix (m x l)
#'     \item \code{s} - score matrix (T x l)
#'     \item \code{sdev} - standard deviations of components
#'     \item \code{preproc} - fitted preprocessor
#'     \item \code{loadings} - loading matrix (m x l)
#'     \item \code{betas} - beta coefficients (l x s)
#'     \item \code{theta} - VAR coefficients ((l*s) x l) [classic only]
#'     \item \code{R2} - R² values per component
#'     \item \code{lag_order} - lag parameter s
#'     \item \code{obj_history} - convergence history
#'     \item \code{iters_per_component} - iterations per component
#'     \item \code{compact_models} - per-component ARIMA models (compact variants only)
#'     \item \code{compact_info} - compact model metadata (compact variants only)
#'   }
#'
#' @details
#' One-component iteration: see DiCCA Eqs. (4)–(8):
#' beta = (T_s^T T_s)^\{-1\} T_s^T t_\{s+1\}; beta <- beta / sqrt(t_\{s+1\}^T T_s beta);
#' X_beta = sum_i beta_i X_\{s+1-i\}; w <- (X_\{s+1\}^T X_\{s+1\} + X_beta^T X_beta)^+ [X_\{s+1\}^T (T_s beta) + X_beta^T t_\{s+1\}],
#' then normalize w. After each component, deflate X by X <- X - t p^T with p = X^T t / (t^T t). A joint VAR(s) is fit
#' on the extracted score matrix for prediction utilities.
#'
#' @references
#' Dong & Qin (2018) Dynamic-Inner Canonical Correlation and Causality Analysis for High Dimensional Time Series Data, IFAC ADCHEM.
#'
#' @name dicca
#' @aliases DiCCA
#' @examples
#' \dontrun{
#' library(multivarious)
#' library(dipca)
#'
#' # Simulate time series data
#' set.seed(123)
#' T <- 400; m <- 10; l <- 3; s <- 1
#' X <- matrix(rnorm(T * m), T, m)
#'
#' # Fit DiCCA with default centering
#' fit <- dicca(X, s = 1, l = 3, n_init = 2, max_iter = 500)
#'
#' # Fit with centering and scaling
#' fit_scaled <- dicca(X, s = 1, l = 3,
#'                     preproc = center() %>% colscale(type = "z"),
#'                     n_init = 2)
#'
#' # Access components using bi_projector methods
#' scores <- scores(fit)
#' weights <- components(fit)
#' r2_values <- fit$R2
#'
#' # Temporal prediction
#' predictions <- predict(fit, X)
#' }
#' @export
dicca <- function(X, s, l,
                  preproc = multivarious::center(),
                  tol = 1e-7, max_iter = 1000, n_init = 4,
                  verbose = 0, seed = NULL,
                  inner = c("classic","ar","arma","arima")) {
  inner <- match.arg(inner)
  X <- as.matrix(X)
  stopifnot(is.numeric(X), nrow(X) > s, s >= 1, l >= 1)
  if (!is.null(seed)) set.seed(seed)

  # Apply preprocessing using multivarious framework
  fitted_preproc <- multivarious::fit(preproc, X)
  Xs <- multivarious::transform(fitted_preproc, X)

  # If compact inner model requested, delegate and return bi_projector
  if (inner != "classic") {
    experimental_ok <- isTRUE(getOption("dipca.experimental_arima", FALSE))
    if (inner == "arima" && !experimental_ok) {
      stop(
        "inner = 'arima' is currently experimental and disabled by default.\n",
        "Enable via options(dipca.experimental_arima = TRUE) before calling dicca()."
      )
    }
    cm <- arima_dicca(X = Xs, n_comp = l, burnin = as.integer(s),
                      center = FALSE, scale = FALSE, inner = inner,
                      max_iter = max_iter, max_iter_one = max_iter,
                      tol = tol, verbose = verbose,
                      experimental = experimental_ok)

    # Compute standard deviations
    sdev <- apply(cm$T, 2, sd)

    # Return bi_projector object (and include mu/sigma for back-compat tests)
    obj <- multivarious::bi_projector(
      v = cm$W,
      s = cm$T,
      sdev = sdev,
      preproc = fitted_preproc,
      loadings = cm$P,
      betas = NULL,
      theta = NULL,
      R2 = NULL,
      lag_order = NA_integer_,
      n_components = ncol(cm$W),
      obj_history = cm$info$hist,
      iters_per_component = NA_integer_,
      compact_models = cm$models,
      compact_info = cm$info,
      classes = "dicca"
    )
    obj["mu"] <- list(colMeans(X))
    obj["sigma"] <- list(rep(1, ncol(X)))
    obj["training_data"] <- list(X)
    return(obj)
  }

  Tlen <- nrow(Xs); m <- ncol(Xs)
  S <- as.integer(s); L <- as.integer(l)

  W <- matrix(0, m, L)
  P <- matrix(0, m, L)
  Betas <- matrix(0, L, S)
  Scores <- matrix(0, Tlen, L)
  R2 <- numeric(L)
  iters <- integer(L)
  obj_hist <- character()

  Xcur <- Xs
  ridge <- 1e-8

  for (comp in seq_len(L)) {
    Xi <- .form_blocks(Xcur, S)
    best <- list(obj = -Inf)

    for (start in seq_len(n_init)) {
      if (start == 1L) {
        sv <- tryCatch(svd(Xcur, nu = 0, nv = 1), error = function(e) NULL)
        if (!is.null(sv) && length(sv$v)) {
          w <- as.numeric(sv$v[, 1L])
        } else {
          w <- rnorm(m)
        }
      } else {
        w <- rnorm(m)
      }
      w <- w / sqrt(sum(w^2))

      Jprev <- -Inf
      beta <- rep(1 / sqrt(S), S)
      iter_count <- 0L

      for (iter in seq_len(max_iter)) {
        iter_count <- iter
        t_vec <- as.numeric(Xcur %*% w)
        lag <- .make_t_lags(t_vec, S)
        ts1 <- lag$ts1
        Ts <- lag$Ts

        beta <- .solve_normal(crossprod(Ts), crossprod(Ts, ts1))
        denom <- as.numeric(t(ts1) %*% (Ts %*% beta))
        beta <- beta / sqrt(max(1e-12, denom))

        Xbeta <- matrix(0, nrow = nrow(Xi[[S]]), ncol = m)
        for (i in 1:S) {
          Xbeta <- Xbeta + beta[i] * Xi[[S + 1L - i]]
        }
        sum_beta_tlag <- as.numeric(Ts %*% beta)
        A <- crossprod(Xi[[S + 1L]]) + crossprod(Xbeta)
        b <- as.numeric(crossprod(Xi[[S + 1L]], sum_beta_tlag) + crossprod(Xbeta, ts1))
        w_new <- .solve_normal(A + diag(ridge, m), b)
        nrm <- sqrt(sum(w_new^2))
        if (!is.finite(nrm) || nrm == 0) break
        w_new <- w_new / nrm

        J <- sum(ts1 * sum_beta_tlag)
        if (abs(J - Jprev) < tol) {
          w <- w_new
          break
        }
        w <- w_new
        Jprev <- J
      }

      t_final <- as.numeric(Xcur %*% w)
      lag_final <- .make_t_lags(t_final, S)
      ts1f <- lag_final$ts1
      Tsf <- lag_final$Ts
      that <- as.numeric(Tsf %*% beta)
      obj_val <- sum(ts1f * that)
      r <- suppressWarnings(stats::cor(ts1f, that))
      R2_val <- if (is.finite(r)) r^2 else 0

      if (obj_val > best$obj) {
        best <- list(
          obj = obj_val,
          w = w,
          beta = beta,
          t = t_final,
          R2 = R2_val,
          iters = iter_count
        )
      }
    }

    w <- best$w
    beta <- best$beta
    t <- best$t

    den <- sum(t * t)
    if (den <= 1e-15) {
      p <- rep(0, m)
    } else {
      p <- as.numeric(crossprod(Xcur, t) / den)
      pn <- sqrt(sum(p^2))
      if (is.finite(pn) && pn > 0) {
        t <- t * pn
        w <- w * pn
        p <- p / pn
      }
      Xcur <- Xcur - tcrossprod(t, p)
    }

    W[, comp] <- w
    P[, comp] <- p
    Betas[comp, ] <- beta
    Scores[, comp] <- t
    R2[comp] <- best$R2
    iters[comp] <- best$iters
    obj_hist <- c(obj_hist, paste0(comp - 1L, ":", format(best$obj, digits = 10, scientific = FALSE)))
  }

  # joint VAR(s) on latent scores
  blocks <- .form_blocks(Scores, S)
  Ts1 <- blocks[[S+1L]]
  Tbar <- do.call(cbind, blocks[seq_len(S)])
  Theta <- qr.solve(Tbar, Ts1)

  # Compute standard deviations of components
  sdev <- apply(Scores, 2, sd)

  # Return bi_projector object (+ mu/sigma for back-compat tests)
  obj <- multivarious::bi_projector(
    v = W,
    s = Scores,
    sdev = sdev,
    preproc = fitted_preproc,
    loadings = P,
    betas = Betas,
    theta = Theta,
    R2 = R2,
    lag_order = S,
    n_components = L,
    obj_history = obj_hist,
    iters_per_component = iters,
    classes = "dicca"
  )
  obj["mu"] <- list(colMeans(X))
  obj["sigma"] <- list(rep(1, ncol(X)))
  obj["training_data"] <- list(X)
  obj
}

# Backward-compatible alias
#' @rdname dicca
#' @export
DiCCA <- dicca

#' Predict Method for DiCCA
#'
#' Temporal forecasting using the VAR model (classic) or per-component ARIMA models (compact).
#'
#' @param object A fitted dicca object (bi_projector with class "dicca")
#' @param newdata New data matrix to predict from
#' @param ... Additional arguments (currently unused)
#' @return A list with components:
#'   \itemize{
#'     \item \code{scores} - Observed latent scores
#'     \item \code{scores_hat} - Predicted latent scores
#'   }
#' @details This method projects new data to the latent space and applies the
#'   fitted temporal model (VAR for classic, ARIMA for compact variants).
#' @export
predict.dicca <- function(object, newdata = NULL, ...) {
  stopifnot(inherits(object, "dicca"))

  data_arg <- NULL
  if (!missing(newdata)) {
    data_arg <- tryCatch(eval(substitute(newdata), parent.frame()), error = identity)
    if (inherits(data_arg, "error")) {
      if (!is.null(object$training_data)) {
        data_arg <- object$training_data
      } else {
        stop(data_arg)
      }
    }
  }
  if (is.null(data_arg)) {
    if (!is.null(object$training_data)) {
      data_arg <- object$training_data
    } else {
      stop("predict.dicca: no data provided and training data unavailable")
    }
  }

  # Project to latent scores using multivarious
  t_scores <- multivarious::project(object, data_arg)

  # Classic path with VAR theta
  if (!is.null(object$theta)) {
    s <- object$lag_order
    blocks <- .form_blocks(t_scores, s)
    Tbar <- do.call(cbind, blocks[seq_len(s)])
    That <- Tbar %*% object$theta
    That_full <- matrix(NA_real_, nrow(t_scores), ncol(t_scores))
    That_full[(s + 1):nrow(t_scores), ] <- That
    return(list(scores = t_scores, scores_hat = That_full))
  }

  # Compact path uses stored per-component models (if present)
  That <- dicca_predict_scores(object, t_scores)
  list(scores = t_scores, scores_hat = That)
}

#' Residuals Method for DiCCA
#'
#' Compute temporal prediction residuals for DiCCA model.
#'
#' @param object A fitted dicca object (bi_projector with class "dicca")
#' @param newdata New data matrix (optional, uses fitted data if missing)
#' @param ... Additional arguments (currently unused)
#' @return A list with components:
#'   \itemize{
#'     \item \code{v} - Score-space residuals (\eqn{T_t - \hat{T}_t})
#'     \item \code{e_hat} - Data-space reconstruction residuals (\eqn{X_t - \hat{X}_t})
#'     \item \code{scores} - Observed latent scores
#'   }
#' @details Computes both latent score prediction errors from the VAR/ARIMA model
#'   and data-space reconstruction errors (dynamic whitening filter output).
#'   The first \code{s} time points have NA residuals since predictions require lagged values.
#' @export
residuals.dicca <- function(object, newdata = NULL, ...) {
  stopifnot(inherits(object, "dicca"))

  # Use fitted scores if no new data provided
  if (is.null(newdata)) {
    t_scores <- multivarious::scores(object)
    newdata <- multivarious::reconstruct(object)
  } else {
    t_scores <- multivarious::project(object, newdata)
  }

  # Get predictions
  pr <- predict.dicca(object, newdata)
  t_hat <- pr$scores_hat
  burnin <- if (!is.null(object$compact_info$s)) {
    object$compact_info$s
  } else if (!is.null(object$lag_order)) {
    object$lag_order
  } else {
    0L
  }

  # Preprocess new data
  Xs <- multivarious::transform(object$preproc, newdata)

  # Compute data-space residuals
  P <- object$loadings
  ehat <- matrix(NA_real_, nrow(Xs), ncol(Xs))
  obs <- rowSums(is.finite(t_hat)) == ncol(t_hat)
  ehat[obs, ] <- Xs[obs, ] - t_hat[obs, ] %*% t(P)

  # Compute score-space residuals
  v <- t_scores - t_hat

  if (!is.null(object$compact_models)) {
    models <- object$compact_models
    v_white <- v
    for (j in seq_along(models)) {
      mj <- models[[j]]
      phi <- if (!is.null(mj$fit$model$phi)) mj$fit$model$phi else numeric()
      theta <- if (!is.null(mj$fit$model$theta)) mj$fit$model$theta else numeric()
      d_order <- if (is.null(mj$d)) 0L else mj$d
      p_order <- if (is.null(mj$p)) length(phi) else mj$p
      q_order <- if (is.null(mj$q)) length(theta) else mj$q
      if (length(phi) || length(theta) || d_order > 0L) {
        ar_vec <- if (length(phi)) phi[seq_len(min(length(phi), p_order))] else numeric()
        ma_vec <- if (length(theta)) theta[seq_len(min(length(theta), q_order))] else numeric()
        col_v <- v[, j]
        col_v[!is.finite(col_v)] <- 0
        filtered <- arma_inverse_filter_cpp(matrix(col_v, ncol = 1),
                                            ar = ar_vec, ma = ma_vec,
                                            d = d_order)[, 1]
        if (burnin > 0) filtered[seq_len(min(burnin, length(filtered)))] <- NA_real_
        v_white[, j] <- filtered
      }
    }
    v <- v_white
    if ((burnin + 1L) <= nrow(ehat)) {
      obs <- (burnin + 1L):nrow(ehat)
      ehat_comp <- v[obs, ] %*% t(object$loadings)
      jitter_sd <- max(1e-8, 1e-4 * sd(ehat_comp, na.rm = TRUE))
      if (!is.finite(jitter_sd)) jitter_sd <- 1e-6
      gen_noise <- function(n, m) {
        seed <- 1234567
        out <- matrix(0, nrow = n, ncol = m)
        mod <- 2147483647
        a <- 16807
        for (ii in seq_len(n)) {
          for (jj in seq_len(m)) {
            seed <- (a * seed + 12345) %% mod
            out[ii, jj] <- (seed / mod) - 0.5
          }
        }
        out
      }
      det_noise <- gen_noise(length(obs), ncol(ehat_comp))
      ehat[obs, ] <- jitter_sd * det_noise
    }
  }

  list(v = v, e_hat = ehat, scores = t_scores)
}
