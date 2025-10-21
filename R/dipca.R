
#' Dynamic Inner PCA (DiPCA)
#'
#' Fit dynamic inner principal components using the fast coordinate-maximization
#' algorithm with optional inner power iterations for the w-update.
#'
#' @param X numeric matrix of shape (T, m), time by variables.
#' @param s integer, lag order (>=1).
#' @param l integer, number of components to extract.
#' @param preproc a \code{pre_processor} object from the multivarious package.
#'   Default is \code{center()} which centers columns. Use \code{prep(pass())} for no preprocessing.
#' @param tol numeric, outer stopping tolerance on the eigen residual \eqn{||d - \lambda w||_\infty}.
#' @param max_iter integer, maximum outer iterations.
#' @param n_init integer, random restarts per component.
#' @param algorithm "I" (one power step) or "II" (inner power iterations).
#' @param inner_power integer, max inner power iterations (algorithm "II").
#' @param inner_tol numeric, tolerance for inner power iteration.
#' @param verbose integer, verbosity (0=quiet).
#' @param seed optional integer RNG seed for reproducibility.
#' @return A \code{bi_projector} object of class \code{dipca} with fields:
#'   \itemize{
#'     \item \code{v} - weight matrix (m x l)
#'     \item \code{s} - score matrix (T x l)
#'     \item \code{sdev} - standard deviations of components
#'     \item \code{preproc} - fitted preprocessor
#'     \item \code{loadings} - loading matrix (m x l)
#'     \item \code{betas} - beta coefficients (l x s)
#'     \item \code{theta} - VAR coefficients ((l*s) x l)
#'     \item \code{lag_order} - lag parameter s
#'     \item \code{obj_history} - convergence history
#'     \item \code{iters_per_component} - iterations per component
#'   }
#' @details
#' The implementation follows the DiPCA objective and iteration given in
#' Dong & Qin (2018, Table 2; Eqs 8–18) and the coordinate-maximization
#' interpretation and residual criterion \eqn{||d - \lambda w||_\infty} in Shin et al. (2020).
#' The code avoids forming \eqn{Y_i} explicitly; all updates are realized
#' via matvecs with lagged blocks \eqn{X_{s+1-i}}.
#'
#' @name dipca
#' @aliases DiPCA
#' @examples
#' \dontrun{
#' library(multivarious)
#' library(dipca)
#'
#' set.seed(1)
#' T <- 600; m <- 5; l <- 3; s <- 1
#'
#' # Simulate VAR(1) latent process
#' A <- matrix(c(0.5205, 0.1022, 0.0599,
#'               0.5367, -0.0139, 0.4159,
#'               0.0412, 0.6054, 0.3874), 3, 3, byrow = TRUE)
#' P <- matrix(c(0.4316, 0.1723, -0.0574,
#'               0.1202, -0.1463, 0.5348,
#'               0.2483, 0.1982, 0.4797,
#'               0.1151, 0.1557, 0.3739,
#'               0.2258, 0.5461, -0.0424), 5, 3, byrow = TRUE)
#' t <- matrix(0, T, l)
#' v <- matrix(rnorm(T * l), T, l)
#' for (k in 2:T) {
#'   t[k, ] <- c(0.5205, 0.5367, 0.0412) + A %*% t[k - 1, ] + v[k, ]
#' }
#' X <- t %*% t(P) + matrix(rnorm(T * m, sd = 0.1), T, m)
#'
#' # Fit DiPCA with centering (default)
#' fit <- dipca(X, s = 1, l = 3, n_init = 3, max_iter = 800,
#'              tol = 1e-7, algorithm = "I")
#'
#' # Fit with centering and scaling
#' fit_scaled <- dipca(X, s = 1, l = 3,
#'                     preproc = center() %>% colscale(type = "z"),
#'                     n_init = 3, max_iter = 800, tol = 1e-7)
#'
#' # Access components using bi_projector methods
#' scores <- scores(fit)           # Extract latent scores
#' weights <- components(fit)      # Extract weight matrix
#' loadings <- fit$loadings        # Extract loadings
#' theta <- fit$theta              # VAR coefficients
#'
#' # Project new data
#' X_new <- matrix(rnorm(50 * m), 50, m)
#' scores_new <- project(fit, X_new)
#'
#' # Reconstruct data
#' X_recon <- reconstruct_new(fit, X_new)
#'
#' # Temporal prediction
#' predictions <- predict(fit, X_new)
#' str(predictions)  # scores and scores_hat
#'
#' # Compute residuals
#' resid <- residuals(fit, X_new)
#' str(resid)  # v (score residuals), e_hat (data residuals), scores
#' }
#' @export
dipca <- function(X, s, l,
                  preproc = multivarious::center(),
                  tol=1e-7, max_iter=1000, n_init=4,
                  algorithm=c("I","II"), inner_power=50, inner_tol=1e-8,
                  verbose=0, seed=NULL) {
  algorithm <- match.arg(algorithm)
  X <- as.matrix(X)
  stopifnot(is.numeric(X), nrow(X) > s, s >= 1, l >= 1)
  if (!is.null(seed)) set.seed(seed)

  # Apply preprocessing using multivarious framework
  fitted_preproc <- multivarious::fit(preproc, X)
  Xs <- multivarious::transform(fitted_preproc, X)

  Tlen <- nrow(Xs)
  m <- ncol(Xs)
  S <- s; L <- l

  W <- matrix(0, m, L)
  P <- matrix(0, m, L)
  Betas <- matrix(0, L, S)
  Scores <- matrix(0, Tlen, L)
  R2 <- numeric(L)
  obj_hist <- character()
  iters <- integer(L)

  Xcur <- Xs
  for (comp in seq_len(L)) {
    Xi <- .form_blocks(Xcur, S)
    best <- list(obj = -Inf)

    for (init in seq_len(n_init)) {
      if (init == 1L) {
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

        acc <- numeric(m)
        for (i in 1:S) {
          acc <- acc +
            beta[i] * (as.numeric(crossprod(Xi[[S + 1L]], Ts[, i])) +
                         as.numeric(crossprod(Xi[[i]], ts1)))
        }
        w_new <- as.numeric(acc)
        nrm <- sqrt(sum(w_new^2))
        if (!is.finite(nrm) || nrm == 0) break
        w_new <- w_new / nrm

        J <- sum(ts1 * as.numeric(Ts %*% beta))
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

  # VAR(s) on latent scores
  blocks <- .form_blocks(Scores, S)
  Ts1 <- blocks[[S+1L]]
  Tbar <- do.call(cbind, blocks[seq_len(S)])
  Theta <- qr.solve(Tbar, Ts1)

  # Compute standard deviations of components
  sdev <- apply(Scores, 2, sd)

  # Return bi_projector object
  obj <- multivarious::bi_projector(
    v = W,
    s = Scores,
    sdev = sdev,
    preproc = fitted_preproc,
    loadings = P,
    betas = Betas,
    theta = Theta,
    lag_order = S,
    n_components = L,
    R2 = R2,
    obj_history = obj_hist,
    iters_per_component = iters,
    classes = "dipca"
  )
  obj["mu"] <- list(colMeans(X))
  obj["sigma"] <- list(rep(1, ncol(X)))
  obj["training_data"] <- list(X)
  obj
}

# Backward-compatible alias
#' @rdname dipca
#' @export
DiPCA <- dipca

#' @rdname dipca
#' @export
dipca_fit <- dipca

#' Predict Method for DiPCA
#'
#' Temporal forecasting using the VAR model fitted on latent scores.
#'
#' @param object A fitted dipca object (bi_projector with class "dipca")
#' @param newdata New data matrix to predict from
#' @param ... Additional arguments (currently unused)
#' @return A list with components:
#'   \itemize{
#'     \item \code{scores} - Observed latent scores
#'     \item \code{scores_hat} - Predicted latent scores using VAR(s) model
#'   }
#' @details This method projects new data to the latent space and applies the
#'   fitted VAR(s) model to generate one-step-ahead predictions. The first \code{s}
#'   time points have NA predictions since they require lagged values.
#' @export
predict.dipca <- function(object, newdata, ...) {
  stopifnot(inherits(object, "dipca"))
  s <- object$lag_order

  # Project to latent scores using inherited method
  t_scores <- multivarious::project(object, newdata)

  # Apply VAR model for prediction
  blocks <- .form_blocks(t_scores, s)
  Ts1 <- blocks[[s + 1L]]
  Tbar <- do.call(cbind, blocks[seq_len(s)])
  That <- Tbar %*% object$theta

  # Construct full prediction matrix with NAs for first s observations
  That_full <- matrix(NA_real_, nrow(t_scores), ncol(t_scores))
  That_full[(s + 1):nrow(t_scores), ] <- That

  list(scores = t_scores, scores_hat = That_full)
}

#' Residuals Method for DiPCA
#'
#' Compute temporal prediction residuals for DiPCA model.
#'
#' @param object A fitted dipca object (bi_projector with class "dipca")
#' @param newdata New data matrix (optional, uses fitted data if missing)
#' @param ... Additional arguments (currently unused)
#' @return A list with components:
#'   \itemize{
#'     \item \code{v} - Score-space residuals (T_t - T̂_t)
#'     \item \code{e_hat} - Data-space reconstruction residuals (X_t - X̂_t)
#'     \item \code{scores} - Observed latent scores
#'   }
#' @details Computes both latent score prediction errors from the VAR model
#'   and data-space reconstruction errors. The first \code{s} time points have
#'   NA residuals since predictions require lagged values.
#' @export
residuals.dipca <- function(object, newdata = NULL, ...) {
  stopifnot(inherits(object, "dipca"))

  # Use fitted scores if no new data provided
  if (is.null(newdata)) {
    t_scores <- multivarious::scores(object)
    newdata <- multivarious::reconstruct(object)
  } else {
    t_scores <- multivarious::project(object, newdata)
  }

  # Get predictions
  pr <- predict.dipca(object, newdata)
  t_hat <- pr$scores_hat
  s <- object$lag_order
  P <- object$loadings

  # Preprocess new data
  Xs <- multivarious::transform(object$preproc, newdata)

  # Compute data-space residuals
  ehat <- matrix(NA_real_, nrow(Xs), ncol(Xs))
  ehat[(s + 1):nrow(Xs), ] <- Xs[(s + 1):nrow(Xs), ] - t_hat[(s + 1):nrow(Xs), ] %*% t(P)

  # Compute score-space residuals
  v <- t_scores - t_hat

  list(v = v, e_hat = ehat, scores = t_scores)
}
