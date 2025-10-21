#' Dynamic-inner Partial Least Squares (DiPLS)
#'
#' Fit the Dynamic-inner PLS algorithm of Dong & Qin (2018) to jointly model
#' input block \eqn{X_t} and output block \eqn{Y_t}. The outer loop updates the
#' dynamic weights \eqn{w}, \eqn{q}, and FIR coefficients \eqn{\beta} following
#' Appendix A of the paper, while the inner model is estimated either as a FIR
#' filter (default) or an ARX refinement.
#'
#' @param X numeric matrix with T rows (time) and m columns (inputs).
#' @param Y numeric matrix with T rows and p columns (outputs).
#' @param n_comp integer number of dynamic latent components to extract.
#' @param s non-negative integer lag order for the inner model.
#' @param preproc_x a \code{pre_processor} object from the multivarious package for the X block.
#'   Default is \code{center()} which centers columns. Use \code{prep(pass())} for no preprocessing.
#' @param preproc_y a \code{pre_processor} object from the multivarious package for the Y block.
#'   If NULL (default), uses the same preprocessing as \code{preproc_x}.
#' @param mode character, either \code{"fir"} for the basic DiPLS inner model
#'   or \code{"arx"} to augment with AR terms on \eqn{u}.
#' @param max_iter integer, maximum outer iterations per component.
#' @param tol numeric tolerance on the maximum change of \eqn{w}, \eqn{q},
#'   and \eqn{\beta} between outer iterations.
#' @param verbose integer verbosity level (0 = quiet, 1 = per-component banner,
#'   >1 also logs iteration deltas).
#'
#' @return A \code{cross_projector} object of class \code{dipls} with fields:
#'   \itemize{
#'     \item \code{vx} - X-block weight matrix (m x n_comp)
#'     \item \code{vy} - Y-block weight matrix (p x n_comp)
#'     \item \code{preproc_x} - fitted X-block preprocessor
#'     \item \code{preproc_y} - fitted Y-block preprocessor
#'     \item \code{T} - X-block score matrix (T x n_comp)
#'     \item \code{U} - Y-block score matrix (T x n_comp), first s rows are NA
#'     \item \code{P} - X-block loadings (m x n_comp)
#'     \item \code{C} - Y-block loadings (p x n_comp)
#'     \item \code{inner} - per-component inner model coefficients (list)
#'     \item \code{R} - X-block projection matrix
#'     \item \code{lag_order} - lag parameter s
#'     \item \code{mode} - inner model mode ("fir" or "arx")
#'     \item \code{iterations} - iteration counts per component
#'   }
#'
#' @export
#' @examples
#' \dontrun{
#' set.seed(1)
#' Tn <- 400; m <- 6; p <- 2; s <- 2
#' X <- matrix(rnorm(Tn * m), Tn, m)
#' w0 <- rnorm(m); w0 <- w0 / sqrt(sum(w0^2))
#' t0 <- as.numeric(X %*% w0)
#' beta0 <- c(0.8, -0.3, 0.2)
#' u <- numeric(Tn)
#' u[(s + 1):Tn] <- as.numeric(dipls_build_Ts_cpp(t0, s) %*% beta0)
#' Q0 <- matrix(rnorm(p), p); Q0 <- Q0 / sqrt(sum(Q0^2))
#' Y <- u %o% as.numeric(Q0) + matrix(rnorm(Tn * p, sd = 0.2), Tn, p)
#'
#' fit <- dipls(X, Y, n_comp = 1, s = s, mode = "fir",
#'              preproc_x = multivarious::center(), verbose = 1)
#' Yhat <- predict(fit, X)
#' }
dipls <- function(X, Y, n_comp, s,
                  preproc_x = multivarious::center(),
                  preproc_y = NULL,
                  mode = c("fir", "arx"),
                  max_iter = 200L, tol = 1e-7, verbose = 1L) {
  mode <- match.arg(mode)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  stopifnot(is.numeric(X), is.numeric(Y))

  Tn <- nrow(X)
  if (nrow(Y) != Tn) stop("X and Y must have the same number of rows.")

  s <- as.integer(s)
  if (s < 0L) stop("s must be a non-negative integer.")
  if (Tn <= s) stop("Number of rows must exceed lag order s.")

  n_comp <- as.integer(n_comp)
  if (n_comp < 1L) stop("n_comp must be >= 1.")

  max_iter <- as.integer(max_iter)
  if (max_iter < 1L) stop("max_iter must be >= 1.")
  tol <- as.numeric(tol)

  # If preproc_y not specified, create a separate instance with same spec as preproc_x
  # This is necessary because X and Y may have different numbers of columns
  if (is.null(preproc_y)) {
    # Use the same preprocessing type as X
    # This works because multivarious::center() creates a new instance each time
    preproc_y <- multivarious::center()
  }

  # Apply preprocessing using multivarious framework
  # Fit separate preprocessors to X and Y (they may have different dimensions)
  fitted_preproc_x <- multivarious::fit(preproc_x, X)
  fitted_preproc_y <- multivarious::fit(preproc_y, Y)
  Xwork <- multivarious::transform(fitted_preproc_x, X)
  Ywork <- multivarious::transform(fitted_preproc_y, Y)
  m <- ncol(Xwork)
  p <- ncol(Ywork)

  W_scaled <- matrix(0, m, n_comp)
  Q_scaled <- matrix(0, p, n_comp)
  P_scaled <- matrix(0, m, n_comp)
  C_scaled <- matrix(0, p, n_comp)
  T_scores <- matrix(0, Tn, n_comp)
  U_scores <- matrix(NA_real_, Tn, n_comp)
  inner_models <- vector("list", n_comp)
  iter_counts <- integer(n_comp)

  build_Ys <- function(Ymat, lag) {
    Ymat[seq.int(lag + 1L, nrow(Ymat)), , drop = FALSE]
  }

  comp_names <- paste0("comp", seq_len(n_comp))
  colnames(T_scores) <- comp_names
  colnames(U_scores) <- comp_names
  colnames(W_scaled) <- comp_names
  colnames(Q_scaled) <- comp_names
  colnames(P_scaled) <- comp_names
  colnames(C_scaled) <- comp_names
  if (!is.null(colnames(X))) {
    rownames(W_scaled) <- colnames(X)
    rownames(P_scaled) <- colnames(X)
  }
  if (!is.null(colnames(Y))) {
    rownames(Q_scaled) <- colnames(Y)
    rownames(C_scaled) <- colnames(Y)
  }
  for (h in seq_len(n_comp)) {
    if (verbose) {
      message(sprintf("=== DiPLS(%s): extracting component %d/%d (s=%d) ===",
                      mode, h, n_comp, s))
    }

    Ys <- build_Ys(Ywork, s)
    if (nrow(Ys) <= 1L) stop("Not enough aligned rows for current lag order.")

    beta <- numeric(s + 1L); beta[1L] <- 1
    us <- as.numeric(Ys[, 1L])
    Xlags <- dipls_build_Xlags_cpp(Xwork, s)

    w <- rep(0, m)
    q <- rep(0, p)
    delta <- Inf
    it <- 0L

    while (it < max_iter && delta > tol) {
      it <- it + 1L

      w_new <- dipls_weight_w_cpp(Xlags, beta, us)
      nrmw <- sqrt(sum(w_new^2))
      if (!is.finite(nrmw) || nrmw <= 0) {
        stop("Encountered non-finite w during iteration.")
      }
      w_new <- w_new / nrmw

      t_all <- as.numeric(Xwork %*% w_new)
      Ts <- dipls_build_Ts_cpp(t_all, s)
      v <- as.numeric(Ts %*% beta)

      q_new <- as.numeric(crossprod(Ys, v))
      nrmq <- sqrt(sum(q_new^2))
      if (!is.finite(nrmq) || nrmq <= 0) {
        stop("Encountered non-finite q during iteration.")
      }
      q_new <- q_new / nrmq
      us <- as.numeric(Ys %*% q_new)

      beta_new <- as.numeric(crossprod(Ts, us))
      nrmb <- sqrt(sum(beta_new^2))
      if (!is.finite(nrmb) || nrmb <= 0) {
        stop("Encountered non-finite beta during iteration.")
      }
      beta_new <- beta_new / nrmb

      dw <- sqrt(sum((w - w_new)^2))
      dq <- sqrt(sum((q - q_new)^2))
      db <- sqrt(sum((beta - beta_new)^2))
      delta <- max(dw, dq, db)

      if (verbose > 1L) {
        message(sprintf("  it=%d |dw|=%.3g |dq|=%.3g |db|=%.3g",
                        it, dw, dq, db))
      }

      w <- w_new
      q <- q_new
      beta <- beta_new
    }
    iter_counts[h] <- it

    t_all <- as.numeric(Xwork %*% w)
    Ts <- dipls_build_Ts_cpp(t_all, s)
    us <- as.numeric(build_Ys(Ywork, s) %*% q)

    if (mode == "arx" && s > 0L) {
      pad_us <- c(rep(0, s), us)
      Ulags <- dipls_build_Ts_cpp(pad_us, s)[, -1, drop = FALSE]
      if (ncol(Ulags) == 0L) {
        Z <- Ts
      } else {
        Z <- cbind(Ulags, Ts)
      }
      sol <- dipls_ridge_ls_cpp(Z, matrix(us, ncol = 1L), 1e-10)
      coef_hat <- as.numeric(sol$coef)
      if (ncol(Ulags) == 0L) {
        alpha_hat <- numeric(0L)
      } else {
        alpha_hat <- coef_hat[seq_len(s)]
      }
      beta_hat <- coef_hat[seq.int(length(alpha_hat) + 1L, length(coef_hat))]
      uhat <- as.numeric(sol$fitted)
      inner_models[[h]] <- list(
        mode = "arx",
        s = s,
        alpha = alpha_hat,
        beta = beta_hat
      )
    } else {
      beta_hat <- as.numeric(dipls_ridge_ls_cpp(Ts, matrix(us, ncol = 1L), 1e-10)$coef)
      uhat <- as.numeric(Ts %*% beta_hat)
      inner_models[[h]] <- list(
        mode = "fir",
        s = s,
        beta = beta_hat
      )
    }

    den_t <- sum(t_all^2)
    p_vec <- if (den_t > 0) as.numeric(crossprod(Xwork, t_all) / den_t) else rep(0, m)
    den_u <- sum(us^2)
    c_vec <- if (den_u > 0) as.numeric(crossprod(build_Ys(Ywork, s), us) / den_u) else rep(0, p)

    Xwork <- Xwork - tcrossprod(t_all, p_vec)
    Yhat_s <- uhat %o% q
    Y_rows <- seq.int(s + 1L, Tn)
    Ywork[Y_rows, ] <- Ywork[Y_rows, ] - Yhat_s

    W_scaled[, h] <- w
    Q_scaled[, h] <- q
    P_scaled[, h] <- p_vec
    C_scaled[, h] <- c_vec
    T_scores[, h] <- t_all
    U_fill <- U_scores[, h]
    U_fill[Y_rows] <- us
    U_scores[, h] <- U_fill
  }

  # Compute projection matrix R
  pinv <- function(A, tol = NULL) {
    s <- svd(A)
    if (is.null(tol)) tol <- max(dim(A)) * .Machine$double.eps * max(s$d)
    di <- ifelse(s$d > tol, 1 / s$d, 0)
    s$v %*% (di * t(s$u))
  }
  R <- W_scaled %*% pinv(crossprod(P_scaled, W_scaled))

  # Return cross_projector object
  obj <- multivarious::cross_projector(
    vx = W_scaled,
    vy = Q_scaled,
    preproc_x = fitted_preproc_x,
    preproc_y = fitted_preproc_y,
    T = T_scores,
    U = U_scores,
    P = P_scaled,
    C = C_scaled,
    inner = inner_models,
    R = R,
    lag_order = s,
    mode = mode,
    iterations = iter_counts,
    classes = "dipls"
  )

  obj
}

#' Project X to DiPLS scores
#'
#' @param model fitted DiPLS object from \code{\link{dipls}}.
#' @param Xnew numeric matrix with the same number of columns as the training
#'   \code{X}.
#'
#' @return Matrix of latent scores (rows correspond to observations).
#' @export
dipls_scores <- function(model, Xnew) {
  stopifnot(inherits(model, "dipls"))
  # Reprocess new data using fitted preprocessor
  Xp <- multivarious::reprocess(model, as.matrix(Xnew), source="X")
  scores <- Xp %*% model$R
  if (!is.null(colnames(model$T))) {
    colnames(scores) <- colnames(model$T)
  }
  scores
}

#' Predict Y using a fitted DiPLS model
#'
#' @param object fitted \code{"dipls"} object.
#' @param newdata matrix \code{X} or a list containing an element \code{X}.
#' @param ... unused.
#'
#' @return Matrix of one-step-ahead predictions for \code{Y}.
#' @export
predict.dipls <- function(object, newdata, ...) {
  stopifnot(inherits(object, "dipls"))
  Xnew <- if (is.list(newdata)) newdata$X else newdata
  Xnew <- as.matrix(Xnew)

  # Reprocess using fitted preprocessor
  Xp <- multivarious::reprocess(object, Xnew, source="X")

  Tn <- nrow(Xp)
  p <- nrow(object$vy)
  s <- object$lag_order
  Yhat_p <- matrix(0, Tn, p)

  for (j in seq_along(object$inner)) {
    wj <- object$vx[, j]
    qj <- object$vy[, j]
    tj <- as.numeric(Xp %*% wj)
    Ts <- dipls_build_Ts_cpp(tj, s)
    uh <- rep(NA_real_, Tn)

    mdl <- object$inner[[j]]
    if (mdl$mode == "arx" && s > 0L) {
      alpha <- mdl$alpha
      beta <- mdl$beta
      for (k in seq.int(s + 1L, Tn)) {
        tb <- sum(beta * Ts[k - s, ])
        ub <- 0
        if (length(alpha)) {
          for (i in seq_len(min(s, k - 1L))) {
            ub <- ub + alpha[i] * uh[k - i]
          }
        }
        uh[k] <- tb + ub
      }
    } else {
      beta <- mdl$beta
      uh_indices <- seq.int(s + 1L, Tn)
      uh[uh_indices] <- as.numeric(Ts %*% beta)
    }

    uh[seq_len(min(s, Tn))] <- NA_real_
    contrib <- uh %o% qj
    contrib[is.na(contrib)] <- 0
    Yhat_p <- Yhat_p + contrib
  }

  # Inverse transform Y predictions
  Yhat <- multivarious::inverse_transform(object$preproc_y, Yhat_p)

  if (!is.null(rownames(object$vy))) {
    colnames(Yhat) <- rownames(object$vy)
  }
  if (s > 0L && Tn >= s) {
    Yhat[seq_len(s), ] <- NA_real_
  }
  Yhat
}

#' Compute residuals for DiPLS model
#'
#' @param object fitted \code{"dipls"} object.
#' @param newdata optional list with elements \code{X} and \code{Y}. If NULL,
#'   uses training data (if stored).
#' @param ... unused.
#'
#' @return Matrix of residuals (Y - Yhat).
#' @export
residuals.dipls <- function(object, newdata = NULL, ...) {
  stopifnot(inherits(object, "dipls"))

  if (is.null(newdata)) {
    stop("residuals.dipls requires newdata with both X and Y components")
  }

  if (!is.list(newdata) || is.null(newdata$X) || is.null(newdata$Y)) {
    stop("newdata must be a list with elements 'X' and 'Y'")
  }

  # Get predictions
  Yhat <- predict.dipls(object, newdata$X)

  # Compute residuals
  Y <- as.matrix(newdata$Y)
  resid <- Y - Yhat

  # Rows affected by lag are already NA in Yhat, so residuals will be NA too
  resid
}
