#' Compact DiCCA with AR/ARMA/ARIMA inner models (measurement-space update)
#'
#' Implements the measurement-space DiCCA alternating scheme with an inner
#' AR/ARMA/ARIMA model fitted on each component. The deflation step follows
#' Dong & Qin (2018): \eqn{X \leftarrow X - t p^\top} with
#' \eqn{p = X^\top t / (t^\top t)}.
#'
#' @param X numeric matrix (T x m), rows are time, columns variables.
#' @param n_comp number of dynamic latent variables to extract.
#' @param burnin integer lag order \eqn{s} (>=1). Historically the compact
#'   routine used a burn-in; we treat this as the lag order for consistency with
#'   the main \code{dicca()} function.
#' @param center logical; center columns before fitting.
#' @param scale logical; scale to unit variance if TRUE.
#' @param inner one of \code{"ar"}, \code{"arma"}, or \code{"arima"} (experimental; see \code{experimental} argument).
#' @param order_mode how to handle ARIMA orders across iterations:
#'   \code{"select_once"} (default) selects once per component and reuses;
#'   \code{"fixed"} requires the \code{order} argument; \code{"auto_each"}
#'   re-selects at every iteration.
#' @param max_iter integer maximum outer iterations for each component.
#' @param max_iter_one retained for backwards compatibility; the effective
#'   iteration cap is \code{min(max_iter, max_iter_one)}.
#' @param tol numeric convergence tolerance on the objective.
#' @param order optional list(p=, d=, q=) when \code{order_mode = "fixed"}.
#' @param pmax,qmax upper bounds when auto-selecting orders.
#' @param verbose integer verbosity (>=0).
#' @param fit_every integer; refit the inner model every \code{fit_every} iterations (default: Inf, fit only once).
#' @param experimental logical; when \code{TRUE} enables the experimental
#'   \code{inner = "arima"} pathway, otherwise requesting \code{"arima"} will
#'   error unless the option \code{options(dipca.experimental_arima = TRUE)}
#'   has been set.
#' @return a list with fields: W, P, T, models, info, R, method;
#'   with class \code{'dicca_model'}.
#' @export
arima_dicca <- function(
  X, n_comp, burnin = 0, center = TRUE, scale = FALSE,
  inner = c("arma","arima","ar"),
  order_mode = c("select_once","fixed","auto_each"),
  max_iter = 200L, max_iter_one = 100L, tol = 1e-6,
  order = NULL, pmax = 5L, qmax = 5L, verbose = 1L,
  fit_every = Inf,
  experimental = FALSE
) {
  stopifnot(is.matrix(X), nrow(X) > burnin + 1L, n_comp >= 1)
  s <- as.integer(burnin)
  if (s < 1L) stop("burnin must be >= 1 for compact DiCCA.")

  inner <- match.arg(inner)
  experimental <- isTRUE(experimental) || isTRUE(getOption("dipca.experimental_arima", FALSE))
  if (inner == "arima" && !experimental) {
    stop(
      "inner = 'arima' is currently experimental and disabled by default.\n",
      "Set experimental = TRUE or options(dipca.experimental_arima = TRUE) to enable."
    )
  }
  order_mode <- match.arg(order_mode)
  iter_cap <- min(as.integer(max_iter), as.integer(max_iter_one))

  `%||%` <- function(a, b) if (is.null(a)) b else a

  pre <- .center_scale(X, center = center, scale = scale)
  Xw <- pre$X
  Tn <- nrow(Xw); m <- ncol(Xw)

  W <- matrix(0, m, n_comp)
  P <- matrix(0, m, n_comp)
  Tmat <- matrix(0, Tn, n_comp)
  models <- vector("list", n_comp)
  obj_hist <- vector("list", n_comp)

  ridge <- 1e-8
  n_starts <- 2L

  fit_inner <- function(ts, inner, order, pmax, qmax, cache_order = NULL) {
    ts <- as.numeric(ts)
    has_forecast <- requireNamespace("forecast", quietly = TRUE)

    fitted_from_model <- function(fit) {
      fv <- try(stats::fitted(fit), silent = TRUE)
      if (!inherits(fv, "try-error") && length(fv) > 0) {
        return(as.numeric(fv))
      }
      if (!is.null(fit$residuals)) {
        res <- as.numeric(fit$residuals)
        if (length(res) == length(ts)) {
          return(ts - res)
        }
      }
      rep(NA_real_, length(ts))
    }

    fit_arima_reuse <- function(p, d, q) {
      if (has_forecast) {
        fit <- tryCatch(
          forecast::Arima(ts, order = c(p, d, q), method = "CSS-ML"),
          error = function(e1) tryCatch(
            forecast::Arima(ts, order = c(p, d, q), method = "ML"),
            error = function(e2) forecast::Arima(ts, order = c(p, d, q), method = "CSS")
          )
        )
      } else {
        fit <- tryCatch(
          stats::arima(ts, order = c(p, d, q), method = "ML"),
          error = function(e1) tryCatch(
            stats::arima(ts, order = c(p, d, q), method = "CSS-ML"),
            error = function(e2) stats::arima(ts, order = c(p, d, q), method = "CSS")
          )
        )
      }
      fit
    }

    if (!is.null(order)) {
      p <- order$p %||% 0L; d <- order$d %||% 0L; q <- order$q %||% 0L
      fit <- fit_arima_reuse(p, d, q)
      fitted_vals <- fitted_from_model(fit)
      arma <- if (!is.null(fit$arma)) fit$arma else .as_arma_meta(p = p, d = d, q = q)
      return(list(fit = fit, arma = arma, p = p, d = d, q = q,
                  order = c(p, d, q), fitted = fitted_vals))
    }

    if (!is.null(cache_order)) {
      p <- cache_order[1]; d <- cache_order[2]; q <- cache_order[3]
      fit <- fit_arima_reuse(p, d, q)
      fitted_vals <- fitted_from_model(fit)
      arma <- if (!is.null(fit$arma)) fit$arma else .as_arma_meta(p = p, d = d, q = q)
      return(list(fit = fit, arma = arma, p = p, d = d, q = q,
                  order = c(p, d, q), fitted = fitted_vals))
    }

    if (inner == "ar") {
      arfit <- try(stats::ar(ts, method = "yule-walker", order.max = pmax, aic = TRUE), silent = TRUE)
      p <- if (inherits(arfit, "try-error")) 1L else arfit$order
      fit <- fit_arima_reuse(p, 0L, 0L)
      fitted_vals <- fitted_from_model(fit)
      arma <- .as_arma_meta(p = p, d = 0L, q = 0L)
      return(list(fit = fit, arma = arma, p = p, d = 0L, q = 0L,
                  order = c(p, 0L, 0L), fitted = fitted_vals))
    }

    if (has_forecast) {
      if (inner == "arma") {
        fit <- forecast::auto.arima(
          ts, d = 0, seasonal = FALSE, stepwise = TRUE,
          approximation = FALSE, ic = "bic",
          max.p = pmax, max.q = qmax,
          allowdrift = FALSE, allowmean = FALSE
        )
      } else {
        fit <- forecast::auto.arima(
          ts, seasonal = FALSE, stepwise = TRUE,
          approximation = FALSE, ic = "bic",
          max.p = pmax, max.q = qmax
        )
      }
      arma <- fit$arma
      p <- arma[1]; q <- arma[2]; d <- arma[6]
      fitted_vals <- fitted_from_model(fit)
      return(list(fit = fit, arma = arma, p = p, d = d, q = q,
                  order = c(p, d, q), fitted = fitted_vals))
    }

    best <- list(aic = Inf)
    dgrid <- if (inner == "arma") 0L else 0:1
    for (d in dgrid) for (p in 0:pmax) for (q in 0:qmax) {
      if (p == 0L && q == 0L && d == 0L) next
      fit <- try(fit_arima_reuse(p, d, q), silent = TRUE)
      if (inherits(fit, "try-error")) next
      aic_val <- suppressWarnings(AIC(fit))
      if (length(aic_val) == 0L || !is.finite(aic_val)) next
      if (aic_val < best$aic) best <- list(aic = aic_val, fit = fit, p = p, d = d, q = q)
    }
    fit <- best$fit %||% fit_arima_reuse(1L, 0L, 0L)
    p <- best$p %||% 1L
    d <- best$d %||% 0L
    q <- best$q %||% 0L
    fitted_vals <- fitted_from_model(fit)
    arma <- if (!is.null(fit$arma)) fit$arma else .as_arma_meta(p = p, d = d, q = q)
    list(fit = fit, arma = arma, p = p, d = d, q = q,
         order = c(p, d, q), fitted = fitted_vals)
  }

  align_fitted <- function(fitted_vals, target_len, fallback) {
    if (length(fitted_vals) == 0L) return(fallback)
    if (length(fitted_vals) >= target_len) {
      fitted_vals <- tail(fitted_vals, target_len)
    } else if (length(fitted_vals) > 0L) {
      fitted_vals <- c(rep(NA_real_, target_len - length(fitted_vals)), fitted_vals)
    }
    if (length(fitted_vals) != target_len) fitted_vals <- rep(NA_real_, target_len)
    idx <- !is.finite(fitted_vals)
    fitted_vals[idx] <- fallback[idx]
    fitted_vals
  }

  for (h in seq_len(n_comp)) {
    if (verbose) message(sprintf("=== DiCCA(%s): extracting component %d/%d ===", inner, h, n_comp))
    Xi <- .form_blocks(Xw, s)
    # cache fixed cross-products for efficiency
    XtX_sp1 <- crossprod(Xi[[s + 1L]])
    G_list <- vector("list", s * s)
    for (i in 1:s) {
      for (j in 1:s) {
        G_list[[ (i - 1L) * s + j ]] <- crossprod(Xi[[s + 1L - i]], Xi[[s + 1L - j]])
      }
    }
    best <- list(obj = -Inf)
    component_order <- NULL
    fin_cache <- NULL

    for (start in seq_len(n_starts)) {
      if (start == 1L) {
        sv <- tryCatch(svd(Xw, nu = 0, nv = 1), error = function(e) NULL)
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
      hist_vals <- numeric(0)
      beta <- rep(1 / sqrt(s), s)
      iter_count <- 0L

      for (iter in seq_len(iter_cap)) {
        iter_count <- iter
        t_vec <- as.numeric(Xw %*% w)
        lag <- .make_t_lags(t_vec, s)
        ts1 <- lag$ts1
        Ts <- lag$Ts
        # Fit inner model sparingly: first iter and then every fit_every
        if (iter == 1L || (is.finite(fit_every) && (iter %% fit_every) == 0L)) {
          fin <- switch(
            order_mode,
            fixed = {
              if (is.null(order)) stop("order_mode='fixed' requires order=list(p=,d=,q=).")
              fit_inner(ts1, inner = inner, order = order, pmax = pmax, qmax = qmax, cache_order = NULL)
            },
            select_once = {
              fit_inner(ts1, inner = inner, order = NULL, pmax = pmax, qmax = qmax, cache_order = component_order)
            },
            auto_each = {
              fit_inner(ts1, inner = inner, order = NULL, pmax = pmax, qmax = qmax, cache_order = NULL)
            }
          )
          fin_cache <- fin
          if (order_mode == "select_once" && is.null(component_order)) component_order <- fin$order
        } else {
          fin <- fin_cache
        }

        beta <- .solve_normal(crossprod(Ts), crossprod(Ts, ts1))
        denom <- as.numeric(t(ts1) %*% (Ts %*% beta))
        beta <- beta / sqrt(max(1e-12, denom))

        sum_beta_tlag <- as.numeric(Ts %*% beta)
        # use fitted predictor if available, else AR proxy
        that_use <- if (!is.null(fin) && length(fin$fitted) > 0) align_fitted(fin$fitted, length(ts1), sum_beta_tlag) else sum_beta_tlag

        # A = X_{s+1}'X_{s+1} + Xbeta'Xbeta computed via cached G_list
        XbetaTXbeta <- matrix(0, m, m)
        for (i in 1:s) for (j in 1:s) {
          XbetaTXbeta <- XbetaTXbeta + (beta[i] * beta[j]) * G_list[[ (i - 1L) * s + j ]]
        }
        A <- XtX_sp1 + XbetaTXbeta
        # b = X_{s+1}' (Ts beta) + X_beta' ts1
        xb_t <- numeric(m)
        for (i in 1:s) xb_t <- xb_t + beta[i] * as.numeric(crossprod(Xi[[s + 1L - i]], ts1))
        b <- as.numeric(crossprod(Xi[[s + 1L]], sum_beta_tlag) + xb_t)
        w_new <- .solve_normal(A + diag(ridge, m), b)
        nrm <- sqrt(sum(w_new^2))
        if (!is.finite(nrm) || nrm == 0) break
        w_new <- w_new / nrm

        J <- sum(ts1 * that_use)
        hist_vals <- c(hist_vals, J)
        if (abs(J - Jprev) < tol) {
          w <- w_new
          break
        }
        w <- w_new
        Jprev <- J
      }

      t_final <- as.numeric(Xw %*% w)
      lag_final <- .make_t_lags(t_final, s)
      ts1f <- lag_final$ts1
      Tsf <- lag_final$Ts
      sum_beta_tlag <- as.numeric(Tsf %*% beta)
      r <- suppressWarnings(stats::cor(ts1f, sum_beta_tlag))
      R2_val <- if (is.finite(r)) r^2 else 0
      obj_val <- sum(ts1f * sum_beta_tlag)

      if (obj_val > best$obj) {
        best <- list(
          obj = obj_val,
          w = w,
          beta = beta,
          t = t_final,
          R2 = R2_val,
          iters = iter_count,
          hist = hist_vals,
          order = if (order_mode == "auto_each") fin$order else component_order
        )
      }
    }

    w <- best$w
    beta <- best$beta
    t <- best$t

    denom <- sum(t * t)
    if (denom <= 1e-15) {
      p_vec <- rep(0, m)
    } else {
      p_vec <- as.numeric(crossprod(Xw, t) / denom)
      pn <- sqrt(sum(p_vec^2))
      if (is.finite(pn) && pn > 0) {
        t <- t * pn
        w <- w * pn
        p_vec <- p_vec / pn
      }
      Xw <- Xw - tcrossprod(t, p_vec)
    }

    lag_store <- .make_t_lags(t, s)
    final_order <- best$order
    fin <- fit_inner(
      lag_store$ts1,
      inner = inner,
      order = if (order_mode == "fixed") order else NULL,
      pmax = pmax,
      qmax = qmax,
      cache_order = final_order
    )

    W[, h] <- w
    P[, h] <- p_vec
    Tmat[, h] <- t
    models[[h]] <- list(
      arma = if (!is.null(fin$arma)) fin$arma else .as_arma_meta(fin$p, fin$d, fin$q),
      p = fin$p,
      d = fin$d,
      q = fin$q,
      fit = fin$fit,
      fitted = fin$fitted,
      scores = t,
      innov = tryCatch(as.numeric(stats::residuals(fin$fit)), error = function(e) numeric())
    )
    obj_hist[[h]] <- best$hist
  }

  PtW <- crossprod(P, W)
  k <- ncol(PtW)
  Rmat <- W %*% .solve_normal(PtW, diag(1, k))

  # Unscale back to original units
  W_un <- sweep(W, 1, pre$scale, "/")
  P_un <- sweep(P, 1, 1 / pre$scale, "*")
  PtW_un <- crossprod(P_un, W_un)
  k_un <- ncol(PtW_un)
  R_un <- W_un %*% .solve_normal(PtW_un, diag(1, k_un))

  out <- list(
    W = W_un,
    P = P_un,
    T = Tmat,
    models = models,
    info = list(
      burnin = s,
      center = pre$center,
      scale = pre$scale,
      s = s,
      hist = obj_hist,
      inner = inner,
      order_mode = order_mode
    ),
    R = R_un,
    method = paste0("DiCCA(", inner, ")")
  )
  class(out) <- c("dicca_model", class(out))
  out
}

#' Transform new data to Di(L)V scores using a fitted compact DiCCA model
#' @param model a fitted DiCCA model object
#' @param Xnew numeric matrix of new data to transform
#' @return matrix of scores
#' @export
dicca_scores <- function(model, Xnew) {
  stopifnot(inherits(model, "dicca_model"))
  Xnew <- sweep(as.matrix(Xnew), 2, model$info$center, "-")
  Xnew <- sweep(Xnew, 2, model$info$scale, "/")
  Xnew %*% model$R
}

#' One-step-ahead prediction of DLVs (diagonal G(B))
#' @param model a fitted DiCCA model object with compact inner models
#' @param scores numeric matrix of observed scores to predict from
#' @return matrix of predicted scores
#' @export
dicca_predict_scores <- function(model, scores) {
  Tmat <- as.matrix(scores)
  Tn <- nrow(Tmat); H <- ncol(Tmat)
  That <- matrix(NA_real_, Tn, H)

  `%||%` <- function(a, b) if (is.null(a)) b else a

  predict_with_recursion <- function(model_j, new_scores) {
    fit <- model_j$fit
    if (is.null(fit) || !inherits(fit, c("Arima","ARIMA","arima"))) {
      return(rep(0, length(new_scores)))
    }
    y_train <- model_j$scores
    if (length(y_train) == 0) return(rep(0, length(new_scores)))
    p <- model_j$p %||% 0L
    q <- model_j$q %||% 0L
    d <- model_j$d %||% 0L
    phi <- if (p > 0) fit$model$phi[seq_len(p)] else numeric(0)
    theta <- if (q > 0) fit$model$theta[seq_len(q)] else numeric(0)
    y_hist <- as.numeric(y_train)
    diff_series <- function(y, d) {
      if (d <= 0) return(as.numeric(y))
      out <- as.numeric(y)
      for (i in seq_len(d)) {
        out <- diff(out)
      }
      out
    }
    z_hist <- diff_series(y_hist, d)
    e_hist <- model_j$innov
    if (length(e_hist) < length(z_hist)) {
      e_hist <- c(rep(0, length(z_hist) - length(e_hist)), e_hist)
    } else if (length(e_hist) > length(z_hist)) {
      e_hist <- tail(e_hist, length(z_hist))
    }
    if (!length(e_hist)) e_hist <- rep(0, length(z_hist))
    y_pred <- numeric(length(new_scores))
    for (i in seq_along(new_scores)) {
      if (length(z_hist) < p) {
        ar_term <- if (length(z_hist) == 0) 0 else sum(phi[seq_len(length(z_hist))] * rev(z_hist))
      } else {
        ar_term <- if (p > 0) sum(phi * rev(tail(z_hist, p))) else 0
      }
      ma_term <- if (q > 0 && length(e_hist) >= q) sum(theta * rev(tail(e_hist, q))) else 0
      z_hat <- ar_term + ma_term
      if (d == 0) {
        y_hat <- z_hat
      } else {
        xi <- tail(y_hist, d)
        y_hat <- stats::diffinv(c(rep(0, d), z_hat), differences = d, xi = xi)[d + 1L]
      }
      y_pred[i] <- y_hat

      # update histories with actual observation
      y_next <- new_scores[i]
      y_hist <- c(y_hist, y_next)
      z_actual <- if (d == 0) {
        y_next
      } else {
        diff_series(tail(y_hist, d + 1L), d)[1L]
      }
      z_hist <- c(z_hist, z_actual)
      e_new <- z_actual - z_hat
      e_hist <- c(e_hist, e_new)
    }
    y_pred
  }

  if (!is.null(model$models) && !is.null(model$info)) {
    models <- model$models
    info <- model$info
  } else if (!is.null(model$compact_models) && !is.null(model$compact_info)) {
    models <- model$compact_models
    info <- model$compact_info
  } else {
    stop("dicca_predict_scores: no compact models found on object")
  }

  burnin <- info$s %||% info$burnin %||% 0L
  for (j in seq_len(H)) {
    model_j <- models[[j]]
    fit <- model_j$fit
    fitted_vals <- model_j$fitted

    use_training_fit <- !is.null(fitted_vals) && length(fitted_vals) == Tn && Tn == length(model_j$scores)
    if (!use_training_fit) {
      fitted_vals <- predict_with_recursion(model_j, Tmat[, j])
    }

    if ((length(fitted_vals) == 0 || all(is.na(fitted_vals))) && !is.null(fit) && inherits(fit, c("Arima","ARIMA","arima"))) {
      fitted_try <- try(as.numeric(stats::fitted(fit)), silent = TRUE)
      if (!inherits(fitted_try, "try-error") && length(fitted_try) > 0) {
        fitted_vals <- fitted_try
      }
    }

    if (length(fitted_vals) > 0) {
      len <- min(length(fitted_vals), Tn - burnin)
      if (len > 0) {
        idx <- (Tn - len + 1):Tn
        That[idx, j] <- tail(fitted_vals, len)
        finite_idx <- which(is.finite(That[idx, j]))
        if (length(finite_idx) > 0) {
          first_idx <- idx[finite_idx[1]]
          if (first_idx > (burnin + 1L)) {
            That[(burnin + 1L):(first_idx - 1L), j] <- That[first_idx, j]
          }
        }
      }
    } else if ((burnin + 1L) <= Tn) {
      That[(burnin + 1L):Tn, j] <- 0
    }
    if ((burnin + 1L) <= Tn) {
      tail_idx <- (burnin + 1L):Tn
      nf <- !is.finite(That[tail_idx, j])
      if (any(nf)) That[tail_idx[nf], j] <- 0
    }
    if (burnin > 0) {
      That[seq_len(min(burnin, Tn)), j] <- NA_real_
    }
  }
  That
}

# Fallback R implementations if C++ symbols are absent -----------------------
if (!exists("arma_inverse_filter_cpp")) {
  arma_inverse_filter_cpp <- function(U, ar, ma, d) {
    X <- as.matrix(U)
    n <- nrow(X); p <- ncol(X)
    if (d > 0) {
      for (k in seq_len(d)) {
        Z <- matrix(0, n, p)
        for (j in seq_len(p)) {
          prev <- 0
          for (i in seq_len(n)) {
            curr <- X[i, j]
            Z[i, j] <- curr - prev
            prev <- curr
          }
        }
        X <- Z
      }
    }
    P <- length(ar); Q <- length(ma)
    numer <- c(1, -ar)
    denom <- c(1,  ma)
    Y <- matrix(0, n, p)
    for (col in seq_len(p)) {
      for (t in seq_len(n)) {
        ti <- t - 1L
        acc <- 0
        for (i in 0:P) {
          idx <- ti - i
          if (idx >= 0) acc <- acc + numer[i + 1L] * X[idx + 1L, col]
        }
        for (j in 1:Q) {
          yidx <- ti - j
          if (yidx >= 0) acc <- acc - denom[j + 1L] * Y[yidx + 1L, col]
        }
        Y[t, col] <- acc
      }
    }
    Y
  }
}

if (!exists("smallest_eigenvector_crossprod_cpp")) {
  smallest_eigenvector_crossprod_cpp <- function(U) {
    ev <- eigen(crossprod(U), symmetric = TRUE)
    v <- ev$vectors[, which.min(ev$values)]
    v / sqrt(sum(v^2))
  }
}
