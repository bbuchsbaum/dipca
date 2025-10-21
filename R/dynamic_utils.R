# Internal helpers shared across DiPCA/DiCCA implementations

#' @keywords internal
.center_scale <- function(X, center = TRUE, scale = FALSE) {
  X <- as.matrix(X)
  cm <- if (isTRUE(center)) colMeans(X) else rep(0, ncol(X))
  Xc <- sweep(X, 2, cm, "-")
  cs <- rep(1, ncol(X))
  if (isTRUE(scale)) {
    cs <- sqrt(colSums(Xc^2) / max(1, (nrow(X) - 1)))
    cs[cs == 0] <- 1
    Xc <- sweep(Xc, 2, cs, "/")
  }
  list(X = Xc, center = cm, scale = cs)
}

#' @keywords internal
.form_blocks <- function(X, s) {
  tryCatch(
    form_blocks_cpp(X, as.integer(s)),
    error = function(e) {
      X <- as.matrix(X); s <- as.integer(s)
      if (s <= 0L) return(list(X))
      Tn <- nrow(X); m <- ncol(X); N <- Tn - s
      if (N <= 0L) stop(".form_blocks: not enough rows for lag order")
      lapply(1:(s + 1L), function(i) X[i:(i + N - 1L), , drop = FALSE])
    }
  )
}

#' @keywords internal
.make_t_lags <- function(t, s) {
  Tn <- length(t)
  stopifnot(Tn > s)
  N <- Tn - s
  ts1 <- t[(s + 1L):Tn]
  Ts <- matrix(NA_real_, nrow = N, ncol = s)
  for (i in 1:s) {
    Ts[, i] <- t[(s + 1L - i):(Tn - i)]
  }
  list(ts1 = ts1, Ts = Ts)
}

#' @keywords internal
.solve_normal <- function(A, b, ridge = 0) {
  if (ridge > 0) {
    A <- A + diag(ridge, nrow(A))
  }
  out <- tryCatch(
    solve(A, b),
    error = function(e) {
      sv <- svd(A)
      tol <- max(dim(A)) * .Machine$double.eps * max(sv$d)
      di <- ifelse(sv$d > tol, 1 / sv$d, 0)
      sv$v %*% (di * t(sv$u) %*% b)
    }
  )
  if (is.matrix(b)) {
    out
  } else {
    as.numeric(out)
  }
}

#' @keywords internal
.orth_resid <- function(X, t) {
  p <- as.numeric(crossprod(X, t) / sum(t^2))
  X - tcrossprod(t, p)
}

#' @keywords internal
.as_arma_meta <- function(p = 0L, d = 0L, q = 0L, P = 0L, Q = 0L, m = 1L, D = 0L) {
  c(as.integer(p), as.integer(q), as.integer(P), as.integer(Q), as.integer(m),
    as.integer(d), as.integer(D))
}
