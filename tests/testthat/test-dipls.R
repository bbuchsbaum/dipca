
test_that("dipls produces coherent projector objects", {
  set.seed(123)
  Tn <- 120; m <- 5; p <- 3; s <- 2
  X <- matrix(rnorm(Tn * m), Tn, m)
  Y <- matrix(rnorm(Tn * p), Tn, p)

  fit <- suppressWarnings(
    dipls(
      X, Y,
      n_comp = 2, s = s,
      preproc_x = multivarious::pass(),
      preproc_y = multivarious::pass(),
      mode = "fir", verbose = 0
    )
  )

  expect_s3_class(fit, c("dipls", "cross_projector"))
  expect_equal(dim(fit$vx), c(m, 2))
  expect_equal(dim(fit$vy), c(p, 2))
  expect_equal(fit$lag_order, s)
  expect_equal(length(fit$inner), 2)

  scores <- suppressWarnings(dipls_scores(fit, X))
  expect_equal(dim(scores), c(Tn, 2))
  expect_lt(max(abs(scores - fit$T), na.rm = TRUE), 1e-6)

  preds <- suppressWarnings(predict(fit, X))
  expect_equal(dim(preds), c(Tn, p))
  expect_true(all(is.na(preds[seq_len(s), ])))

  res <- suppressWarnings(residuals(fit, list(X = X, Y = Y)))
  expect_equal(dim(res), c(Tn, p))
  expect_true(all(is.na(res[seq_len(s), ])))
})

test_that("dipls (FIR) recovers dynamic mapping on synthetic data", {
  skip_on_cran()
  set.seed(42)

  Tn <- 250; m <- 5; p <- 2; s <- 3
  X <- matrix(rnorm(Tn * m), Tn, m)

  w_true <- rnorm(m); w_true <- w_true / sqrt(sum(w_true^2))
  t_true <- as.numeric(X %*% w_true)
  beta_true <- c(0.7, -0.2, 0.15, 0.05)
  Ts <- dipca:::dipls_build_Ts_cpp(t_true, s)
  idx <- seq.int(s + 1L, Tn)

  u <- numeric(Tn)
  u[idx] <- as.numeric(Ts %*% beta_true)

  q_true <- rnorm(p); q_true <- q_true / sqrt(sum(q_true^2))
  Y <- u %o% q_true
  Y[idx, ] <- Y[idx, ] + matrix(rnorm(length(idx) * p, sd = 0.01), length(idx), p)
  Y[seq_len(s), ] <- Y[seq_len(s), ] + matrix(rnorm(s * p, sd = 0.01), s, p)

  fit <- suppressWarnings(
    dipls(
      X, Y,
      n_comp = 1, s = s,
      preproc_x = multivarious::pass(),
      preproc_y = multivarious::pass(),
      mode = "fir", verbose = 0
    )
  )

  expect_gt(abs(cor(fit$vx[, 1], w_true)), 0.95)
  expect_gt(abs(cor(fit$vy[, 1], q_true)), 0.95)
  expect_gt(abs(cor(fit$inner[[1]]$beta, beta_true)), 0.98)

  Yhat <- suppressWarnings(predict(fit, X))
  mspe_val <- mean((Yhat[idx, ] - Y[idx, ])^2)
  expect_lt(mspe_val, 0.01)
})

test_that("dipls (ARX) captures latent recursion", {
  skip_on_cran()
  set.seed(99)

  Tn <- 250; m <- 6; p <- 3; s <- 2
  X <- matrix(rnorm(Tn * m), Tn, m)

  w_true <- rnorm(m); w_true <- w_true / sqrt(sum(w_true^2))
  t_true <- as.numeric(X %*% w_true)
  Ts <- dipca:::dipls_build_Ts_cpp(t_true, s)
  idx <- seq.int(s + 1L, Tn)

  beta_true <- c(0.6, -0.25, 0.1)
  alpha_true <- c(0.4, -0.15)
  u <- numeric(Tn)
  for (k in idx) {
    tb <- sum(beta_true * Ts[k - s, ])
    ub <- sum(alpha_true * rev(u[(k - s):(k - 1)]))
    u[k] <- tb + ub + rnorm(1, sd = 0.01)
  }

  q_true <- rnorm(p); q_true <- q_true / sqrt(sum(q_true^2))
  Y <- u %o% q_true + matrix(rnorm(Tn * p, sd = 0.01), Tn, p)

  fit <- suppressWarnings(
    dipls(
      X, Y,
      n_comp = 1, s = s,
      preproc_x = multivarious::pass(),
      preproc_y = multivarious::pass(),
      mode = "arx", verbose = 0
    )
  )

  expect_gt(abs(cor(fit$inner[[1]]$alpha, alpha_true)), 0.9)
  expect_gt(abs(cor(fit$inner[[1]]$beta, beta_true)), 0.95)

  Yhat <- suppressWarnings(predict(fit, X))
  mspe_val <- mean((Yhat[idx, ] - Y[idx, ])^2)
  expect_lt(mspe_val, 0.2)
})

test_that("dipls recovers multiple components (FIR)", {
  skip_on_cran()
  set.seed(1235)

  Tn <- 220; m <- 7; p <- 4; s <- 2; h <- 2
  X <- matrix(rnorm(Tn * m), Tn, m)

  W_true <- qr.Q(qr(matrix(rnorm(m * h), m, h)))
  T_latent <- X %*% W_true

  beta_true <- rbind(
    c(0.6, -0.2, 0.1),
    c(0.5, 0.25, -0.15)
  )

  U_latent <- matrix(0, Tn, h)
  idx <- seq.int(s + 1L, Tn)
  for (j in seq_len(h)) {
    Ts <- dipca:::dipls_build_Ts_cpp(T_latent[, j], s)
    U_latent[idx, j] <- as.numeric(Ts %*% beta_true[j, ])
  }

  Q_true <- qr.Q(qr(matrix(rnorm(p * h), p, h)))
  Y <- matrix(0, Tn, p)
  for (j in seq_len(h)) {
    Y <- Y + U_latent[, j] %o% Q_true[, j]
  }
  Y <- Y + matrix(rnorm(Tn * p, sd = 0.02), Tn, p)

  all_perms <- function(idx) {
    if (length(idx) == 1L) return(matrix(idx, nrow = 1L))
    out <- NULL
    for (i in seq_along(idx)) {
      rest <- idx[-i]
      sub <- all_perms(rest)
      out <- rbind(out, cbind(idx[i], sub))
    }
    out
  }

  fit <- suppressWarnings(
    dipls(
      X, Y,
      n_comp = h, s = s,
      preproc_x = multivarious::pass(),
      preproc_y = multivarious::pass(),
      mode = "fir", verbose = 0
    )
  )

  perm_mat <- all_perms(seq_len(h))
  vx_cor <- abs(cor(fit$vx, W_true))
  scores <- apply(perm_mat, 1, function(p) sum(vx_cor[cbind(seq_len(h), p)]))
  best_perm <- perm_mat[which.max(scores), ]
  if (h == 1) best_perm <- 1L

  vx_aligned <- fit$vx[, best_perm, drop = FALSE]
  vy_aligned <- fit$vy[, best_perm, drop = FALSE]
  inner_aligned <- fit$inner[best_perm]

  vx_sv <- svd(t(vx_aligned) %*% W_true)$d
  vy_sv <- svd(t(vy_aligned) %*% Q_true)$d
  expect_gt(min(vx_sv), 0.9)
  expect_gt(min(vy_sv), 0.9)

  for (j in seq_len(h)) {
    expect_gt(abs(cor(inner_aligned[[j]]$beta, beta_true[j, ])), 0.9)
  }

  Yhat <- suppressWarnings(predict(fit, X))
  mspe_val <- mean((Yhat[idx, ] - Y[idx, ])^2)
  expect_lt(mspe_val, 0.02)
})

test_that("dipls handles differing preprocessing pipelines", {
  set.seed(321)
  Tn <- 80; m <- 4; p <- 3; s <- 1
  X <- matrix(rnorm(Tn * m), Tn, m)
  Y <- matrix(rnorm(Tn * p), Tn, p)
  Y <- sweep(Y, 2, seq_len(p) * 5, "+")

  fit <- suppressWarnings(
    dipls(
      X, Y,
      n_comp = 1, s = s,
      preproc_x = multivarious::pass(),
      preproc_y = multivarious::center(),
      mode = "fir", verbose = 0
    )
  )

  Yhat <- suppressWarnings(predict(fit, X))
  resid <- suppressWarnings(residuals(fit, list(X = X, Y = Y)))

  expect_equal(dim(Yhat), c(Tn, p))
  expect_equal(dim(resid), c(Tn, p))
  expect_true(all(is.na(Yhat[seq_len(s), ])))
  expect_equal(resid, Y - Yhat)

  new_X <- X[seq_len(10), , drop = FALSE]
  scores_new <- suppressWarnings(dipls_scores(fit, new_X))
  expect_equal(nrow(scores_new), nrow(new_X))
  expect_equal(ncol(scores_new), 1)
})
