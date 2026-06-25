align_signed_columns <- function(reference, candidate, start_row = 1L) {
  reference <- as.matrix(reference)
  candidate <- as.matrix(candidate)
  out <- candidate
  rows <- seq.int(start_row, nrow(reference))

  for (j in seq_len(ncol(reference))) {
    keep <- is.finite(reference[rows, j]) & is.finite(candidate[rows, j])
    if (!any(keep)) next
    rho <- suppressWarnings(stats::cor(reference[rows, j][keep], candidate[rows, j][keep]))
    if (is.finite(rho) && rho < 0) {
      out[, j] <- -out[, j]
    }
  }

  out
}

constant_shift_matrix <- function(n, shift) {
  matrix(rep(shift, each = n), nrow = n)
}

test_that("DiPCA is invariant to additive column offsets under centering", {
  set.seed(11)
  sim <- sim_VAR(N = 120, m = 6, h = 2, A_diag = c(0.8, 0.5))
  X_shift <- sweep(sim$X, 2, c(3, -2, 1, 4, -1, 2), "+")

  fit_ref <- dipca(sim$X, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 99)
  fit_shift <- dipca(X_shift, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 99)

  scores_ref <- multivarious::scores(fit_ref)
  scores_shift <- align_signed_columns(scores_ref, multivarious::scores(fit_shift))
  pred_ref <- predict(fit_ref, sim$X)$scores_hat
  pred_shift <- align_signed_columns(pred_ref, predict(fit_shift, X_shift)$scores_hat, start_row = 3L)
  resid_ref <- residuals(fit_ref, sim$X)$e_hat
  resid_shift <- residuals(fit_shift, X_shift)$e_hat

  expect_equal(scores_shift, scores_ref, tolerance = 1e-10)
  expect_equal(pred_shift, pred_ref, tolerance = 1e-10)
  expect_equal(resid_shift, resid_ref, tolerance = 1e-10)
})

test_that("DiCCA is invariant to additive column offsets under centering", {
  set.seed(11)
  sim <- sim_VAR(N = 120, m = 6, h = 2, A_diag = c(0.8, 0.5))
  X_shift <- sweep(sim$X, 2, c(3, -2, 1, 4, -1, 2), "+")

  fit_ref <- dicca(sim$X, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 99, inner = "classic")
  fit_shift <- dicca(X_shift, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 99, inner = "classic")

  scores_ref <- multivarious::scores(fit_ref)
  scores_shift <- align_signed_columns(scores_ref, multivarious::scores(fit_shift))
  pred_ref <- predict(fit_ref, sim$X)$scores_hat
  pred_shift <- align_signed_columns(pred_ref, predict(fit_shift, X_shift)$scores_hat, start_row = 3L)
  resid_ref <- residuals(fit_ref, sim$X)$e_hat
  resid_shift <- residuals(fit_shift, X_shift)$e_hat

  expect_equal(scores_shift, scores_ref, tolerance = 1e-10)
  expect_equal(pred_shift, pred_ref, tolerance = 1e-6)
  expect_equal(resid_shift, resid_ref, tolerance = 1e-6)
  expect_equal(fit_shift$R2, fit_ref$R2, tolerance = 1e-10)
})

test_that("DiPCA remains finite on rank-deficient inputs", {
  set.seed(4)
  X_base <- matrix(rnorm(150 * 4), 150, 4)
  X <- cbind(
    X_base[, 1],
    X_base[, 1],
    X_base[, 2] + X_base[, 3],
    X_base[, 2] + X_base[, 3],
    1
  )

  fit <- dipca(X, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 10)
  pred <- predict(fit, X)
  resid <- residuals(fit, X)

  expect_true(all(is.finite(fit$v)))
  expect_true(all(is.finite(fit$loadings)))
  expect_true(all(is.finite(fit$theta)))
  expect_true(all(is.finite(pred$scores[-(1:2), ])))
  expect_true(all(is.finite(pred$scores_hat[-(1:2), ])))
  expect_true(all(is.finite(resid$e_hat[-(1:2), ])))
})

test_that("DiCCA remains finite on rank-deficient inputs", {
  set.seed(4)
  X_base <- matrix(rnorm(150 * 4), 150, 4)
  X <- cbind(
    X_base[, 1],
    X_base[, 1],
    X_base[, 2] + X_base[, 3],
    X_base[, 2] + X_base[, 3],
    1
  )

  fit <- dicca(X, s = 2, l = 2, n_init = 2, max_iter = 300, seed = 10, inner = "classic")
  pred <- predict(fit, X)
  resid <- residuals(fit, X)

  expect_true(all(is.finite(fit$v)))
  expect_true(all(is.finite(fit$loadings)))
  expect_true(all(is.finite(fit$theta)))
  expect_true(all(is.finite(fit$betas)))
  expect_true(all(is.finite(pred$scores[-(1:2), ])))
  expect_true(all(is.finite(pred$scores_hat[-(1:2), ])))
  expect_true(all(is.finite(resid$e_hat[-(1:2), ])))
})

test_that("DiPLS FIR mode is shift-equivariant under centered preprocessing", {
  set.seed(7)
  Tn <- 100
  m <- 5
  p <- 3
  s <- 2
  X <- matrix(rnorm(Tn * m), Tn, m)
  Y <- matrix(rnorm(Tn * p), Tn, p)
  x_shift <- c(10, -3, 5, 1, -4)
  y_shift <- c(2, -6, 4)
  X_shift <- sweep(X, 2, x_shift, "+")
  Y_shift <- sweep(Y, 2, y_shift, "+")

  fit_ref <- dipls(X, Y, n_comp = 1, s = s, mode = "fir", verbose = 0)
  fit_shift <- dipls(X_shift, Y_shift, n_comp = 1, s = s, mode = "fir", verbose = 0)

  scores_ref <- dipls_scores(fit_ref, X)
  scores_shift <- align_signed_columns(scores_ref, dipls_scores(fit_shift, X_shift))
  pred_ref <- predict(fit_ref, X)
  pred_shift <- predict(fit_shift, X_shift)
  resid_ref <- residuals(fit_ref, list(X = X, Y = Y))
  resid_shift <- residuals(fit_shift, list(X = X_shift, Y = Y_shift))

  expect_equal(scores_shift, scores_ref, tolerance = 1e-10)
  expect_equal(pred_shift, pred_ref + constant_shift_matrix(Tn, y_shift), tolerance = 1e-10)
  expect_equal(resid_shift, resid_ref, tolerance = 1e-10)
})

test_that("DiPLS supports s = 0 on rank-deficient centered blocks", {
  set.seed(4)
  Tn <- 120
  X_base <- matrix(rnorm(Tn * 3), Tn, 3)
  Y_base <- matrix(rnorm(Tn * 2), Tn, 2)
  X <- cbind(X_base[, 1], X_base[, 1], X_base[, 2] + X_base[, 3], 1)
  Y <- cbind(Y_base[, 1], Y_base[, 1], 2)

  fit <- dipls(
    X, Y,
    n_comp = 1, s = 0,
    preproc_x = multivarious::center(),
    preproc_y = multivarious::center(),
    mode = "fir", verbose = 0
  )

  scores <- dipls_scores(fit, X)
  pred <- predict(fit, X)
  resid <- residuals(fit, list(X = X, Y = Y))

  expect_equal(scores, fit$T, tolerance = 1e-6)
  expect_length(fit$inner[[1]]$beta, 1)
  expect_true(all(is.finite(fit$vx)))
  expect_true(all(is.finite(fit$vy)))
  expect_true(all(is.finite(fit$R)))
  expect_true(all(is.finite(pred)))
  expect_true(all(is.finite(resid)))
  expect_false(any(is.na(pred)))
  expect_false(any(is.na(resid)))
})
