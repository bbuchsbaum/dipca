# DiCCA Diagnostic Tests
# Tests canonical correlation objective: descending predictability (R²)
# Based on Dong & Qin (2018) DiCCA algorithm and monotone R² property

context("DiCCA diagnostics (AR inner)")

test_that("DiCCA (AR): descending predictability and whitening", {
  skip_on_cran()

  set.seed(123)  # For reproducibility
  sim <- sim_VAR(N = 240, m = 8, h = 4, A_diag = c(0.9, 0.7, 0.5, 0.2))

  # Fit DiCCA with AR(4) inner model using multivarious preprocessing
  # Use n_init=10 for better convergence to global optimum
  fit <- dicca(sim$X, s = 4, l = 4, inner = "classic", n_init = 10,
              preproc = multivarious::center(),
              max_iter = 500, tol = 1e-7, seed = 456)

  # Access fields via bi_projector interface
  T_scores <- multivarious::scores(fit)
  P_load   <- fit$loadings
  betas_matrix <- fit$betas
  s <- fit$lag_order

  # (a) Monotone nonincreasing R² across components
  # DiCCA canonical correlation objective ensures this property
  R2s <- fit$R2

  # Check all are valid R² values
  expect_true(all(R2s >= 0 & R2s <= 1))

  # Check monotone nonincreasing
  expect_true(all(diff(R2s) <= 1e-6),
             label = "R² values should be nonincreasing across components")

  # First component should have reasonable predictability
  # Note: stochastic optimization may not always reach global optimum
  expect_gt(R2s[1], 0.5,
           label = "First component should be reasonably predictable")

  # (b) Residual whitening vs raw data
  # DiCCA should reduce dynamic correlations
  t_hat <- one_step_predict_scores_AR(T_scores, betas_matrix, s)
  Xhat <- t_hat %*% t(P_load)

  # Compute residuals using multivarious preprocessing
  Xs <- multivarious::transform(fit$preproc, sim$X)
  E_preprocessed <- Xs - Xhat
  E <- multivarious::inverse_transform(fit$preproc, E_preprocessed)

  raw_acf <- acf_energy(sim$X, lag.max = 10)
  res_acf <- acf_energy(E, lag.max = 10)

  # Residuals should have reduced autocorrelation (relaxed threshold for stochastic optimization)
  expect_lt(res_acf, 0.75 * raw_acf,
           label = "Residuals should have less lag-correlation than raw data")
})

test_that("DiCCA (AR): R² computation matches direct calculation", {
  skip_on_cran()

  sim <- sim_VAR(N = 200, m = 6, h = 3, A_diag = c(0.85, 0.6, 0.3))

  fit <- dicca(sim$X, s = 3, l = 3, inner = "classic", n_init = 2,
              preproc = multivarious::center(),
              max_iter = 400)

  T_scores <- multivarious::scores(fit)
  betas_matrix <- fit$betas
  s <- fit$lag_order

  # Compute one-step predictions
  t_hat <- one_step_predict_scores_AR(T_scores, betas_matrix, s)

  # Compute R² directly for each component
  R2_direct <- vapply(1:3, function(j) {
    valid_idx <- (s + 1):nrow(T_scores)
    cor(T_scores[valid_idx, j], t_hat[valid_idx, j])^2
  }, numeric(1))

  # Compare with stored R² values
  expect_equal(fit$R2, R2_direct, tolerance = 1e-6,
              label = "Stored R² should match direct calculation")
})

test_that("DiCCA (AR): prediction and residuals methods", {
  skip_on_cran()

  sim <- sim_VAR(N = 150, m = 5, h = 2, A_diag = c(0.8, 0.5))
  X <- scale(sim$X, center = TRUE, scale = FALSE)

  fit <- dicca(X, s = 2, l = 2, inner = "classic", n_init = 2)

  # Test predict() S3 method
  pred <- predict(fit, X)
  expect_equal(names(pred), c("scores", "scores_hat"))
  expect_equal(dim(pred$scores), c(150, 2))
  expect_equal(dim(pred$scores_hat), c(150, 2))

  # First s rows should be NA
  expect_true(all(is.na(pred$scores_hat[1:2, ])))

  # Test residuals() S3 method
  resid <- residuals(fit, X)
  expect_equal(names(resid), c("v", "e_hat", "scores"))

  # Residuals should be finite (except first s time points)
  expect_true(all(is.finite(resid$v[(fit$lag_order + 1):150, ])))
  expect_true(all(is.finite(resid$e_hat[(fit$lag_order + 1):150, ])))
})

test_that("DiCCA (AR): handles different lag orders", {
  skip_on_cran()

  sim <- sim_VAR(N = 180, m = 6, h = 2, A_diag = c(0.7, 0.5))
  X <- scale(sim$X, center = TRUE, scale = FALSE)

  # Fit with different lag orders
  fit_s1 <- dicca(X, s = 1, l = 2, inner = "classic", n_init = 2, max_iter = 300)
  fit_s4 <- dicca(X, s = 4, l = 2, inner = "classic", n_init = 2, max_iter = 300)

  # Both should have valid R² values
  expect_true(all(fit_s1$R2 >= 0 & fit_s1$R2 <= 1))
  expect_true(all(fit_s4$R2 >= 0 & fit_s4$R2 <= 1))

  # Higher lag order might capture more dynamics (higher R²)
  # But this is not guaranteed, so just check they're reasonable
  expect_gt(fit_s4$R2[1], 0.3)

  # Betas should have correct dimensions
  expect_equal(dim(fit_s1$betas), c(2, 1))
  expect_equal(dim(fit_s4$betas), c(2, 4))
})

test_that("DiCCA (AR): convergence properties", {
  skip_on_cran()

  sim <- sim_VAR(N = 200, m = 6, h = 2, A_diag = c(0.75, 0.5))

  # Use multivarious preprocessing
  fit <- dicca(sim$X, s = 2, l = 2, inner = "classic", n_init = 2,
              preproc = multivarious::center(),
              max_iter = 500, tol = 1e-7)

  # Check convergence
  expect_true(any(fit$iters_per_component < 500),
             label = "Algorithm should converge before max_iter")

  # Objective history should exist
  expect_true(length(fit$obj_history) > 0)

  # Parse and check objective values
  for (obj_str in fit$obj_history) {
    obj_val <- as.numeric(strsplit(obj_str, ":")[[1]][2])
    expect_true(is.finite(obj_val))
    expect_gt(obj_val, 0)
  }
})
