# DiPCA Diagnostic Tests
# Tests core theoretical properties: orthogonality, whitening, predictability
# Based on Dong & Qin (2018) DiPCA objective and whitening-filter view

context("DiPCA diagnostics")

test_that("DiPCA: orthogonality and whitening on small VAR", {
  skip_on_cran()

  set.seed(123)  # For reproducibility
  sim <- sim_VAR(N = 240, m = 8, h = 3, A_diag = c(0.85, 0.6, 0.3))

  # Fit DiPCA with AR(3) inner model using multivarious preprocessing
  # Use n_init=10 for better convergence
  fit <- dipca(sim$X, s = 3, l = 3, n_init = 10,
              preproc = multivarious::center(),
              max_iter = 500, tol = 1e-7, seed = 456)

  # Field mapping (multivarious bi_projector structure)
  T_scores <- multivarious::scores(fit)  # N x h
  P_load   <- fit$loadings               # m x h
  betas_matrix <- fit$betas              # l x s matrix (NOT list)
  s <- fit$lag_order

  expect_equal(ncol(T_scores), 3)
  expect_equal(ncol(P_load), 3)

  # 1) Scores nearly orthogonal
  # DiPCA extracts orthogonal latent components via deflation
  C <- cor(T_scores)
  diag(C) <- 0
  expect_lt(max(abs(C)), 0.1,
           label = "Latent scores should be nearly orthogonal")

  # 2) One-step prediction error has less autocorrelation than raw X
  # This tests the dynamic whitening filter property
  t_hat <- one_step_predict_scores_AR(T_scores, betas_matrix, s)
  Xhat  <- t_hat %*% t(P_load)

  # Compute residuals in original scale
  X_preprocessed <- multivarious::transform(fit$preproc, sim$X)
  E_preprocessed <- X_preprocessed - Xhat
  E <- multivarious::inverse_transform(fit$preproc, E_preprocessed)

  raw_acf <- acf_energy(sim$X, lag.max = 10)
  res_acf <- acf_energy(E, lag.max = 10)

  # Relaxed threshold for stochastic optimization
  expect_lt(res_acf, 0.75 * raw_acf,
           label = "Residuals should have less lag-correlation than raw data")

  # 3) First component should have reasonable predictability
  # DiPCA extracts components in order of predictability
  # Note: stochastic optimization may not always reach global optimum
  t1    <- T_scores[, 1]
  t1hat <- t_hat[, 1]
  R2_1  <- cor(t1, t1hat)^2

  expect_gt(R2_1, 0.5,
           label = "First component should be reasonably predictable")
})

test_that("DiPCA: convergence and objective values", {
  skip_on_cran()

  sim <- sim_VAR(N = 200, m = 6, h = 2, A_diag = c(0.8, 0.5))
  X <- scale(sim$X, center = TRUE, scale = FALSE)

  fit <- dipca(X, s = 2, l = 2, n_init = 2, max_iter = 300, tol = 1e-7)

  # Check that algorithm converged (iters < max_iter for at least one component)
  expect_true(any(fit$iters_per_component < 300),
             label = "Algorithm should converge before max_iter")

  # Objective history should be parseable
  expect_true(length(fit$obj_history) > 0)

  # Parse last objective value
  last_obj_str <- tail(fit$obj_history, 1)
  obj_val <- as.numeric(strsplit(last_obj_str, ":")[[1]][2])
  expect_true(is.finite(obj_val))
  expect_gt(obj_val, 0)
})

test_that("DiPCA: Algorithm I vs II produce similar results", {
  skip_on_cran()

  sim <- sim_VAR(N = 180, m = 6, h = 2, A_diag = c(0.75, 0.5))
  X <- scale(sim$X, center = TRUE, scale = FALSE)

  # Algorithm I: single power step
  f1 <- dipca(X, s = 2, l = 2, algorithm = "I", n_init = 2,
             max_iter = 400, tol = 1e-7, seed = 123)

  # Algorithm II: inner power iterations
  f2 <- dipca(X, s = 2, l = 2, algorithm = "II", n_init = 2,
             max_iter = 400, tol = 1e-7, inner_power = 100,
             inner_tol = 1e-9, seed = 123)

  # Extract final objectives
  obj1 <- as.numeric(strsplit(tail(f1$obj_history, 1), ":")[[1]][2])
  obj2 <- as.numeric(strsplit(tail(f2$obj_history, 1), ":")[[1]][2])

  # Should be close (within 5%)
  expect_lt(abs(obj1 - obj2) / max(1, abs(obj2)), 0.05,
           label = "Algorithm I and II should reach similar objectives")
})

test_that("DiPCA: prediction and residuals methods work correctly", {
  skip_on_cran()

  sim <- sim_VAR(N = 150, m = 5, h = 2, A_diag = c(0.7, 0.4))
  X <- scale(sim$X, center = TRUE, scale = FALSE)

  fit <- dipca(X, s = 2, l = 2, n_init = 2, max_iter = 300)

  # Test predict() S3 method
  pred <- predict(fit, X)
  expect_equal(names(pred), c("scores", "scores_hat"))
  expect_equal(dim(pred$scores), c(150, 2))
  expect_equal(dim(pred$scores_hat), c(150, 2))

  # First s rows of scores_hat should be NA
  expect_true(all(is.na(pred$scores_hat[1:2, ])))
  expect_false(any(is.na(pred$scores_hat[3:150, ])))

  # Test residuals() S3 method
  resid <- residuals(fit, X)
  expect_equal(names(resid), c("v", "e_hat", "scores"))
  expect_equal(dim(resid$v), c(150, 2))
  expect_equal(dim(resid$e_hat), c(150, 5))
  expect_equal(dim(resid$scores), c(150, 2))

  # Score residuals v = scores - scores_hat
  expect_equal(resid$v, pred$scores - pred$scores_hat)
})
