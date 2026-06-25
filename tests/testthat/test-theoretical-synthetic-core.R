score_prediction_mse <- function(fit, Xnew, component = 1L) {
  pred <- predict(fit, Xnew)
  s <- fit$lag_order
  idx <- seq.int(s + 1L, nrow(pred$scores))
  mean((pred$scores_hat[idx, component] - pred$scores[idx, component])^2)
}

simulate_exact_ar3_observed <- function(seed = 320L,
                                        N = 240L,
                                        phi = c(0.8, -0.55, 0.35),
                                        noise_sd = 0.15,
                                        n_noise = 3L) {
  set.seed(seed)
  t_latent <- numeric(N)
  eps <- rnorm(N, sd = noise_sd)
  for (k in 4:N) {
    t_latent[k] <- sum(phi * rev(t_latent[(k - 3L):(k - 1L)])) + eps[k]
  }

  X <- cbind(t_latent, matrix(rnorm(N * n_noise, sd = 0.01), N, n_noise))
  list(X = X, t = t_latent, phi = phi)
}

test_that("DiPCA and DiCCA recover latent factors in persistence order", {
  skip_on_cran()

  set.seed(101)
  sim <- sim_VAR(
    N = 320, m = 8, h = 4,
    A_diag = c(0.95, 0.7, 0.35, 0.0),
    sigma_v = 0.2, sigma_e = 0.05
  )
  true_r2 <- c(0.95, 0.7, 0.35, 0.0)^2

  check_persistence_fit <- function(fit, min_diag_cor) {
    C <- abs(cor(multivarious::scores(fit), sim$T))
    diag_cor <- diag(C)

    expect_equal(max.col(C[1:3, , drop = FALSE], ties.method = "first"), 1:3)
    expect_true(all(diag_cor[1:3] > min_diag_cor))
    expect_equal(fit$R2[1:3], true_r2[1:3], tolerance = 0.1)
    expect_lt(fit$R2[4], 0.05)
  }

  fit_dipca <- dipca(sim$X, s = 1, l = 4, n_init = 4, max_iter = 500, seed = 1)
  fit_dicca <- dicca(sim$X, s = 1, l = 4, n_init = 4, max_iter = 500, seed = 1, inner = "classic")

  check_persistence_fit(fit_dipca, c(0.99, 0.95, 0.95))
  check_persistence_fit(fit_dicca, c(0.98, 0.95, 0.94))
})

test_that("DiPCA prefers the true AR(3) lag on an exact latent process", {
  skip_on_cran()

  sim <- simulate_exact_ar3_observed()
  train_idx <- 1:160
  test_idx <- 161:nrow(sim$X)
  Xtr <- sim$X[train_idx, , drop = FALSE]
  Xte <- sim$X[test_idx, , drop = FALSE]

  fit_s1 <- dipca(Xtr, s = 1, l = 1, n_init = 2, max_iter = 300, seed = 7)
  fit_s3 <- dipca(Xtr, s = 3, l = 1, n_init = 2, max_iter = 300, seed = 7)
  fit_s5 <- dipca(Xtr, s = 5, l = 1, n_init = 2, max_iter = 300, seed = 7)

  mse_s1 <- score_prediction_mse(fit_s1, Xte)
  mse_s3 <- score_prediction_mse(fit_s3, Xte)
  mse_s5 <- score_prediction_mse(fit_s5, Xte)
  scores_s3 <- predict(fit_s3, Xte)$scores[, 1]

  expect_gt(abs(cor(scores_s3, sim$t[test_idx])), 0.99)
  expect_lt(mse_s3, 0.85 * mse_s1)
  expect_lt(mse_s3, mse_s5)
})

test_that("DiPLS uses the forward temporal direction for purely lagged coupling", {
  skip_on_cran()

  set.seed(404)
  Tn <- 220
  m <- 6
  p <- 2
  s <- 2
  X <- matrix(rnorm(Tn * m), Tn, m)

  w_true <- rnorm(m)
  w_true <- w_true / sqrt(sum(w_true^2))
  t_true <- as.numeric(X %*% w_true)
  Ts <- dipca:::dipls_build_Ts_cpp(t_true, s)
  beta_true <- c(0.0, 0.9, -0.45)
  idx <- seq.int(s + 1L, Tn)

  u <- numeric(Tn)
  u[idx] <- as.numeric(Ts %*% beta_true)

  q_true <- rnorm(p)
  q_true <- q_true / sqrt(sum(q_true^2))
  Y <- u %o% q_true + matrix(rnorm(Tn * p, sd = 0.02), Tn, p)

  fit_forward <- dipls(
    X, Y,
    n_comp = 1, s = s,
    preproc_x = multivarious::pass(),
    preproc_y = multivarious::pass(),
    mode = "fir", verbose = 0
  )
  fit_reversed <- dipls(
    X[Tn:1, , drop = FALSE], Y,
    n_comp = 1, s = s,
    preproc_x = multivarious::pass(),
    preproc_y = multivarious::pass(),
    mode = "fir", verbose = 0
  )

  Yhat_forward <- predict(fit_forward, X)
  Yhat_reversed <- predict(fit_reversed, X[Tn:1, , drop = FALSE])
  mse_forward <- mean((Yhat_forward[idx, ] - Y[idx, ])^2)
  mse_reversed <- mean((Yhat_reversed[idx, ] - Y[idx, ])^2)

  expect_gt(abs(cor(fit_forward$inner[[1]]$beta, beta_true)), 0.99)
  expect_lt(mse_forward, 0.05)
  expect_gt(mse_reversed, 10 * mse_forward)
})
