local_align_sign <- function(reference, candidate) {
  if (suppressWarnings(stats::cor(reference, candidate)) < 0) {
    -candidate
  } else {
    candidate
  }
}

local_exact_ar_observation <- function(seed = 1L, n = 180L, phi = 0.82,
                                       noise_sd = 0.12, obs_noise = 0.015) {
  set.seed(seed)
  t <- numeric(n)
  eps <- rnorm(n, sd = noise_sd)
  for (i in 2:n) {
    t[i] <- phi * t[i - 1L] + eps[i]
  }

  p <- c(0.7, -0.4, 0.2, 0.5, -0.1)
  p <- p / sqrt(sum(p^2))
  X <- t %o% p + matrix(rnorm(n * length(p), sd = obs_noise), n, length(p))
  list(X = X, t = t, p = p, phi = phi)
}

local_ar_mse <- function(scores, scores_hat, s) {
  idx <- seq.int(s + 1L, nrow(scores))
  mean((scores[idx, 1L] - scores_hat[idx, 1L])^2)
}

test_that("DiPCA recovers a one-dimensional AR latent process and its predictive law", {
  sim <- local_exact_ar_observation(seed = 10)
  fit <- dipca(
    sim$X, s = 1, l = 1,
    n_init = 3, max_iter = 250, seed = 99,
    preproc = multivarious::center()
  )

  scores <- local_align_sign(sim$t, drop(multivarious::scores(fit)))
  loading <- local_align_sign(sim$p, drop(fit$loadings))
  pred <- predict(fit, sim$X)
  ar_coef <- drop(.lm.fit(cbind(scores[-length(scores)]), scores[-1L])$coefficients)
  null_mse <- mean((scores[-1L] - mean(scores[-1L]))^2)
  model_mse <- local_ar_mse(pred$scores, pred$scores_hat, fit$lag_order)

  expect_gt(abs(stats::cor(scores, sim$t)), 0.985)
  expect_gt(abs(stats::cor(loading, sim$p)), 0.985)
  expect_equal(ar_coef, sim$phi, tolerance = 0.08)
  expect_gt(fit$R2[1L], 0.55)
  expect_lt(model_mse, 0.55 * null_mse)
})

test_that("Classic DiCCA recovers a one-dimensional AR latent process and stores direct R2", {
  sim <- local_exact_ar_observation(seed = 11, phi = 0.76)
  fit <- dicca(
    sim$X, s = 1, l = 1, inner = "classic",
    n_init = 3, max_iter = 250, seed = 99,
    preproc = multivarious::center()
  )

  pred <- predict(fit, sim$X)
  idx <- seq.int(fit$lag_order + 1L, nrow(pred$scores))
  direct_r2 <- stats::cor(pred$scores[idx, 1L], pred$scores_hat[idx, 1L])^2
  scores <- local_align_sign(sim$t, drop(multivarious::scores(fit)))
  ar_coef <- drop(.lm.fit(cbind(scores[-length(scores)]), scores[-1L])$coefficients)

  expect_gt(abs(stats::cor(scores, sim$t)), 0.985)
  expect_equal(ar_coef, sim$phi, tolerance = 0.08)
  expect_equal(fit$R2[1L], direct_r2, tolerance = 1e-8)
  expect_gt(fit$R2[1L], 0.45)
})

test_that("DiPCA deflation extracts components in descending dynamic predictability", {
  set.seed(12)
  sim <- sim_VAR(
    N = 220, m = 7, h = 3,
    A_diag = c(0.88, 0.55, 0.18),
    sigma_v = 0.18, sigma_e = 0.04
  )

  fit <- dipca(
    sim$X, s = 1, l = 3,
    n_init = 4, max_iter = 300, seed = 12,
    preproc = multivarious::center()
  )
  C <- abs(stats::cor(multivarious::scores(fit), sim$T))

  expect_equal(max.col(C, ties.method = "first"), 1:3)
  expect_true(all(diag(C) > c(0.96, 0.9, 0.7)))
  expect_true(all(diff(fit$R2) <= 1e-8))
  expect_equal(fit$R2[1:2], c(0.88, 0.55)^2, tolerance = 0.18)
})

test_that("DiPCA Algorithm II is active and reaches a comparable predictive optimum", {
  set.seed(13)
  sim <- sim_VAR(N = 180, m = 6, h = 2, A_diag = c(0.78, 0.45), sigma_e = 0.06)

  fit_i <- dipca(
    sim$X, s = 2, l = 2, algorithm = "I",
    n_init = 3, max_iter = 250, seed = 100,
    preproc = multivarious::center()
  )
  fit_ii <- dipca(
    sim$X, s = 2, l = 2, algorithm = "II", inner_power = 50,
    n_init = 3, max_iter = 250, seed = 100,
    preproc = multivarious::center()
  )

  obj_i <- as.numeric(sub("^[0-9]+:", "", tail(fit_i$obj_history, 1L)))
  obj_ii <- as.numeric(sub("^[0-9]+:", "", tail(fit_ii$obj_history, 1L)))
  C <- abs(stats::cor(multivarious::scores(fit_i), multivarious::scores(fit_ii)))

  expect_true(all(is.finite(fit_i$v)))
  expect_true(all(is.finite(fit_ii$v)))
  expect_gt(diag(C)[1L], 0.9)
  expect_lte(abs(obj_i - obj_ii) / max(abs(obj_i), abs(obj_ii), 1), 0.08)
})

test_that("Compact AR DiCCA agrees with classic DiCCA on an AR(1) latent signal", {
  sim <- local_exact_ar_observation(seed = 14, phi = 0.8)

  classic <- dicca(
    sim$X, s = 1, l = 1, inner = "classic",
    n_init = 3, max_iter = 250, seed = 14,
    preproc = multivarious::center()
  )
  compact <- dicca(
    sim$X, s = 1, l = 1, inner = "ar",
    n_init = 3, max_iter = 250, seed = 14,
    preproc = multivarious::center()
  )

  scores_classic <- drop(multivarious::scores(classic))
  scores_compact <- local_align_sign(scores_classic, drop(multivarious::scores(compact)))
  pred <- predict(compact, sim$X)

  expect_gt(abs(stats::cor(scores_classic, scores_compact)), 0.96)
  expect_gt(abs(stats::cor(scores_compact, sim$t)), 0.96)
  expect_true(all(is.finite(pred$scores[-1L, , drop = FALSE])))
  expect_true(all(is.finite(pred$scores_hat[-1L, , drop = FALSE])))
})

test_that("DiPCA and DiCCA remain stable across small randomized stress cases", {
  for (seed in 21:25) {
    set.seed(seed)
    Xbase <- matrix(rnorm(70 * 4), 70, 4)
    X <- cbind(
      Xbase[, 1L],
      Xbase[, 1L] + 1e-8 * rnorm(70),
      Xbase[, 2L] - Xbase[, 3L],
      10 * Xbase[, 4L],
      rep(1, 70)
    )

    fit_dipca <- dipca(X, s = 2, l = 2, n_init = 2, max_iter = 120, seed = seed)
    fit_dicca <- dicca(X, s = 2, l = 2, n_init = 2, max_iter = 120, seed = seed, inner = "classic")

    expect_true(all(is.finite(fit_dipca$v)))
    expect_true(all(is.finite(fit_dipca$betas)))
    expect_true(all(is.finite(fit_dicca$v)))
    expect_true(all(is.finite(fit_dicca$betas)))
    expect_true(all(fit_dipca$R2 >= -1e-10 & fit_dipca$R2 <= 1 + 1e-10))
    expect_true(all(fit_dicca$R2 >= -1e-10 & fit_dicca$R2 <= 1 + 1e-10))
  }
})
