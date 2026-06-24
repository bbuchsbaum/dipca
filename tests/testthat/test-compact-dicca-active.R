
test_that("arima_dicca rejects invalid burnin and experimental arima", {
  X <- matrix(rnorm(60), 20, 3)
  expect_error(arima_dicca(X, n_comp = 1, burnin = 0), "burnin must be >= 1")

  old <- getOption("dipca.experimental_arima")
  on.exit(options(dipca.experimental_arima = old), add = TRUE)
  options(dipca.experimental_arima = FALSE)
  expect_error(
    arima_dicca(X, n_comp = 1, burnin = 2, inner = "arima"),
    "experimental and disabled by default"
  )
})

test_that("arima_dicca fits AR inner model and supports score utilities", {
  set.seed(42)
  sim <- sim_VAR(N = 200, m = 6, h = 2, A_diag = c(0.7, 0.5))

  fit <- arima_dicca(
    sim$X, n_comp = 2, burnin = 2,
    inner = "ar", center = TRUE, scale = FALSE,
    max_iter = 100L, verbose = 0L
  )

  expect_s3_class(fit, "dicca_model")
  expect_equal(dim(fit$W), c(6, 2))
  expect_equal(dim(fit$T), c(200, 2))
  expect_length(fit$models, 2)

  scores_new <- dicca_scores(fit, sim$X)
  expect_equal(dim(scores_new), c(200, 2))
  expect_true(all(is.finite(scores_new)))

  pred <- dicca_predict_scores(fit, fit$T)
  expect_equal(dim(pred), c(200, 2))
  expect_true(any(is.finite(pred[-(1:2), ])))
})

test_that("arima_dicca with scaling preprocesses data", {
  set.seed(7)
  X <- scale(matrix(rnorm(150), 50, 3), center = TRUE, scale = TRUE)

  fit <- arima_dicca(
    X, n_comp = 1, burnin = 2,
    inner = "ar", center = TRUE, scale = TRUE,
    max_iter = 80L, verbose = 0L
  )

  expect_equal(fit$info$scale, rep(1, ncol(X)))
  expect_true(all(is.finite(fit$T)))
})

test_that("dicca compact inner models integrate with S3 methods", {
  set.seed(11)
  sim <- sim_VAR(N = 180, m = 5, h = 2, A_diag = c(0.7, 0.5))

  fit <- dicca(sim$X, s = 2, l = 2, inner = "ar", n_init = 1, max_iter = 120, verbose = 0)
  expect_s3_class(fit, "dicca")
  expect_true(!is.null(fit$compact_models))
  expect_equal(length(fit$compact_models), 2)

  pred <- predict(fit, sim$X)
  expect_equal(names(pred), c("scores", "scores_hat"))
  expect_equal(dim(pred$scores), c(180, 2))

  resid <- residuals(fit, sim$X)
  expect_equal(names(resid), c("v", "e_hat", "scores"))
  expect_equal(dim(resid$e_hat), dim(sim$X))
})

test_that("dicca inner = arima uses forecast auto.arima pathway", {
  skip_if_not_installed("forecast")
  set.seed(3)
  sim <- sim_VAR(N = 160, m = 5, h = 1, A_diag = 0.7)

  old <- getOption("dipca.experimental_arima")
  on.exit(options(dipca.experimental_arima = old), add = TRUE)
  options(dipca.experimental_arima = TRUE)

  fit <- dicca(
    sim$X, s = 2, l = 1, inner = "arima",
    n_init = 1, max_iter = 80, verbose = 0
  )

  expect_s3_class(fit, "dicca")
  expect_true(!is.null(fit$compact_models))
  expect_s3_class(fit$compact_models[[1]]$fit, c("Arima", "ARIMA", "arima"))
})
