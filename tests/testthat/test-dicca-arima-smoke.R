test_that("dicca arima mode requires experimental opt-in", {
  old <- getOption("dipca.experimental_arima")
  on.exit(options(dipca.experimental_arima = old), add = TRUE)
  options(dipca.experimental_arima = FALSE)

  X <- matrix(rnorm(120 * 5), 120, 5)
  expect_error(
    dicca(X, s = 2, l = 2, inner = "arima", n_init = 1, max_iter = 50),
    "experimental and disabled by default"
  )
})

test_that("dicca arima mode fits and predicts when enabled", {
  old <- getOption("dipca.experimental_arima")
  on.exit(options(dipca.experimental_arima = old), add = TRUE)
  options(dipca.experimental_arima = TRUE)

  set.seed(123)
  n <- 160
  l <- 2
  m <- 6

  t1 <- as.numeric(arima.sim(list(order = c(1, 1, 1), ar = 0.6, ma = -0.3), n = n))
  t2 <- as.numeric(arima.sim(list(order = c(1, 0, 1), ar = 0.5, ma = 0.2), n = n))
  t1 <- t1[seq_len(n)]
  t2 <- t2[seq_len(n)]
  Ttrue <- cbind(t1, t2)

  Ptrue <- matrix(rnorm(m * l), m, l)
  Ptrue <- sweep(Ptrue, 2, sqrt(colSums(Ptrue^2)), "/")
  X <- Ttrue %*% t(Ptrue) + matrix(rnorm(n * m, sd = 0.2), n, m)

  fit <- dicca(X, s = 2, l = l, inner = "arima", n_init = 1, max_iter = 200)
  expect_s3_class(fit, "dicca")
  expect_true(!is.null(fit$compact_models))
  expect_equal(length(fit$compact_models), l)

  pred <- predict(fit, X)
  expect_equal(names(pred), c("scores", "scores_hat"))
  expect_equal(dim(pred$scores), c(n, l))
  expect_equal(dim(pred$scores_hat), c(n, l))

  res <- residuals(fit, X)
  expect_equal(names(res), c("v", "e_hat", "scores"))
  expect_equal(dim(res$v), c(n, l))
  expect_equal(dim(res$e_hat), c(n, m))
})
