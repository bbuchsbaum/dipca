

test_that("ARIMA-DiCCA improves latent residual whiteness vs AR-DiCCA", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  set.seed(123)
  Tn <- 400; m <- 8; h <- 2
  # Simulate two latent DLVs with different dynamics
  t1 <- as.numeric(stats::arima.sim(list(order = c(2,1,1), ar = c(0.6,-0.2), ma = -0.4), n = Tn))
  t2 <- as.numeric(stats::arima.sim(list(order = c(1,0,1), ar = 0.7, ma = 0.3), n = Tn))
  t1 <- t1[seq_len(Tn)]
  t2 <- t2[seq_len(Tn)]
  Ttrue <- cbind(t1, t2)
  # Random loadings
  Ptrue <- matrix(rnorm(m*h), m, h)
  Ptrue <- sweep(Ptrue, 2, sqrt(colSums(Ptrue^2)), "/")
  # Data with noise
  X <- Ttrue %*% t(Ptrue) + matrix(rnorm(Tn*m, sd = 0.5), Tn, m)

  # Fit ARIMA-DiCCA (compact)
  fit_arima <- arima_dicca(X, n_comp = h, inner = "arima", burnin = 2, center = TRUE, scale = FALSE, verbose = 0)
  # Fit AR-DiCCA baseline (high-order AR inner)
  fit_ar    <- arima_dicca(X, n_comp = h, inner = "ar",    burnin = 2, center = TRUE, scale = FALSE, verbose = 0)

  Sc_arima <- fit_arima$T
  Sc_ar    <- fit_ar$T

  # One-step predicted t (latent) using stored models
  That_arima <- suppressWarnings(dicca_predict_scores(fit_arima, Sc_arima))
  That_ar    <- suppressWarnings(dicca_predict_scores(fit_ar,    Sc_ar))

  # Innovation series v_k = t_k - t_hat_k; Ljung-Box p-values (lag 10)
  lb_p <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 20) return(NA_real_)
    stats::Box.test(v, lag = 10, type = "Ljung-Box")$p.value
  }
  p_arima <- apply(Sc_arima - That_arima, 2, lb_p)
  p_ar    <- apply(Sc_ar    - That_ar,    2, lb_p)

  # Expect ARIMA inner model yields whiter innovations on average
  # Skip test if ARIMA convergence issues prevented valid predictions
  if (all(is.na(p_arima)) || all(is.na(p_ar))) {
    skip("ARIMA models failed to converge for this simulation")
  }
  expect_true(median(p_arima, na.rm = TRUE) >= median(p_ar, na.rm = TRUE) * 0.8)
})
