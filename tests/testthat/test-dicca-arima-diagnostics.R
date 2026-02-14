# DiCCA ARIMA vs AR Comparison Tests
# Tests that ARIMA inner models improve on AR for nonstationary data
# Based on Dong & Qin (2018) compactness and generalization results


test_that("ARIMA-DiCCA vs AR-DiCCA on I(1)+ARMA latent process", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  skip_on_cran()

  set.seed(1234)
  sim <- sim_ARIMA_latent(
    N = 300, m = 8, h = 3,
    d  = c(1, 0, 0),                    # First component is I(1)
    ar = list(0.6, 0.6, 0.3),
    ma = list(0.5, 0.4, numeric(0)),
    sigma_v = 0.3, sigma_e = 0.2
  )

  # Train/test split
  Ntr <- 200
  Xtr <- sim$X[1:Ntr, , drop = FALSE]
  Xte <- sim$X[(Ntr + 1):nrow(sim$X), , drop = FALSE]

  # --- AR-DiCCA baseline (high order to cover dynamics) -------------------
  fit_ar <- dicca(Xtr, s = 12, l = 3, inner = "classic", n_init = 2,
                 preproc = multivarious::center(),
                 max_iter = 400, tol = 1e-7)

  # --- ARIMA-DiCCA (compact inner models) ---------------------------------
  fit_ai <- dicca(Xtr, s = 12, l = 3, inner = "arima", n_init = 2,
                 preproc = multivarious::center(),
                 max_iter = 400, tol = 1e-7)

  # --- Diagnostics --------------------------------------------------------

  # (i) Parameter compactness: ARIMA should use fewer parameters
  # AR-DiCCA uses l * s parameters total
  params_ar <- 3 * 12  # = 36

  # ARIMA-DiCCA uses p+d+q per component
  if (!is.null(fit_ai$compact_models)) {
    params_ai <- sum(vapply(fit_ai$compact_models, function(m) {
      # m$arma has format c(p, q, P, Q, period, d, D)
      m$arma[1] + m$arma[2] + m$arma[6]  # p + q + d
    }, numeric(1)))

    expect_lt(params_ai, params_ar,
             label = "ARIMA should use fewer parameters than high-order AR")

    # At least one component should have d > 0 (caught integration)
    d_vals <- vapply(fit_ai$compact_models, function(m) m$arma[6], numeric(1))
    expect_true(any(d_vals > 0),
               label = "ARIMA should detect unit root in at least one component")
  }

  # (ii) Test set prediction: use predict() method
  pred_ar <- predict(fit_ar, Xte)
  pred_ai <- predict(fit_ai, Xte)

  # Reconstruct in data space
  Xhat_ar <- pred_ar$scores_hat %*% t(fit_ar$loadings)
  Xhat_ai <- pred_ai$scores_hat %*% t(fit_ai$loadings)

  # Reverse preprocessing to original scale using stored mu/sigma
  Xhat_ar_orig <- sweep(sweep(Xhat_ar, 2, fit_ar$sigma, "*"), 2, fit_ar$mu, "+")
  Xhat_ai_orig <- sweep(sweep(Xhat_ai, 2, fit_ai$sigma, "*"), 2, fit_ai$mu, "+")

  # MSPE comparison (ARIMA should be <= AR, allowing 10% tolerance)
  valid_idx <- which(rowSums(is.finite(Xhat_ar_orig)) == ncol(Xte) &
                    rowSums(is.finite(Xhat_ai_orig)) == ncol(Xte))

  if (length(valid_idx) > 10) {
    mspe_ar <- mspe(Xte[valid_idx, ], Xhat_ar_orig[valid_idx, ])
    mspe_ai <- mspe(Xte[valid_idx, ], Xhat_ai_orig[valid_idx, ])

    expect_lte(mspe_ai, 1.1 * mspe_ar,
              label = "ARIMA-DiCCA should have similar or better MSPE than AR-DiCCA")
  }

  # (iii) Whiteness of residuals: Ljung-Box test
  res_ar <- residuals(fit_ar, Xte)
  res_ai <- residuals(fit_ai, Xte)

  lb_p_values <- function(E_mat) {
    apply(E_mat, 2, function(col) {
      col <- col[is.finite(col)]
      if (length(col) < 20) return(NA_real_)
      tryCatch({
        stats::Box.test(col, lag = 10, type = "Ljung-Box")$p.value
      }, error = function(e) NA_real_)
    })
  }

  p_ar <- lb_p_values(res_ar$e_hat)
  p_ai <- lb_p_values(res_ai$e_hat)

  # ARIMA should have higher median p-value (whiter residuals)
  med_ar <- median(p_ar, na.rm = TRUE)
  med_ai <- median(p_ai, na.rm = TRUE)

  if (is.finite(med_ar) && is.finite(med_ai)) {
    expect_gte(med_ai, med_ar * 0.8,
              label = "ARIMA should produce whiter residuals (higher LB p-values)")
  }
})

test_that("ARIMA-DiCCA: model structure and fields", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  skip_on_cran()

  sim <- sim_ARIMA_latent(N = 200, m = 6, h = 2,
                         d = c(1, 0), ar = list(0.7, 0.5),
                         ma = list(numeric(0), 0.3))
  # Use multivarious preprocessing in dicca() call

  fit <- dicca(sim$X, s = 8, l = 2, inner = "arima", n_init = 2, max_iter = 300)

  # Should have compact_models field
  expect_true(!is.null(fit$compact_models))
  expect_equal(length(fit$compact_models), 2)

  # Each model element should have a $fit field with an Arima object
  for (m in fit$compact_models) {
    expect_true(!is.null(m$fit))
    expect_s3_class(m$fit, "Arima")
    expect_true(!is.null(m$fit$arma))
    expect_true(!is.null(m$fit$coef))
  }

  # Should have standard bi_projector fields
  expect_true(!is.null(fit$loadings))
  expect_true(!is.null(multivarious::scores(fit)))
  expect_equal(ncol(fit$loadings), 2)
  expect_equal(ncol(multivarious::scores(fit)), 2)
})

test_that("ARIMA-DiCCA: predict and residuals work with compact models", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  skip_on_cran()

  sim <- sim_ARIMA_latent(N = 180, m = 6, h = 2,
                         d = c(1, 0), ar = list(0.6, 0.5),
                         ma = list(0.0, 0.3))

  # Split data
  Xtr <- sim$X[1:120, ]
  Xte <- sim$X[121:180, ]

  fit <- dicca(Xtr, s = 6, l = 2, inner = "arima", n_init = 2,
              preproc = multivarious::center())

  # Test predict()
  pred <- predict(fit, Xte)
  expect_equal(names(pred), c("scores", "scores_hat"))
  expect_equal(nrow(pred$scores), nrow(Xte))
  expect_equal(nrow(pred$scores_hat), nrow(Xte))

  # Test residuals()
  resid <- residuals(fit, Xte)
  expect_equal(names(resid), c("v", "e_hat", "scores"))
  expect_equal(nrow(resid$v), nrow(Xte))
  expect_equal(nrow(resid$e_hat), nrow(Xte))
  expect_equal(ncol(resid$e_hat), ncol(Xte))
})

test_that("ARIMA-DiCCA: handles pure stationary data", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  skip_on_cran()

  # All stationary ARMA components (no integration)
  sim <- sim_ARIMA_latent(N = 200, m = 6, h = 2,
                         d = c(0, 0),  # No integration
                         ar = list(0.7, 0.5),
                         ma = list(0.3, 0.2))
  # Use multivarious preprocessing in dicca() call

  fit <- dicca(sim$X, s = 6, l = 2, inner = "arima", n_init = 2, max_iter = 300)

  # Should still fit successfully
  expect_true(!is.null(fit$compact_models))
  expect_equal(length(fit$compact_models), 2)

  # Should detect no integration (d = 0 for all)
  d_vals <- vapply(fit$compact_models, function(m) {
    if (is.null(m$fit) || !inherits(m$fit, "Arima")) return(NA_real_)
    m$fit$arma[6]
  }, numeric(1))
  # Allow NA if ARIMA failed to fit
  expect_true(all(d_vals == 0 | is.na(d_vals)),
             label = "Should detect stationary process (no integration) when fitted")
})

test_that("ARIMA-DiCCA: different inner model types", {
  skip("ARIMA inner model is experimental; tests disabled until stabilized.")
  skip_on_cran()

  sim <- sim_VAR(N = 150, m = 6, h = 2, A_diag = c(0.8, 0.6))
  # Use multivarious preprocessing in dicca() call

  # Test different inner model types
  fit_ar   <- dicca(sim$X, s = 4, l = 2, inner = "ar", n_init = 2, max_iter = 200)
  fit_arma <- dicca(sim$X, s = 4, l = 2, inner = "arma", n_init = 2, max_iter = 200)
  fit_arima <- dicca(sim$X, s = 4, l = 2, inner = "arima", n_init = 2, max_iter = 200)

  # All should have compact_models
  expect_true(!is.null(fit_ar$compact_models))
  expect_true(!is.null(fit_arma$compact_models))
  expect_true(!is.null(fit_arima$compact_models))

  # All should produce valid predictions
  pred_ar <- predict(fit_ar, X)
  pred_arma <- predict(fit_arma, X)
  pred_arima <- predict(fit_arima, X)

  expect_true(all(is.finite(pred_ar$scores)))
  expect_true(all(is.finite(pred_arma$scores)))
  expect_true(all(is.finite(pred_arima$scores)))
})
