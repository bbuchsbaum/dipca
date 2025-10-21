
context("dicca core")

simulate_latent_VAR <- function(seed=0, T=1200, l=3, m=5, s=1, noise_x=0.1, noise_t=1.0) {
  set.seed(seed)
  A <- matrix(c(0.5205, 0.1022, 0.0599,
                0.5367, -0.0139, 0.4159,
                0.0412, 0.6054, 0.3874), 3, 3, byrow=TRUE)
  c0 <- c(0.5205, 0.5367, 0.0412)
  P <- matrix(c(0.4316, 0.1723, -0.0574,
                0.1202, -0.1463, 0.5348,
                0.2483, 0.1982, 0.4797,
                0.1151, 0.1557, 0.3739,
                0.2258, 0.5461, -0.0424), 5, 3, byrow=TRUE)
  t <- matrix(0, T, l); v <- matrix(rnorm(T*l, sd=sqrt(noise_t)), T, l)
  for (k in 2:T) t[k,] <- c0 + A %*% t[k-1,] + v[k,]
  e <- matrix(rnorm(T*m, sd=sqrt(noise_x)), T, m)
  X <- t %*% t(P) + e
  list(X=X, t=t, P=P)
}

test_that("DiCCA recovers dynamic subspace on VAR(1)", {
  sim <- simulate_latent_VAR(T=1000)
  fit <- dicca(sim$X, s=1, l=3, n_init=2, max_iter=600, tol=1e-7, inner="classic")
  t_est <- multivarious::scores(fit)  # Use multivarious accessor
  U <- qr.Q(qr(sim$t[201:1000,]))
  V <- qr.Q(qr(t_est[201:1000,]))
  svals <- svd(t(U) %*% V, nu=0, nv=0)$d
  expect_gt(svals[1], 0.8)
})

test_that("dicca predict and errors shapes", {
  sim <- simulate_latent_VAR(T=250)
  fit <- dicca(sim$X, s=1, l=2, inner="classic")
  pr <- predict(fit, sim$X)  # Use S3 predict method
  expect_equal(dim(pr$scores), c(250, 2))
  errs <- residuals(fit, sim$X)  # Use S3 residuals method
  expect_equal(dim(errs$e_hat), dim(sim$X))
})

