
context("dipca core")

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

test_that("recovers dynamic subspace", {
  sim <- simulate_latent_VAR(T=1200)
  fit <- dipca(sim$X, s=1, l=3, n_init=3, max_iter=800, tol=1e-7, algorithm="I")
  t_est <- multivarious::scores(fit)  # Use multivarious accessor
  U <- qr.Q(qr(sim$t[201:1200,]))
  V <- qr.Q(qr(t_est[201:1200,]))
  svals <- svd(t(U) %*% V, nu=0, nv=0)$d
  expect_gt(svals[1], 0.9)
})

test_that("scores are near-orthogonal", {
  sim <- simulate_latent_VAR(T=900)
  fit <- dipca(sim$X, s=1, l=3, n_init=2, max_iter=600, tol=1e-7)
  S <- multivarious::scores(fit)  # Use multivarious accessor
  G <- crossprod(S)/nrow(S)
  diag(G) <- 0
  expect_lt(max(abs(G)), 1e-2)
})

autocorr <- function(x, lag) {
  x <- x - mean(x)
  if (lag >= length(x)) return(0)
  sum(head(x, -lag) * tail(x, -lag)) / sum(x*x)
}

test_that("prediction errors whitening tendency", {
  sim <- simulate_latent_VAR(T=900)
  fit <- dipca(sim$X, s=1, l=3, n_init=2, max_iter=600, tol=1e-7)
  pe <- residuals(fit, sim$X)  # Use S3 residuals method
  E <- pe$e_hat[(fit$lag_order+1):nrow(sim$X),]
  acs <- apply(E, 2, function(col) abs(autocorr(col, 1)))
  expect_lt(max(acs, na.rm=TRUE), 0.4)  # practical threshold
})

test_that("algorithm I vs II close on small problem", {
  sim <- simulate_latent_VAR(T=400)
  f1 <- dipca(sim$X, s=1, l=2, algorithm="I", n_init=2, max_iter=500, tol=1e-7)
  f2 <- dipca(sim$X, s=1, l=2, algorithm="II", n_init=2, max_iter=500, tol=1e-7, inner_power=100, inner_tol=1e-10)
  obj1 <- as.numeric(strsplit(tail(f1$obj_history, 1), ":")[[1]][2])
  obj2 <- as.numeric(strsplit(tail(f2$obj_history, 1), ":")[[1]][2])
  expect_lt(abs(obj1 - obj2)/max(1, abs(obj2)), 0.05)
})

test_that("predict and transform shapes", {
  sim <- simulate_latent_VAR(T=250)
  fit <- dipca(sim$X, s=1, l=2)
  sc <- multivarious::project(fit, sim$X)  # Use multivarious project
  expect_equal(dim(sc), c(250, 2))
  pe <- residuals(fit, sim$X)  # Use S3 residuals method
  expect_equal(dim(pe$e_hat), dim(sim$X))
})
