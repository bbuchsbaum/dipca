
test_that(".center_scale handles centering, scaling, and zero-variance columns", {
  X <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, ncol = 3)

  centered <- dipca:::.center_scale(X, center = TRUE, scale = FALSE)
  expect_equal(centered$center, colMeans(X))
  expect_equal(centered$scale, rep(1, ncol(X)))
  expect_equal(as.numeric(centered$X), as.numeric(scale(X, center = TRUE, scale = FALSE)))

  scaled_only <- dipca:::.center_scale(X, center = FALSE, scale = TRUE)
  expect_equal(scaled_only$center, rep(0, ncol(X)))
  col_sds <- sqrt(colSums(X^2) / max(1, nrow(X) - 1))
  expect_equal(scaled_only$X, sweep(X, 2, col_sds, "/"), tolerance = 1e-10)

  both <- dipca:::.center_scale(X, center = TRUE, scale = TRUE)
  expect_equal(both$center, colMeans(X))
  expect_true(all(abs(colMeans(both$X)) < 1e-10))
  expect_true(all(abs(apply(both$X, 2, sd) - 1) < 1e-10))

  const_col <- cbind(X[, 1], rep(5, nrow(X)))
  const_scaled <- dipca:::.center_scale(const_col, center = TRUE, scale = TRUE)
  expect_equal(const_scaled$scale[2], 1)
})

test_that(".form_blocks builds lag blocks and falls back when C++ errors", {
  X <- matrix(rnorm(30), 10, 3)
  s <- 2L

  blocks <- dipca:::.form_blocks(X, s)
  expect_length(blocks, s + 1L)
  expect_equal(nrow(blocks[[1]]), nrow(X) - s)
  expect_equal(ncol(blocks[[1]]), ncol(X))

  expect_error(dipca:::.form_blocks(X, nrow(X)), "not enough rows")

  local_mocked_bindings(
    form_blocks_cpp = function(X, s) stop("forced R fallback"),
    .package = "dipca"
  )
  fallback <- dipca:::.form_blocks(X, s)
  expect_length(fallback, s + 1L)
  expect_equal(dim(fallback[[1]]), dim(blocks[[1]]))
})

test_that(".make_t_lags constructs aligned lag matrices", {
  t <- seq_len(10)
  lag <- dipca:::.make_t_lags(t, 3L)
  expect_equal(lag$ts1, t[4:10])
  expect_equal(nrow(lag$Ts), 7L)
  expect_equal(ncol(lag$Ts), 3L)
  expect_equal(lag$Ts[1, 1], t[3])
  expect_equal(lag$Ts[1, 3], t[1])
})

test_that(".solve_normal solves systems with ridge and SVD fallback", {
  A <- diag(c(2, 3))
  b <- c(4, 9)
  expect_equal(dipca:::.solve_normal(A, b), c(2, 3))

  ridge <- dipca:::.solve_normal(A, b, ridge = 0.5)
  expect_equal(ridge, solve(A + diag(0.5, 2), b))

  B <- matrix(c(1, 2, 3, 4), 2, 2)
  mat_sol <- dipca:::.solve_normal(A, B)
  expect_equal(mat_sol, solve(A, B))

  singular <- matrix(c(1, 1, 1, 1), 2, 2)
  sv_sol <- dipca:::.solve_normal(singular, c(2, 2))
  expect_length(sv_sol, 2L)
  expect_true(all(is.finite(sv_sol)))
})

test_that(".orth_resid removes projection onto t", {
  set.seed(1)
  X <- matrix(rnorm(20), 10, 2)
  t <- rnorm(10)
  resid <- dipca:::.orth_resid(X, t)
  p <- as.numeric(crossprod(X, t) / sum(t^2))
  expect_equal(resid, X - tcrossprod(t, p))
})

test_that(".as_arma_meta packs ARIMA order metadata", {
  expect_equal(
    dipca:::.as_arma_meta(p = 2L, d = 1L, q = 3L, P = 0L, Q = 0L, m = 1L, D = 0L),
    c(2L, 3L, 0L, 0L, 1L, 1L, 0L)
  )
})
