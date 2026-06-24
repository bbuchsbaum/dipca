
test_that("form_blocks_cpp matches R fallback structure", {
  X <- matrix(rnorm(40), 10, 4)
  s <- 2L
  cpp_blocks <- dipca:::form_blocks_cpp(X, s)

  local_mocked_bindings(
    form_blocks_cpp = function(X, s) stop("forced R fallback"),
    .package = "dipca"
  )
  r_blocks <- dipca:::.form_blocks(X, s)

  expect_length(cpp_blocks, s + 1L)
  for (i in seq_along(cpp_blocks)) {
    expect_equal(dim(cpp_blocks[[i]]), dim(r_blocks[[i]]))
    expect_equal(as.matrix(cpp_blocks[[i]]), r_blocks[[i]], tolerance = 1e-10)
  }
})

test_that("dipca_component_cpp converges for algorithm I and II", {
  set.seed(5)
  sim <- sim_VAR(N = 200, m = 4, h = 2, A_diag = c(0.7, 0.5))
  X <- sim$X
  s <- 2L
  w0 <- rnorm(ncol(X))
  w0 <- w0 / sqrt(sum(w0^2))

  res_i <- dipca:::dipca_component_cpp(
    X, s, w0, "I", tol = 1e-6, max_iter = 200L,
    inner_power = 50L, inner_tol = 1e-8, verbose = 0L
  )
  expect_named(res_i, c("w", "beta", "obj", "iters", "R2", "obj_hist"))
  expect_length(res_i$w, ncol(X))
  expect_length(res_i$beta, s)
  expect_true(res_i$obj > 0)
  expect_true(res_i$iters >= 1L)

  res_ii <- dipca:::dipca_component_cpp(
    X, s, w0, "II", tol = 1e-6, max_iter = 200L,
    inner_power = 30L, inner_tol = 1e-8, verbose = 0L
  )
  expect_true(res_ii$obj > 0)
  expect_true(all(is.finite(res_ii$w)))
})

test_that("dicca_component_cpp returns finite weights on synthetic data", {
  set.seed(8)
  sim <- sim_VAR(N = 180, m = 5, h = 2, A_diag = c(0.8, 0.6))
  X <- sim$X
  s <- 2L
  w0 <- rnorm(ncol(X))
  w0 <- w0 / sqrt(sum(w0^2))

  res <- dipca:::dicca_component_cpp(
    X, s, w0, tol = 1e-6, max_iter = 150L, verbose = 0L
  )
  expect_named(res, c("w", "beta", "obj", "iters", "R2", "obj_hist"))
  expect_true(all(is.finite(res$w)))
  expect_true(res$obj > 0)
})

test_that("dipca forwards component extraction to the C++ core", {
  set.seed(12)
  calls <- list()
  local_mocked_bindings(
    dipca_component_cpp = function(X, s, w0, algorithm, tol, max_iter,
                                   inner_power, inner_tol, verbose) {
      calls[[length(calls) + 1L]] <<- list(
        algorithm = algorithm,
        inner_power = inner_power,
        inner_tol = inner_tol
      )
      list(
        w = as.numeric(w0) / sqrt(sum(w0^2)),
        beta = rep(1 / sqrt(s), s),
        obj = length(calls),
        iters = 2L,
        R2 = 0.25,
        obj_hist = c(0.5, 1)
      )
    },
    .package = "dipca"
  )

  fit <- dipca(
    matrix(rnorm(80), 20, 4),
    s = 2, l = 1, n_init = 2,
    algorithm = "II", inner_power = 7L, inner_tol = 1e-5,
    preproc = multivarious::pass()
  )

  expect_length(calls, 2L)
  expect_equal(vapply(calls, function(x) x$algorithm, character(1)), rep("II", 2))
  expect_equal(vapply(calls, function(x) x$inner_power, integer(1)), rep(7L, 2))
  expect_equal(vapply(calls, function(x) x$inner_tol, numeric(1)), rep(1e-5, 2))
  expect_equal(fit$iters_per_component, 2L)
})

test_that("classic dicca forwards component extraction to the C++ core", {
  set.seed(13)
  calls <- 0L
  local_mocked_bindings(
    dicca_component_cpp = function(X, s, w0, tol, max_iter, verbose) {
      calls <<- calls + 1L
      list(
        w = as.numeric(w0) / sqrt(sum(w0^2)),
        beta = rep(1 / sqrt(s), s),
        obj = calls,
        iters = 3L,
        R2 = 0.5,
        obj_hist = c(0.5, 1)
      )
    },
    .package = "dipca"
  )

  fit <- dicca(
    matrix(rnorm(80), 20, 4),
    s = 2, l = 1, n_init = 3,
    preproc = multivarious::pass(),
    inner = "classic"
  )

  expect_equal(calls, 3L)
  expect_equal(fit$iters_per_component, 3L)
  expect_equal(fit$R2, 0.5)
})
