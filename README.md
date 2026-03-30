# dipca

`dipca` is an R package for extracting dynamic latent structure from
multivariate time series.

It implements three related methods:

| Method | Function | Use case |
| --- | --- | --- |
| Dynamic Inner PCA | `dipca()` | Unsupervised dynamic components with predictable latent scores |
| Dynamic Inner CCA | `dicca()` | Components ordered by predictability / canonical correlation |
| Dynamic Inner PLS | `dipls()` | Supervised dynamic modeling from an input block `X` to an output block `Y` |

The package is built on top of
[`multivarious`](https://bbuchsbaum.github.io/multivarious/) so fitted models
work with standard projection-style workflows such as `scores()`, `project()`,
`predict()`, `reconstruct()`, and `residuals()`.

## Why use it?

Standard PCA is static: it finds directions of high variance, but it ignores
time. `dipca` focuses on latent directions that are dynamically structured and
predictable from their own past, which is often the quantity of interest in
process monitoring, temporal factor extraction, and dynamic regression.

The package currently provides:

- fast Rcpp-backed implementations of DiPCA, DiCCA, and DiPLS
- classic autoregressive and compact inner-model variants for DiCCA
- prediction and residual methods for dynamic whitening and forecasting
- pkgdown articles with worked examples and ground-truth simulations

## Installation

`dipca` is not on CRAN. Install it from GitHub together with its non-CRAN
dependency `multivarious`:

```r
install.packages("remotes")

remotes::install_github("bbuchsbaum/multivarious")
remotes::install_github("bbuchsbaum/dipca")
```

For local development:

```r
remotes::install_local("path/to/multivarious")
remotes::install_local("path/to/dipca")
```

## Quick Start

```r
library(dipca)
library(multivarious)

set.seed(42)
N <- 300
m <- 6

# Two latent AR(1) processes
t1 <- arima.sim(list(ar = 0.8), n = N)
t2 <- arima.sim(list(ar = 0.4), n = N)

# Mix them into six observed variables
P <- qr.Q(qr(matrix(rnorm(m * 2), m, 2)))
X <- cbind(t1, t2) %*% t(P) + matrix(rnorm(N * m, sd = 0.3), N, m)

# Fit dynamic inner PCA
fit <- dipca(X, s = 1, l = 2, n_init = 3)

# Latent scores
scores(fit)

# One-step-ahead latent predictions
pred <- predict(fit, X)

# Dynamic residuals / whitening errors
res <- residuals(fit, X)
```

## Main Interfaces

### `dipca()`

Use `dipca()` when you want unsupervised dynamic factors that balance
variance and temporal predictability.

```r
fit <- dipca(X, s = 1, l = 3, n_init = 5, max_iter = 1000, tol = 1e-7)
```

### `dicca()`

Use `dicca()` when you want components ranked by predictability.

```r
fit <- dicca(X, s = 1, l = 3, inner = "classic", n_init = 5)
fit$R2
```

Compact inner models are also available:

```r
fit_ar <- dicca(X, s = 1, l = 3, inner = "ar")
fit_arma <- dicca(X, s = 1, l = 3, inner = "arma")
```

### `dipls()`

Use `dipls()` when you have an input block `X` and an output block `Y` and
want dynamic latent regression.

```r
fit <- dipls(X, Y, n_comp = 2, s = 2, mode = "fir", verbose = 0)
Y_hat <- predict(fit, X)
```

## Documentation

- Package site: <https://bbuchsbaum.github.io/dipca/>
- Getting started article: <https://bbuchsbaum.github.io/dipca/articles/dipca.html>
- Dynamic-components walkthrough: <https://bbuchsbaum.github.io/dipca/articles/extracting-dynamic-components.html>
- Function reference: <https://bbuchsbaum.github.io/dipca/reference/index.html>

## References

- Dong, Y., and Qin, S. J. (2018). Dynamic-inner principal component analysis
  for dynamic process monitoring. *Journal of Process Control*.
- Dong, Y., and Qin, S. J. (2018). Dynamic-inner canonical correlation and
  causality analysis for high-dimensional time series data. *IFAC-PapersOnLine*.
- Dong, Y., and Qin, S. J. (2018). Dynamic-inner partial least squares for
  dynamic regression and quality prediction. *Journal of Process Control*.
