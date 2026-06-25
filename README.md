# dipca

Extract **dynamic latent structure** from multivariate time series with
Dynamic Inner PCA (DiPCA), CCA (DiCCA), and PLS (DiPLS).

| Method | Function | Use case |
| --- | --- | --- |
| DiPCA | `dipca()` | Unsupervised components with predictable latent scores |
| DiCCA | `dicca()` | Components ordered by predictability / canonical correlation |
| DiPLS | `dipls()` | Supervised dynamic modeling from input `X` to output `Y` |

Fitted models integrate with
[`multivarious`](https://github.com/bbuchsbaum/multivarious) for
`scores()`, `project()`, `predict()`, `reconstruct()`, and `residuals()`.

## Why use it?

Standard PCA finds high-variance directions but ignores temporal structure.
`dipca` targets latent directions that are **predictable from their own past**,
which is often what matters for process monitoring, temporal factor extraction,
and dynamic regression.

- Rcpp-backed DiPCA, DiCCA, and DiPLS
- Classic VAR and compact AR/ARMA inner models for DiCCA
- Prediction, reconstruction, and dynamic whitening via `predict()` and `residuals()`
- Vignettes with worked examples and ground-truth simulations

## Installation

```r
install.packages("remotes")
remotes::install_github("bbuchsbaum/dipca")
```

Requires [`multivarious`](https://github.com/bbuchsbaum/multivarious) (>= 0.3.0),
available on CRAN.

## Quick start

```r
library(dipca)
library(multivarious)

set.seed(42)
N <- 300; m <- 6

t1 <- arima.sim(list(ar = 0.8), n = N)
t2 <- arima.sim(list(ar = 0.4), n = N)
P <- qr.Q(qr(matrix(rnorm(m * 2), m, 2)))
X <- cbind(t1, t2) %*% t(P) + matrix(rnorm(N * m, sd = 0.3), N, m)

fit <- dipca(X, s = 1, l = 2, n_init = 3)

scores(fit)                 # latent scores
predict(fit, X)             # one-step-ahead score forecasts
stats::residuals(fit, X)    # dynamic whitening errors
```

## Main interfaces

**`dipca()`** — unsupervised dynamic factors balancing variance and temporal predictability:

```r
fit <- dipca(X, s = 1, l = 3, n_init = 5, max_iter = 1000, tol = 1e-7)
```

**`dicca()`** — components ranked by R² (classic inner model):

```r
fit <- dicca(X, s = 1, l = 3, inner = "classic", n_init = 5)
fit$R2
```

Compact per-component inner models: `inner = "ar"` or `"arma"`.

**`dipls()`** — dynamic latent regression from `X` to `Y`:

```r
fit <- dipls(X, Y, n_comp = 2, s = 2, mode = "fir", verbose = 0)
predict(fit, X)
```

## Documentation

- `vignette("dipca")` — getting started
- `vignette("extracting-dynamic-components")` — ground-truth simulation walkthrough
- `?dipca`, `?dicca`, `?dipls` — function reference

## References

- Dong, Y., & Qin, S. J. (2018). A new dynamic PCA method for dynamic data modeling and process monitoring. *Journal of Process Control*, 67, 1-11. [doi:10.1016/j.jprocont.2017.05.002](https://doi.org/10.1016/j.jprocont.2017.05.002)
- Dong, Y., & Qin, S. J. (2018). Dynamic-inner canonical correlation and causality analysis for high dimensional time series data. *IFAC-PapersOnLine*, 51(18), 476-481. [doi:10.1016/j.ifacol.2018.09.379](https://doi.org/10.1016/j.ifacol.2018.09.379)
- Dong, Y., & Qin, S. J. (2018). Dynamic-inner partial least squares for dynamic data modeling and monitoring. *Journal of Process Control*, 69, 1-12. [doi:10.1016/j.jprocont.2018.04.006](https://doi.org/10.1016/j.jprocont.2018.04.006)

## License

GPL (>= 3)
