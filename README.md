
# dipca (R)

A high-performance R implementation of **Dynamic Inner PCA (DiPCA)** with an Rcpp core.

## Install

```r
# from a local folder
install.packages("Rcpp")
remotes::install_local("dipca")
# or
devtools::install("dipca")
```

## Usage

```r
library(dipca)
fit <- dipca(X, s=1, l=3, n_init=5, max_iter=1000, tol=1e-7, algorithm="I")
pred <- dipca_predict(fit, X)
pe <- dipca_prediction_errors(fit, X)
```

## Notes
- Iteration and deflation follow Dong & Qin (2018, Table 2; Eqs 8–18).
- Convergence monitoring and Algorithm I/II follow Shin et al. (2020).
