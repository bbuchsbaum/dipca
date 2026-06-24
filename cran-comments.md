## Test environments

* Local: macOS Sonoma 14.3, R 4.5.1 (aarch64), Homebrew clang 20.1.8.
* R-hub (see below).

## R CMD check results

0 errors | 0 warnings | 0 notes

Local macOS check (`R CMD check --as-cran`): 0 errors, 1 warning, 1 NOTE.

* **Warning**: compiler warning from R headers with Homebrew clang (`-Wfixed-enum-extension`); environment-specific, not emitted on CRAN builders.
* **NOTE**: "New submission" from CRAN incoming feasibility — expected for a first submit.

## Cross-platform checks

* `rhub::rhub_check()` — submitted from this branch (results pending).
* win-builder and macOS builder: tarball built locally; submit manually if R-hub is unavailable.

## Test coverage

* `testthat` suite: 262 tests passed (6 skipped for experimental ARIMA paths).
* Package coverage: ~90%.

## Downstream dependencies

This is a new submission with no reverse dependencies.

## Additional notes

* Depends on `multivarious` (>= 0.3.0), which is available on CRAN.
* Experimental `inner = "arima"` for compact DiCCA is disabled by default and guarded by `options(dipca.experimental_arima = TRUE)`.
