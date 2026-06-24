## Test environments

* Local: macOS Sonoma 14.3, R 4.5.1 (aarch64), Homebrew clang 20.1.8.

## R CMD check results

0 errors | 0 warnings | 0 notes

Local macOS check (`R CMD check --as-cran`, commit `53d1f82`): 0 errors, 1 warning, 1 NOTE.

* **Warning**: compiler warning from R headers with Homebrew clang (`-Wfixed-enum-extension`); environment-specific, not emitted on CRAN builders.
* **NOTE**: "New submission" from CRAN incoming feasibility — expected for a first submit.

## Cross-platform checks

* **R-hub**: `.github/workflows/rhub.yaml` added locally. After `git push`, run `rhub::rhub_check(platforms = c("windows", "macos-arm64", "ubuntu-release", "ubuntu-clang"))`. Alternatively, `rhub::rc_new_token()` then `rhub::rc_submit(path = "dipca_0.1.0.tar.gz", confirmation = TRUE)`.
* **win-builder**: upload `dipca_0.1.0.tar.gz` to <https://win-builder.r-project.org/upload.html> or email to win-builder@r-project.org.
* **macOS builder**: upload `dipca_0.1.0.tar.gz` at <https://mac.r-project.org/mac-mini4/submit.html>.

## Test coverage

* `testthat` suite: 262 tests passed (6 skipped for experimental ARIMA paths).
* Package coverage: ~90%.

## Downstream dependencies

This is a new submission with no reverse dependencies.

## Additional notes

* Depends on `multivarious` (>= 0.3.0), which is available on CRAN.
* Experimental `inner = "arima"` for compact DiCCA is disabled by default and guarded by `options(dipca.experimental_arima = TRUE)`.
