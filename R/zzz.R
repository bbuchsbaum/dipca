.onLoad <- function(libname, pkgname) {
  if (requireNamespace("multivarious", quietly = TRUE)) {
    ns <- getNamespace("multivarious")
    registerS3method("predict", "dicca", predict.dicca, envir = ns)
    registerS3method("residuals", "dicca", residuals.dicca, envir = ns)
  }
}

