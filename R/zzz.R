#' @useDynLib scPAS, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Dynamic library initialization
}

.onUnload <- function(libpath) {
  # Clean up when package is unloaded
  library.dynam.unload("scPAS", libpath)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "scPAS v", utils::packageVersion(pkgname), 
    " - Single-Cell Phenotype-Associated Subpopulation Identifier"
  )
}
