#' @title Print summary of flash object
#'
#' @description \code{print} method for class \code{'flash'}.
#'
#' @param x A flash object.
#'
#' @param \dots Further arguments passed to or from other methods.
#'
#' @export
#'
print.flash = function(x, ...) {
  cat("Summary of flash object:\n")
  cat(sprintf("  Number of factor/loading pairs: %d\n", x$nfactors))

  cat("  Proportion of variance explained:\n")
  rows_suppressed = FALSE
  for (k in 1:length(x$pve)) {
    if (x$pve[k] < 0.001) {
      rows_suppressed = TRUE
    } else {
      cat(sprintf("    Factor/loading %d: %0.3f\n", k, x$pve[k]))
    }
  }
  if (rows_suppressed) {
    cat(paste("      *Factor/loadings with PVE < 0.001 are omitted from",
              "this summary.\n"))
  }

  if (!is.na(x$objective)) {
    cat(sprintf("  Value of objective function: %0.3f\n", x$objective))
  }
}
