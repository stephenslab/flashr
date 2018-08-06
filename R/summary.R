#' @title Summarize flash object
#'
#' @description \code{summary} method for class \code{'flash'}.
#'
#' @param f A flash fit object.
#'
#' @param data A matrix or flash data object created using
#'   \code{flash_set_data}.
#'
#' @param \dots Further arguments passed to or from other methods.
#'
#' @return A list with the following elements.
#'
#'   \item{nfactors}{The number of factors in the flash object. See
#'     \code{\link{flash_get_nfactors}}.}
#'
#'   \item{pve}{The proportion of variance explained by each factor. See
#'     \code{\link{flash_get_pve}}.}
#'
#'   \item{fitted.values}{The estimated LF' matrix. See
#'     \code{\link{flash_get_fitted_values}}.}
#'
#'   \item{ldf}{Standardized loadings, factors, and weights. See
#'     \code{\link{flash_get_ldf}}.}
#'
#'   \item{objective}{The value of the objective function attained by the
#'     fit. See \code{\link{flash_get_objective}}.}
#'
#' @export
#' 
summary.flash = function(f, data, ...) {
  if (class(f) != "flash") {
    stop("Input argument f must be an instance of class \"flash\".")
  }
  if (missing(data)) {
    objective = NA
  } else {
    objective = flash_get_objective(data, f)
  }

  out = list(nfactors = flash_get_nfactors(f),
             pve = flash_get_pve(f, drop_zero_factors=TRUE),
             fitted.values = flash_get_fitted_values(f),
             ldf = flash_get_ldf(f, drop_zero_factors=TRUE),
             objective = objective)

  class(out) <- c("summary.flash","list")
  return(out)
}


#' @importFrom methods is
#'
#' @export
#' 
print.summary.flash = function(summary, digits = 3, ...) {
  if (!is(summary, "summary.flash")) {
    stop("Input must be an instance of class \"summary.flash\".")
  }

  with (summary, {
    cat("Summary of flash object:\n")
    cat(sprintf("  Number of factor/loading pairs: %d\n", nfactors))

    cat("  Proportion of variance explained:\n")
    rows_suppressed = FALSE
    for (k in 1:length(pve)) {
      if (pve[k] < 0.001) {
        rows_suppressed = TRUE
      } else {
        cat(sprintf("    Factor/loading %d: %0.3f\n", k, pve[k]))
      }
    }
    if (rows_suppressed) {
      cat(paste("      *Factor/loadings with PVE < 0.001 are omitted from",
                "this summary.\n"))
    }

    if (!is.na(objective)) {
      cat(sprintf("  Value of objective function: %0.3f\n", objective))
    }
    cat(paste("\n*See help(summary.flash) for details about how to\n",
              "extract other useful information from flash objects.\n"))
  })

  return(invisible(summary))
}

print.flash <- function (f, data, digits = 3, ...) {
  print(summary(f, data), digits, ...)
}
