#' @title Set flash data object
#'
#' @description Set up data for reading into flash.
#'
#' @param Y An n by p data matrix.
#'
#' @param S An n by p matrix of the standard errors of the observations in
#'   Y. (Can be a scalar if all standard errors are equal.) If
#'   \code{S = NULL}, then the standard errors will be estimated during
#'   fitting.
#'
#' @details Y can have missing data, but no column or row can be
#' entirely missing. The flash data object contains flags for dealing
#' with missing data and a (naively) imputed version of the original
#' data matrix so that i) some of the initialization methods used by
#' flash (e.g., \code{svd}) do not fail; ii) \code{data$Y * data$missing}
#' is zero if the original data were missing.
#'
#' @return A flash data object.
#'
#' @export
#'
flash_set_data = function(Y, S = NULL) {

  data = list(Y = Y, S = S, anyNA = anyNA(Y))

  if (data$anyNA) {
    data$missing = is.na(Y)
  }

  # TODO: move this to other parameter checks
  if (!is.null(S) && is.matrix(S)) {
    if (nrow(S) != nrow(Y) || ncol(S) != ncol(Y)) {
      stop(paste("If S is a matrix, dimensions of S must match",
                 "dimensions of Y."))
    }
  } else if (!is.null(S) && length(S) != 1) {
    stop("S must be a matrix or a scalar.")
  } else {
    if (requireNamespace("ebnm", quietly = TRUE) &&
        packageVersion("ebnm") < "0.1.13") {
      # Earlier versions of ebnm do not support scalar arguments for S
      data$S = matrix(S, nrow = nrow(Y), ncol = ncol(Y))
    }
  }

  class(data) = "flash_data"

  return(data)
}

get_Yorig = function(data) {
  return(data$Y)
}


# @title Transpose a flash data object
#
# @param f The flash data object.
#
# @return A new flash data object, with the matrices of the original
#   flash data object transposed.
#
flash_transpose_data = function(data) {
  if (is.matrix(data$Y)) {
    data$Y = t(data$Y)
  }
  if (is.matrix(data$missing)) {
    data$missing = t(data$missing)
  }
  if (is.matrix(data$S)) {
    data$S = t(data$S)
  }

  return(data)
}


# @title Subset a flash data object
#
# @param f A flash fit object.
#
# @param row_subset The subset of rows to be retained.
#
# @param col_subset The subset of columns to be retained.
#
# @return A subsetted flash data object.
#
flash_subset_data = function(data, row_subset = NULL, col_subset = NULL) {
  if (is.null(row_subset)) {
    row_subset = 1:nrow(data$Y)
  }
  if (is.null(col_subset)) {
    col_subset = 1:ncol(data$Y)
  }

  subdata = data
  subdata$Y = subdata$Y[row_subset, col_subset, drop = F]
  subdata$anyNA = anyNA(subdata$Y)

  if (subdata$anyNA) {
    subdata$missing = subdata$missing[row_subset, col_subset, drop = F]
  } else {
    subdata$missing = NULL
  }

  class(subdata) = "flash_data"

  return(subdata)
}
