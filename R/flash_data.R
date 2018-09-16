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

    # Initialize data.
    data = list(Yorig = Y, S = S, anyNA = anyNA(Y), missing = is.na(Y))

    if (anyNA(Y)) {
        # Replace missing data with 0s.
        if (any(rowSums(!data$missing) == 0)) {
            stop("Data must not have all missing rows.")
        }
        if (any(colSums(!data$missing) == 0)) {
            stop("Data must not have all missing columns.")
        }
        Y[data$missing] = 0
    }

    data$Y = Y

    if (!is.null(S) && is.matrix(S)) {
      if (nrow(S) != nrow(Y) || ncol(S) != ncol(Y)) {
        stop(paste("If S is a matrix, dimensions of S must match",
                   "dimensions of Y."))
      }
    }

    class(data) = "flash_data"

    return(data)
}


get_Yorig = function(data) {
    if (data$anyNA) {
        return(data$Yorig)
    }

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
  if (is.matrix(data$Yorig)) {
    data$Yorig = t(data$Yorig)
  }
  if (is.matrix(data$missing)) {
    data$missing = t(data$missing)
  }
  if (is.matrix(data$Y)) {
    data$Y = t(data$Y)
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
  subdata$Yorig = subdata$Yorig[row_subset, col_subset, drop = F]
  subdata$anyNA = anyNA(subdata$Yorig)
  subdata$missing = subdata$missing[row_subset, col_subset, drop = F]
  subdata$Y = subdata$Y[row_subset, col_subset, drop = F]
  class(subdata) = "flash_data"

  return(subdata)
}
