# Functions related to the flash data object.

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

    # initialize data
    data = list(Yorig = Y, S = S, anyNA = anyNA(Y), missing = is.na(Y))

    if (anyNA(Y)) {
        # replace missing data with 0s
        if (any(rowSums(!data$missing) == 0)) {
            stop("data must not have all missing rows")
        }
        if (any(colSums(!data$missing) == 0)) {
            stop("data must not have all missing columns")
        }
        Y[data$missing] = 0
    }

    data$Y = Y

    if (!is.null(S)) {
      data$S = matrix(S, nrow=nrow(Y), ncol=ncol(Y))
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
