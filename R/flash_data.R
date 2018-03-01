# Functions related to the flash data object.

#' @title Set up data for reading into flash.
#'
#' @param Y an n by p data matrix
#'
#' @param S An n by p matrix of the standard errors of observations of
#'   Y. (Can be a scalar if all elements of the matrix are equal.)
#'   Currently S is ignored as we have not yet implemented methods for
#'   this.
#'
#' @details Y can have missing data, but no column or row can be
#' entirely missing. The flash data object contains flags for dealing
#' with missing data and a (naively) imputed version of the original
#' data matrix so that i) some of the initializatinon methods used by
#' flash (eg svd) do not fail ii) data$Y * data$missing is 0 if the
#' original data were missing.
#'
#' @return A flash data object.
#'
#' @export
#'
flash_set_data = function(Y, S = 0) {

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
    class(data) = "flash_data"
    return(data)
}

get_Yorig = function(data) {
    if (data$anyNA) {
        return(data$Yorig)
    }
    return(data$Y)
}
