# functions related to the flash data object

#' @title set up data for reading into flash
#' @param Y an n by p data matrix
#' @details Y can have missing data, but no column or row can be entirely missing
#' The flash data object contains flags for dealing with missing data
#' and a (naively) imputed version of the original data matrix so that
#' i) some of the initializatinon methods used by flash (eg svd) do not fail
#' ii) data$Y * data$missing is 0 if the original data were missing
#' @return a flash data object
#' @export
set_flash_data = function(Y, init=c("softimpute","mean")){
  init = match.arg(init)
  data = list(Y=NULL, anyNA=FALSE, missing = FALSE) # initialize data

  if(anyNA(Y)){ # deal with missing data: set flags and impute
    data$missing = is.na(Y)
    data$anyNA=TRUE
    data$Yorig = Y # save original data with missingness

    if(any(rowSums(!data$missing)==0)){
      stop("data must not have all missing rows")}
    if(any(colSums(!data$missing)==0)){
      stop("data must not have all missing columns")}
    Y[data$missing] = 0
  }
  data$Y=Y
  return(data)
}


