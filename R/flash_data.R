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
    if(any(rowSums(!data$missing)==0)){
      stop("data must not have all missing rows")}
    if(any(colSums(!data$missing)==0)){
      stop("data must not have all missing columns")}

    data$anyNA=TRUE
    if(init=="mean"){Y[data$missing] = mean(Y,na.rm=TRUE)}
    if(init=="softimpute"){
      m= softImpute::softImpute(Y, rank.max = 1,type = "als",lambda = 30)
      simple_impute = outer(as.vector(m$u),as.vector(m$d * m$v))
      Y[data$missing] = simple_impute[data$missing] }
  }
  data$Y=Y
  return(data)
}


