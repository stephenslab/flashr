#' @title  Update precision parameter
#' @details Updates precision estimate to increase the objective F
#' @param data a flash data object
#' @param f a flash object
#' @param var_type indicates what type of variance structure to assume
#' @return an updated flash object
#' @export
flash_update_precision = function(data,f,var_type=c("by_column","constant","by_row","kroneker")){
  if(data$S!=0){stop("not yet implemented")}

  R2 = get_R2(data,f)
  f$tau = compute_precision(R2,data$missing,var_type)
  return(f)
}

compute_precision = function(R2,missing,var_type=c("by_column","constant","by_row","kroneker")){
  R2[missing] = NA
  var_type= match.arg(var_type)
  if(var_type=="by_column"){tau = mle_precision_by_column(R2)}
  else if(var_type=="constant"){tau = mle_precision_constant(R2)}
  else(stop("that var_type not yet implemented"))
  tau[missing] = 0
  return(abs(tau))
}


#' @title  mle for precision (separate parameter for each column)
#' @param R2 n by p matrix of squared residuals (with NAs for missing)
#' @return n by p matrix of precisions (separate value for each column)
mle_precision_by_column = function(R2){
  sigma2 = colMeans(R2,na.rm =TRUE) # a p vector
  return(outer(rep(1,nrow(R2)), 1/sigma2)) # an n by p matrix
}

#' @title  mle for precision (separate parameter for each column)
#' @param R2 n by p matrix of squared residuals (with NAs for missing)
#' @return n by p matrix of precisions (separate value for each column)
mle_precision_constant = function(R2){
  sigma2 = mean(R2,na.rm =TRUE) # a scalar
  return(matrix(1/sigma2,nrow=nrow(R2),ncol=ncol(R2))) # an n by p matrix
}


