# functions for extracting useful information about the object

#' @title Return the estimated LF' matrix
#' @param f a flash fit object
#' @return  the estimated value of LF', an n by p matrix
#' @export
flash_get_lf = function(f){
  return(f$EL %*% t(f$EF))
}

#' @title  Get the residuals from a flash data and fit object,
#' excluding factor k
get_Rk = function(data,f,k){
  return( data$Y - f$EL[,-k,drop=FALSE] %*% t(f$EF[,-k,drop=FALSE]) )
}

#' @title  Get the residuals from data and a flash fit object
get_R = function(data,f){
  return( data$Y - flash_get_lf(f) )
}


#' @title  Get the expected squared residuals from a flash data and fit object
get_R2 = function(data,f){
  LF = f$EL %*% t(f$EF)
  return( data$Y^2 - 2 * data$Y * LF
          + LF^2 + f$EL2 %*% t(f$EF2) - f$EL^2 %*% t(f$EF^2) )
}

#' @title is_tiny_fl
#' @details checks whether kth factor/loading combination is tiny
is_tiny_fl=function(f,k,tol=1e-8){
    return(sum(f$EL[,k]^2)*sum(f$EF[,k]^2) < tol)
}

get_l = function(f){f$EL}
get_f = function(f){f$EF}
get_k = function(f){ncol(f$EL)}
get_n = function(f){nrow(f$EL)}
get_p = function(f){nrow(f$EF)}

#' @title flash_get_sizes
#' @details computes the factor contributions (analogous to eigenvalues)
#' for each factor in flash fit f
#' @export
flash_get_sizes = function(f){
  colSums(f$EL^2)*colSums(f$EF^2)
}

#get_conv_criteria = function(f){flash_get_lf(f)}
get_conv_criteria = function(data,f){
  get_F(data,f)
}
