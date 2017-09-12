# functions for extracting useful information about the object

#' @title Return the estimated LF' matrix
#' @param f a flash object
#' @return  the estimated value of LF', an n by p matrix
#' @export
flash_get_lf = function(f){
  return(f$EL %*% t(f$EF))
}

#' @title  Get the residuals from a flash object, excluding factor k
get_Rk = function(Y,f,k){
  return( Y - f$EL[,-k,drop=FALSE] %*% t(f$EF[,-k,drop=FALSE]) )
}

#' @title  Get the residuals from data and a flash object
get_R = function(Y,f){
  return( Y - flash_get_lf(f) )
}


#' @title  Get the expected squared residuals from a flash object
get_R2 = function(Y,f){
  LF = f$EL %*% t(f$EF)
  return( Y^2 - 2 * Y * LF
          + LF^2 + f$EL2 %*% t(f$EF2) - f$EL^2 %*% t(f$EF^2) )
}

#' @details checks whether kth factor/loading combination is tiny
is_tiny_fl=function(f,k,tol=1e-8){
    return(sum(f$EL[,k]^2)*sum(f$EF[,k]^2) < tol)
}

get_l = function(f){f$EL}
get_f = function(f){f$EF}
get_k = function(f){ncol(f$EL)}
get_n = function(f){nrow(f$EL)}
get_p = function(f){nrow(f$EF)}

#' @details computes the factor contributions (analogous to eigenvalues)
#' for each factor in f
#' @export
flash_get_sizes = function(f){
  colSums(f$EL^2)*colSums(f$EF^2)
}

get_conv_criteria = function(f){flash_get_lf(f)}
