# functions for extracting useful information about the object

#' @title Return the estimated LF' matrix
#' @param f a flash fit object
#' @return  the estimated value of LF', an n by p matrix
#' @export
flash_get_lf = function(f){
  if(is.null(f$EL)){return(NULL)}
  return(f$EL %*% t(f$EF))
}

#' @title  Get the residuals from a flash data and fit object,
#' excluding factor k
get_Rk = function(data,f,k){
  if(get_k(f)<k){stop("factor k does not exist")}
  return( data$Y - f$EL[,-k,drop=FALSE] %*% t(f$EF[,-k,drop=FALSE]) )
}

#' @title  Get the residuals from data and a flash fit object
get_R = function(data,f){
  if(is.null(f$EL)){return(data$Y)} # if f is null, return Y
  else{
    return( data$Y - flash_get_lf(f) )
  }
}

#' @title  Get the residuals from data and a flash fit object, with missing data as in original
get_R_withmissing = function(data,f){
  if(is.null(f$EL)){return(get_Yorig(data))} # if f is null, return Y
  else{
    return( get_Yorig(data) - flash_get_lf(f) )
  }
}

#' @title  Get the expected squared residuals from a flash data and fit object
get_R2 = function(data,f){
  if(is.null(f$EL)){return(data$Y^2)}
  else{
    LF = f$EL %*% t(f$EF)
    return( (data$Y - LF)^2 +
          f$EL2 %*% t(f$EF2) - f$EL^2 %*% t(f$EF^2) )
  }
}

#' @title is_tiny_fl
#' @details checks whether kth factor/loading combination is tiny
is_tiny_fl=function(f,k,tol=1e-8){
    return(sum(f$EL[,k]^2)*sum(f$EF[,k]^2) < tol)
}

#' @title flash_get_l
#' @param f a flash fit object
#' @param k indices of loadings to be returned
#' @details returns the loadings from a flash fit (all by default)
#' @export
flash_get_l = function(f,k=NULL){
  if(is.null(k)){k=1:get_k(f)}
  f$EL[,k]}

#' @title flash_get_f
#' @param f a flash fit object
#' @param k indices of factors to be returned
#' @details returns the factors from a flash fit (all by default)
#' @export
flash_get_f = function(f,k=NULL){
  if(is.null(k)){k=1:get_k(f)}
  f$EF[,k]}

#' @title Get number of factors in a fit object
#' @param f a flash fit object
#' @details returns the number of factors in a flash fit
#' @export
flash_get_k = function(f){get_k(f)}


get_l = function(f){f$EL}
get_f = function(f){f$EF}

get_k = function(f){
  k=ncol(f$EL)
  if(is.null(k)){return(0)}
  else{return(k)}
}

get_n = function(f){nrow(f$EL)}
get_p = function(f){nrow(f$EF)}

#' @title flash_get_sizes
#' @details computes the factor contributions (analogous to eigenvalues)
#' for each factor in flash fit f
#' @export
flash_get_sizes = function(f){
  colSums(f$EL^2)*colSums(f$EF^2)
}

#' @title flash_get_pve
#' @details returns the factor contributions ("proportion of variance explained")
#' by each factor/loading combination in flash fit f. Because the factors are
#' not required to be orthogonal this should be interpreted loosely: eg PVE could total more than 1.
#' @export
flash_get_pve = function(f){
  s = flash_get_sizes(f)
  s/(sum(s)+sum(1/f$tau))
# wei's version  (sapply(seq(1,K),function(x){ sum(f$EL2[,x] %*% t(f$EF2[,x])) }))/sum(Y^2)
}

#' @title flash_get_udv
#' @details returns the normalized factors and weights in the same form as svd
#' That is, as a list, with elements u, d and v, where the columns of u and v are normalized to have length 1
#' @export
flash_get_udv = function(f){
  u = scale(f$EL,center=FALSE,scale=sqrt(colSums(f$EL^2)))
  v = scale(f$EF,center=FALSE,scale=sqrt(colSums(f$EF^2)))
  d = colSums(f$EL^2)*colSums(f$EF^2)
  return(list(u=u,d=d,v=v))
}


#get_conv_criteria = function(f){flash_get_lf(f)}
get_conv_criteria = function(data,f){
  flash_get_F(data,f)
}
