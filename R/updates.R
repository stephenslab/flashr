#' @title  Update a flash loading
#' @details Updates loading k of f to increase the objective F
#' @param Y an n times p data matrix
#' @param f a flash object
#' @param k the index of the loading to update
#' @return an updated flash object
#' @export
flash_update_single_loading = function(Y,f,k){
  s = sqrt(1/(f$tau %*% f$EF2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(Y,f,k) #residuals excluding factor k
    x = ((Rk*f$tau) %*% f$EF[,k]) * s^2
    a = ashr::ash(as.vector(x),as.vector(s),outputlevel=4,method="shrink")
    f$EL[,k] = ashr::get_pm(a)
    f$EL2[,k] = ashr::get_psd(a)^2 + ashr::get_pm(a)^2
    f$gl[[k]] = ashr::get_fitted_g(a)
  }
  return(f)
}

#' @title  Update a flash factor
#' @details Updates factor k of f to increase the objective F
#' @param Y an n times p data matrix
#' @param f a flash object
#' @param k the index of the factor to update
#' @return an updated flash object
#' @export
flash_update_single_factor = function(Y,f,k){
  s = sqrt(1/(t(f$tau) %*% f$EL2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(Y,f,k) #residuals excluding factor k
    x = (t(Rk*f$tau) %*% f$EL[,k]) * s^2
    a = ashr::ash(as.vector(x),as.vector(s),outputlevel=4,method="shrink")
    f$EF[,k] = ashr::get_pm(a)
    f$EF2[,k] = ashr::get_psd(a)^2 + ashr::get_pm(a)^2
    f$gf[[k]] = ashr::get_fitted_g(a)
  }
  return(f)
}

#' @title  Update precision parameter
#' @details Updates precision estimate to increase the objective F
#' @param Y an n times p data matrix
#' @param f a flash object
#' @return an updated flash object
#' @export
flash_update_precision = function(Y,f){
  sigma2 = colMeans(get_R2(Y,f),na.rm=TRUE)
  f$tau = outer(rep(1,get_n(f)), 1/sigma2)
  return(f)
}
