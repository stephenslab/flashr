#' @title  Update a flash loading
#' @details Updates loading k of f to increase the objective F
#' Updates only the loading, once (not the factor)
#' @param data a flash data object
#' @param f a flash fit object
#' @param k the index of the loading to update
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return an updated flash object
#' @export
flash_update_single_loading = function(data,f,k,ash_param=list()){
  ash_param=modifyList(flash_default_ash_param(),ash_param)
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0

  s = sqrt(1/(tau %*% f$EF2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(data,f,k) #residuals excluding factor k
    x = ((Rk*tau) %*% f$EF[,k]) * s^2
    a = do.call(ashr::ash,
                c( list(betahat=as.vector(x),sebetahat=as.vector(s)),
                   ash_param) )
    f$EL[,k] = a$flash_data$postmean
    f$EL2[,k] = a$flash_data$postmean2
    f$gl[[k]] = a$flash_data$fitted_g
  }
  return(f)
}

#' @title  Update a flash factor
#' @details Updates factor k of f to increase the objective F.
#' Updates only the factor, once (not the loading)
#' @inheritParams flash_update_single_loading
#' @return an updated flash object
#' @export
flash_update_single_factor = function(data,f,k,ash_param=list()){
  ash_param=modifyList(flash_default_ash_param(),ash_param)
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0

  s = sqrt(1/(t(tau) %*% f$EL2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(data,f,k) #residuals excluding factor k
    x = (t(Rk*tau) %*% f$EL[,k]) * s^2
    a = do.call(ashr::ash,
                c( list(betahat=as.vector(x),sebetahat=as.vector(s)),
                  ash_param) )

    f$EF[,k] = a$flash_data$postmean
    f$EF2[,k] = a$flash_data$postmean2
    f$gf[[k]] = a$flash_data$fitted_g
  }
  return(f)
}

#' @title Update a single flash factor-loading combination (and precision)
#' @inheritParams flash_update_single_loading
flash_update_single_fl = function(data,f,k,ash_param=list()){
  f = flash_update_precision(data,f)
  f = flash_update_single_factor(data,f,k)
  f = flash_update_single_loading(data,f,k)
  return(f)
}

#' @title  Optimize a flash factor-loading combination
#' @details Iteratively updates factor and loading k of f (as well as residual precision)
#' to convergence of objective (used in the greedy algorithm for example)
#' @param data a flash data object
#' @param f a flash object
#' @param k the index of the factor/loading to optmize
#' @param tol a tolerance for the optimization
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return an updated flash object
flash_optimize_single_fl = function(data,f,k,tol=1e-2,ash_param=list()){
  c = get_conv_criteria(f)
  diff = 1
  while(diff > tol){
    f = flash_update_single_fl(data,f,k)
    cnew = get_conv_criteria(f)
    diff = sqrt(mean((cnew-c)^2))
    c = cnew
  }
  return(f)
}


#' @title  Update precision parameter
#' @details Updates precision estimate to increase the objective F
#' @param data a flash data object
#' @param f a flash object
#' @return an updated flash object
#' @export
flash_update_precision = function(data,f){
  if(data$anyNA){
    sigma2 = colSums(get_R2(data,f) * !(data$missing)) / colSums(!data$missing)
  } else {
    sigma2 = colMeans(get_R2(data,f))
  }
  f$tau = outer(rep(1,get_n(f)), 1/sigma2)
  return(f)
}



