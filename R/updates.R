#' @title  Update a flash loading
#' @details Updates loading k of f to increase the objective F
#' @param data a flash data object
#' @param f a flash object
#' @param k the index of the loading to update
#' @return an updated flash object
#' @export
flash_update_single_loading = function(data,f,k){
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0

  s = sqrt(1/(tau %*% f$EF2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(data,f,k) #residuals excluding factor k
    x = ((Rk*tau) %*% f$EF[,k]) * s^2
    a = ashr::ash(as.vector(x),as.vector(s),outputlevel=5,mixcompdist="normal",method="shrink")
    mat_post = list(comp_postprob=a$flash_data$comp_postprob,
                    comp_postmean=a$flash_data$comp_postmean,
                    comp_postmean2=a$flash_data$comp_postmean2)
    f$EL[,k] = a$flash_data$postmean
    f$EL2[,k] = a$flash_data$postmean2
    f$gl[[k]] = a$flash_data$fitted_g
    f$matl[[k]] = mat_post
  }
  return(f)
}

#' @title  Update a flash factor
#' @details Updates factor k of f to increase the objective F updates only the factor, once (not the loading)
#' @param data a flash data object
#' @param f a flash object
#' @param k the index of the factor to update
#' @return an updated flash object
#' @export
flash_update_single_factor = function(data,f,k){
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0

  s = sqrt(1/(t(tau) %*% f$EL2[,k]))
  if(sum(is.finite(s))>0){ # check some finite values before proceeding
    Rk = get_Rk(data,f,k) #residuals excluding factor k
    x = (t(Rk*tau) %*% f$EL[,k]) * s^2
    a = ashr::ash(as.vector(x),as.vector(s),outputlevel=5,mixcompdist="normal",method="shrink")
    mat_post = list(comp_postprob=a$flash_data$comp_postprob,
                    comp_postmean=a$flash_data$comp_postmean,
                    comp_postmean2=a$flash_data$comp_postmean2)
    f$EF[,k] = a$flash_data$postmean
    f$EF2[,k] = a$flash_data$postmean2
    f$gf[[k]] = a$flash_data$fitted_g
    f$matf[[k]] = mat_post
  }
  return(f)
}

#' @title Update a single flash factor-loading combination (and precision)
flash_update_single_fl = function(data,f,k){
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
#' @return an updated flash object
flash_optimize_single_fl = function(data,f,k,tol=1e-3){
  # c = get_conv_criteria(f)
  f = flash_update_single_fl(data,f,k)
  obj_val = flash_get_LB(data,f)
  obj_track = c(obj_val)
  diff = 1
  while(diff > tol){
    f = flash_update_single_fl(data,f,k)
    #cnew = get_conv_criteria(f)
    obj_val_new = flash_get_LB(data,f)
    #diff = mean((cnew-c)^2)
    diff =abs( obj_val_new - obj_val)
    obj_val = obj_val_new
    obj_track = c(obj_track, obj_val)
  }
  f$obj_track = obj_track
  return(f)
}


#' @title  Update precision parameter
#' @details Updates precision estimate to increase the objective F
#' @param data a flash data object
#' @param f a flash object
#' @return an updated flash object
#' @export
flash_update_precision = function(data,f){
  sigma2 = colMeans(get_R2(data,f),na.rm=TRUE)
  f$tau = outer(rep(1,get_n(f)), 1/sigma2)
  return(f)
}
