#' @title  Update a flash loading
#' @details Updates loading k of f to increase the objective F
#' Updates only the loading, once (not the factor)
#' @param data a flash data object
#' @param f a flash fit object
#' @param k the index of the loading to update
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return an updated flash object
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
    f$ash_param_l[[k]] = ash_param
    f$KL_l[[k]] = a$flash_data$penloglik -
      NM_posterior_e_loglik(x,s,a$flash_data$postmean,a$flash_data$postmean2)
  }
  return(f)
}

#' @title  Update a flash factor
#' @details Updates factor k of f to increase the objective F.
#' Updates only the factor, once (not the loading)
#' @inheritParams flash_update_single_loading
#' @return an updated flash object
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
    f$ash_param_f[[k]] = ash_param
    f$KL_f[[k]] = a$flash_data$penloglik -
      NM_posterior_e_loglik(x,s,a$flash_data$postmean,a$flash_data$postmean2)
  }
  return(f)
}

#' @title Update a single flash factor-loading combination (and precision)
#' @inheritParams flash_update_single_loading
flash_update_single_fl = function(data,f,k,var_type,ash_param=list()){
  f = flash_update_precision(data,f,var_type)
  f = flash_update_single_factor(data,f,k,ash_param)
  f = flash_update_single_loading(data,f,k,ash_param)
  return(f)
}

#' @title  Optimize a flash factor-loading combination
#' @details Iteratively updates factor and loading k of f (as well as residual precision)
#' to convergence of objective (used in the greedy algorithm for example)
#' @param data a flash data object
#' @param f a flash object
#' @param k the index of the factor/loading to optimize
#' @param var_type type of variance structure to assume for residuals.
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @param tol a tolerance for the optimization
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return an updated flash object
flash_optimize_single_fl = function(data,f,k,var_type,nullcheck=TRUE,tol=1e-2,ash_param=list(),verbose=FALSE){
  c = get_conv_criteria(data,f)
  diff = 1
  while(diff > tol){
    f = flash_update_single_fl(data,f,k,var_type,ash_param)
    cnew = get_conv_criteria(data,f)
    diff = sqrt(mean((cnew-c)^2))
    c = cnew
    if(verbose){
      message("objective: ",c)
    }
  }

  if(nullcheck){
    f = perform_nullcheck(data,f,k,var_type,verbose)
  }

  return(f)
}

#' @title  Zeros out factors when that improves the objective
#' @details Sometimes zeroing out a factor can improve the objective.
#' This function iterates over factors with indices in kset
#' and checks whether zeroing it out will improve the objective; if so then that factor
#' is set to 0 (and precision is updated). Returns the final flash fit object obtained when this iterative process stops
#' (ie a complete pass is performed with no factor being zerod)
#' @param data a flash data object
#' @param f a flash object
#' @param kset the indices of the factor/loading to check
#' @param var_type type of variance structure to assume for residuals.
#' @param verbose if TRUE various output progress updates will be printed
#' @return a flash object
perform_nullcheck=function(data,f,kset,var_type,verbose){

  f_changed = TRUE #we are going to iterate until f does not change
  while(f_changed){

    f_changed = FALSE
    for(k in kset){

      f0 = flash_zero_out_factor(data,f,k)
      f0 = flash_update_precision(data,f0,var_type)
      F0 = get_F(data,f0)
      F1 = get_F(data,f)

      if(verbose){
        message("performing nullcheck")
        message("objective from deleting factor:",F0)
        message("objective from keeping factor:",F1)
      }

      if(F0>F1){
        if(verbose){
          message("factor zeroed out")
        }
        f=f0
        f_changed = TRUE
      }

    }
  }

  return(f)
}


