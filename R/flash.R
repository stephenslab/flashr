#' @title Fit the rank1 FLASH model to data
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param init_method specifies how to initialize the factors. Options include svd on data matrix, or random (N(0,1))
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' f = flash_r1(Y)
#' flash_get_sizes(f)
#' @export
flash_r1 = function(data,init_method=c("softImpute","svd","random"),nullcheck=TRUE,tol=1e-2,ash_param=list(),verbose = FALSE){
  if(is.matrix(data)){data = set_flash_data(data)}
  init_method=match.arg(init_method)
  f = flash_init(data,1,init_method)
  f = flash_optimize_single_fl(data,f,1,nullcheck,tol,ash_param,verbose)
  return(f)
}


#' @title Fit the FLASH model to data by a greedy approach
#' @details Fits the model by adding a factor and then optimizing it.
#' It is "greedy" in that it does not return to re-optimize previous factors.
#' Stops when an added factor contributes nothing, or Kmax is reached
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param Kmax the maximum number of factors to consider
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' f = flash_greedy(Y,10)
#' flash_get_sizes(f)
#' @export
flash_greedy = function(data,Kmax, nullcheck=TRUE,init_method=c("softImpute","svd","random"),tol=1e-2,ash_param=list(),verbose=FALSE){
  if(is.matrix(data)){data = set_flash_data(data)}
  init_method=match.arg(init_method)
  f = flash_r1(data,init_method,nullcheck)
  for(k in 2:Kmax){
    f = flash_add_factor(data, f, 1, init_method)
    message("fitting factor/loading ",k)
    f = flash_optimize_single_fl(data,f,k,nullcheck,tol,ash_param,verbose)
    if(is_tiny_fl(f,k)) #test whether the factor/loading combination is effectively 0
      break
  }
  return(f)
}


#' @title Refines a fit of the FLASH model to data by "backfitting"
#' @details Iterates through the factors of a flash object, updating each until convergence
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param f a fitted flash object to be refined
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' fg = flash_greedy(Y,10)
#' fb = flash_backfit(Y,fg) # refines fit from greedy by backfitting
#' flash_get_sizes(fb)
#' @export
flash_backfit = function(data,f,tol=1e-2,ash_param=list(),verbose=FALSE){
  if(is.matrix(data)){data = set_flash_data(data)}
  c = get_conv_criteria(data, f)
  diff = 1
  while(diff > tol){
    for(k in 1:get_k(f)){
        f = flash_update_single_fl(data,f,k,ash_param)
    }
    cnew = get_conv_criteria(data, f)
    diff = sqrt(mean((cnew-c)^2))
    c = cnew
    if(verbose){
      message("objective: ",c)
    }
  }
  return(f)
}


