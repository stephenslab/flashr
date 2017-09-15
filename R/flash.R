#' @title Fit the rank1 FLASH model to data
#' @param data a flash data object
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param init_method specifies how to initialize the factors. Options include svd on data matrix, or random (N(0,1))
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' data = set_flash_data(Y)
#' f = flash_r1(data)
#' flash_get_sizes(f)
#' @export
flash_r1 = function(data,init_method=c("svd","random"),tol=1e-2,ash_param=list()){
  init_method=match.arg(init_method)
  f = flash_init(data,1,init_method)
  f = flash_optimize_single_fl(data,f,1,tol)
  return(f)
}


#' @title Fit the FLASH model to data by a greedy approach
#' @details Fits the model by adding a factor and then optimizing it.
#' It is "greedy" in that it does not return to re-optimize previous factors.
#' Stops when an added factor contributes nothing, or Kmax is reached
#' @param data a flash data object
#' @param Kmax the maximum number of factors to consider
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' data = set_flash_data(Y)
#' f = flash_greedy(data,10)
#' flash_get_sizes(f)
#' @export
flash_greedy = function(data,Kmax, init_method=c("svd","random"),tol=1e-2,ash_param=list()){
  init_method=match.arg(init_method)
  f = flash_r1(data,init_method)
  for(k in 2:Kmax){
    f = flash_add_factor(data, f, 1, init_method)
    message("fitting factor/loading ",k)
    f = flash_optimize_single_fl(data,f,k,tol)
    if(is_tiny_fl(f,k)) #test whether the factor/loading combination is effectively 0
      break
  }
  return(f)
}


#' @title Refines a fit of the FLASH model to data by "backfitting"
#' @details Iterates through the factors of a flash object, updating each until convergence
#' @param data a flash data object
#' @param f a fitted flash object to be refined
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' data = set_flash_data(Y)
#' fg = flash_greedy(data,10)
#' fb = flash_backfit(data,fg) # refines fit from greedy by backfitting
#' flash_get_sizes(fb)
#' @export
flash_backfit = function(data,f,tol=1e-2,ash_param=list()){
  c = get_conv_criteria(f)
  diff = 1
  while(diff > tol){
    for(k in 1:get_k(f)){
        f = flash_update_single_fl(data,f,k)
    }
    cnew = get_conv_criteria(f)
    diff = mean((cnew-c)^2)
    c = cnew
  }
  return(f)
}


