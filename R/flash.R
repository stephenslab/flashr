#' @title Fit the rank1 FLASH model to data
#' @param Y an n by p data matrix
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @return a fitted flash object
#' @export
flash_r1 = function(Y,init_method=c("svd","random"),tol=1e-2){
  init_method=match.arg(init_method)
  f = flash_init(Y,1,init_method)
  f = flash_optimize_single_fl(Y,f,1,tol)
  return(f)
}


#' @title Fit the FLASH model to data by a greedy approach
#' @details Fits the model by adding a factor and then optimizing it.
#' It is "greedy" in that it does not return to re-optimize previous factors.
#' Stops when an added factor contributes nothing, or Kmax is reached
#' @param Y an n by p data matrix
#' @param Kmax the maximum number of factors to consider
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @return a fitted flash object
#' @export
flash_greedy = function(Y,Kmax, init_method=c("svd","random"),tol=1e-2){
  init_method=match.arg(init_method)
  f = flash_r1(Y,init_method)
  for(k in 2:Kmax){
    f = flash_add_factor(Y, f, 1, init_method)
    message("fitting factor/loading ",k)
    f = flash_optimize_single_fl(Y,f,k,tol)
    if(is_tiny_fl(f,k)) #test whether the factor/loading combination is effectively 0
      break
  }
  return(f)
}


#' @title Refines a fit of the FLASH model to data by "backfitting"
#' @details Iterates through the factors of a flash object, updating each until convergence
#' @param Y an n by p data matrix
#' @param f a fitted flash object to be refined
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @return a fitted flash object
#' @export
flash_backfit = function(Y,f,tol=1e-2){
  c = get_conv_criteria(f)
  diff = 1
  while(diff > tol){
    for(k in 1:get_k(f)){
        f = flash_update_single_fl(Y,f,k)
    }
    cnew = get_conv_criteria(f)
    diff = mean((cnew-c)^2)
    c = cnew
  }
  return(f)
}


