#' @title Fit the rank1 FLASH model to data
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param var_type type of variance structure to assume for residuals.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param init_fn function to be used to initialize the factor. This function should take parameters (Y,K)
#' where Y is an n by p matrix of data (or a flash data object) and K is a number of factors.
#' It should output a list with elements (u,d,v) where u is n by K matrix
#' v is a p by K matrix  and d is a K vector. See \code{udv_si} for an example.
#' (If the input data includes missing values then this function must be able
#' to deal with missing values in its input matrix.)
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' f = flash_r1(Y)
#' flash_get_sizes(f)
#' @export
flash_r1 = function(data,var_type = c("by_column","constant"), init_fn = "udv_si",tol=1e-2,ash_param=list(),verbose = FALSE, nullcheck=TRUE){
  if(is.matrix(data)){data = set_flash_data(data)}
  var_type=match.arg(var_type)
  f = flash_init_fn(data,init_fn)
  f = flash_optimize_single_fl(data,f,1,var_type,nullcheck,tol,ash_param,verbose)
  return(f)
}


#' @title Fit the FLASH model to data by a greedy approach
#' @details Fits the model by adding a factor and then optimizing it.
#' It is "greedy" in that it does not return to re-optimize previous factors.
#' The function stops when an added factor contributes nothing, or Kmax is reached.
#' Each new factor is intialized by applying the function `init_fn` to the residuals
#' after removing previously-fitted factors.
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param Kmax the maximum number of factors to consider
#' @param var_type type of variance structure to assume for residuals.
#' @param init_fn function to be used to initialize each factor when added. This function should take as
#' input an n by p matrix of data (or a flash data object)
#' and output a list with elements (u,d,v) where u is an n-vector,
#' v is a p-vector and d is a scalar. See \code{udv_si} for an example,
#' and examples below. (If the input data includes missing values then this function must be able
#' to deal with missing values in its input matrix.)
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @return a fitted flash object
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l,f) + matrix(rnorm(1000),nrow=100)
#' f = flash_greedy(Y,10)
#' flash_get_sizes(f)
#' # example to show how to use a different initialization function
#' f2 = flash_greedy(Y,10,function(x,K=1){softImpute::softImpute(x,K,lambda=10)})
#' @export
flash_greedy = function(data,Kmax=1,var_type = c("by_column","constant"),init_fn="udv_si",tol=1e-2,ash_param=list(),verbose=FALSE,nullcheck=TRUE){
  if(is.matrix(data)){data = set_flash_data(data)}
  var_type=match.arg(var_type)
  message("fitting factor/loading ",1)
  f = flash_r1(data,init_fn,tol,ash_param,verbose,nullcheck)
  if(Kmax>1 & !is_tiny_fl(f,1)){ #if the first factor is big enough
    for(k in 2:Kmax){
      f = flash_add_factor(data, f, init_fn)
      message("fitting factor/loading ",k)
      f = flash_optimize_single_fl(data,f,k,var_type,nullcheck,tol,ash_param,verbose)
      if(is_tiny_fl(f,k)) #test whether the factor/loading combination is effectively 0
        break
    }
  }
  return(f)
}


#' @title Refines a fit of the FLASH model to data by "backfitting"
#' @details Iterates through the factors of a flash object, updating each until convergence
#' @param data an n by p matrix or a flash data object created using \code{set_flash_data}
#' @param f a fitted flash object to be refined
#' @param var_type type of variance structure to assume for residuals.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ash_param parameters to be passed to ashr when optimizing; defaults set by flash_default_ash_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' fg = flash_greedy(Y,10)
#' fb = flash_backfit(Y,fg) # refines fit from greedy by backfitting
#' flash_get_sizes(fb)
#' fsi = flash_init_fn(data,"udv_si",5)
#' fb2 = flash_backfit(Y,fsi)
#' flash_get_sizes(fb2)
#' @export
flash_backfit = function(data,f,var_type = c("by_column","constant"),tol=1e-2,ash_param=list(),verbose=FALSE){
  if(is.matrix(data)){data = set_flash_data(data)}
  var_type=match.arg(var_type)
  c = get_conv_criteria(data, f)
  diff = 1
  while(diff > tol){
    for(k in 1:get_k(f)){
        f = flash_update_single_fl(data,f,k,var_type,ash_param)
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


