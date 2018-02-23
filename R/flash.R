#' @title Fit the rank1 flash model to data
#' @param data an n by p matrix or a flash data object created using \code{flash_set_data}
#' @param f_init if supplied, a flash object to which a single new factor is to be added
#' @param var_type type of variance structure to assume for residuals.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param init_fn function to be used to initialize the factor. This function should take parameters (Y,K)
#' where Y is an n by p matrix of data (or a flash data object) and K is a number of factors.
#' It should output a list with elements (u,d,v) where u is n by K matrix
#' v is a p by K matrix  and d is a K vector. See \code{udv_si} for an example.
#' (If the input data includes missing values then this function must be able
#' to deal with missing values in its input matrix.)
#' @param ebnm_fn function to solve the Empirical Bayes Normal Means problem
#' @param ebnm_param named list containing parameters to be passed to ebnm_fn when optimizing; defaults set by flash_default_ebnm_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @param seed a random number seed to use before running method - for reproducibility. Set to NULL if you don't want seed set.
#' (The seed can affect initialization when there are missing data; otherwise the algorithm is deterministic)
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' f = flash_r1(Y)
#' flash_get_sizes(f)
#' f2 = flash_r1(Y,ebnm_fn=ebnm_pn) # run with the faster ebnm function (uses point-normal prior)
#' # example to show how to pass parameters to ash
#' f3 = flash_r1(Y,ebnm_fn = ebnm_ash, ebnm_param = list(mixcompdist = "normal",method="fdr"))
#' @export
flash_r1 = function(data,f_init=NULL,var_type = c("by_column","constant"), init_fn = "udv_si",tol=1e-2,ebnm_fn = ebnm_ash, ebnm_param=flash_default_ebnm_param(ebnm_fn),verbose = FALSE, nullcheck=TRUE,seed=123){
  if(!is.null(seed)){set.seed(seed)}
  if(is.matrix(data)){data = flash_set_data(data)}
  var_type=match.arg(var_type)
  f = flash_add_factors_from_data(data,f_init = f_init, init_fn=init_fn,K=1)
  f = flash_optimize_single_fl(data,f,get_k(f),var_type,nullcheck,tol,ebnm_fn,ebnm_param,verbose)
  return(f)
}


#' @title Adds factors to a flash object by a greedy approach
#' @details Adds factors iteratively, at each time adding a new factor and then optimizing it.
#' It is "greedy" in that it does not return to re-optimize previous factors.
#' The function stops when an added factor contributes nothing, or Kmax is reached.
#' Each new factor is intialized by applying the function `init_fn` to the residuals
#' after removing previously-fitted factors.
#' @param data an n by p matrix or a flash data object created using \code{flash_set_data}
#' @param Kmax the maximum number of factors to add to f_init
#' @param f_init a flash fit object to start the greedy algorithm: the greedy algorithm iteratively adds up to Kmax factors
#' to this initial fit. (If NULL then the greedy algorithm starts with 0 factors)
#' @param var_type type of variance structure to assume for residuals.
#' @param init_fn function to be used to initialize each factor when added. This function should take as
#' input an n by p matrix of data (or a flash data object)
#' and output a list with elements (u,d,v) where u is an n-vector,
#' v is a p-vector and d is a scalar. See \code{udv_si} for an example,
#' and examples below. (If the input data includes missing values then this function must be able
#' to deal with missing values in its input matrix.)
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ebnm_fn function to solve the Empirical Bayes Normal Means problem
#' @param ebnm_param parameters to be passed to ebnm_fn when optimizing; defaults set by flash_default_ebnm_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @param nullcheck flag whether to check, after running
#' hill-climbing updates, whether the achieved optimum is better than setting factor to 0.
#' If this check is performed and fails then the factor will be set to 0 in the returned fit.
#' @param seed a seed for the random number to be set before running, for reproducibility. Set to NULL if you don't want seed set.
#' (The seed can affect initialization when there are missing data; otherwise the algorithm is deterministic)
#' @return a fitted flash object
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l,f) + matrix(rnorm(1000),nrow=100)
#' f = flash_add_greedy(Y,10)
#' flash_get_sizes(f)
#' # example to show how to use a different initialization function
#' f2 = flash_add_greedy(Y,10,function(x,K=1){softImpute::softImpute(x,K,lambda=10)})
#' @export
flash_add_greedy = function(data,Kmax=1,f_init = NULL,var_type = c("by_column","constant"), init_fn="udv_si",tol=1e-2,ebnm_fn = ebnm_ash, ebnm_param=flash_default_ebnm_param(ebnm_fn),verbose=FALSE,nullcheck=TRUE,seed=123){
  if(!is.null(seed)){set.seed(seed)}
  if(is.matrix(data)){data = flash_set_data(data)}
  var_type=match.arg(var_type)
  f = f_init

  for(k in 1:Kmax){
    message("fitting factor/loading ",k)
    f = flash_r1(data,f,var_type,init_fn,tol,ebnm_fn,ebnm_param,verbose,nullcheck,seed=NULL)
    if(is_tiny_fl(f,get_k(f))) #test whether the factor/loading combination is effectively 0
      break
  }

  return(f)
}


#' @title Refines a fit of the flash model to data by "backfitting"
#' @details Iterates through the factors of a flash object, updating each until convergence
#' @param data an n by p matrix or a flash data object created using \code{flash_set_data}
#' @param f a fitted flash object to be refined
#' @param kset the indices of factors to be optimized (NULL indicates all factors)
#' @param var_type type of variance structure to assume for residuals.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ebnm_fn function to solve the Empirical Bayes Normal Means problem
#' @param ebnm_param parameters to be passed to ebnm_fn when optimizing; defaults set by flash_default_ebnm_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @param maxiter a maximum number of iterations to perform (in the inner loop). To perform just one iteration we suggest setting maxiter=1 and nullcheck=FALSE
#' @return a fitted flash object
#' @examples
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' fg = flash_greedy(Y,10)
#' fb = flash_backfit(Y,fg) # refines fit from greedy by backfitting
#' flash_get_sizes(fb)
#' fsi = flash_init_fn(flash_set_data(Y),"udv_si",4)
#' fb2 = flash_backfit(Y,fsi)
#' flash_get_sizes(fb2)
#' @export
flash_backfit = function(data,f,kset=NULL,var_type = c("by_column","constant"),tol=1e-2,ebnm_fn = ebnm_ash, ebnm_param=flash_default_ebnm_param(ebnm_fn),verbose=FALSE,nullcheck=TRUE,maxiter = 1000){
  if(is.matrix(data)){data = flash_set_data(data)}
  if(is.null(kset)){kset = 1:get_k(f)}
  var_type=match.arg(var_type)

  if(verbose){message("iteration:1")}
  for(k in kset){
    f = flash_update_single_fl(data,f,k,var_type,ebnm_fn,ebnm_param)
  }

  c = flash_get_objective(data,f)
  if(verbose){message("objective: ",c)}

  diff = Inf
  iteration = 2

  while((diff > tol) & (iteration <= maxiter)){
    # There are two steps; first backfit, then null check
    # (if nullcheck removes some factors then the whole process is repeated)
    while((diff > tol) & (iteration <= maxiter)){
      if(verbose){message("iteration:", iteration)}
      for(k in kset){
        f = flash_update_single_fl(data,f,k,var_type,ebnm_fn,ebnm_param)
      }
      cnew = flash_get_objective(data, f)
      diff = cnew-c
      c = cnew
      if(verbose){message("objective: ",c)}
      iteration = iteration + 1
    }

    if(nullcheck){
      #remove factors that actually hurt objective
      kset = 1:get_k(f)
      f = perform_nullcheck(data,f,kset,var_type,verbose)

      # recompute objective; if it changes then whole process will be repeated
      cnew = flash_get_objective(data, f)
      diff = cnew-c
      c = cnew
      iteration = 1
    }

  }

  return(f)
}

#' @title Main flash function
#' @details Performs Empirical Bayes factor analysis with adaptive shrinkage on both factors and loadings.
#' @param data an n by p matrix or a flash data object created using \code{flash_set_data}
#' @param Kmax the maximum total number of factors to use
#' @param column_covariates an n by r1 matrix of covariates (eg could be a column of all 1s to allow an intercept for each column)
#' @param row_covariates a p by r2 matrix of covariates (eg could be a vector of p 1s to allow an intercept for each row)
#' @param init_fn function used to initialize factors and loadings,
#' @param f_init a previously fitted flash object to be added to or refined
#' @param var_type type of variance structure to assume for residuals.
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @param ebnm_fn function to solve the Empirical Bayes Normal Means problem
#' @param ebnm_param parameters to be passed to ebnm_fn when optimizing; defaults set by flash_default_ebnm_param()
#' @param verbose if TRUE various output progress updates will be printed
#' @return a fitted flash object
#' @export
flash = function(data,Kmax,column_covariates = NULL,row_covariates = NULL,f_init = NULL,init_fn=NULL, var_type = c("by_column","constant"),tol=1e-2,ebnm_fn = ebnm_ash, ebnm_param=flash_default_ebnm_param(ebnm_fn),verbose=FALSE){
  if(is.matrix(data)){data = flash_set_data(data)}
  var_type=match.arg(var_type)
  f=f_init
  if(is.null(f)){f = flash_init_null()}


  if(!is.null(column_covariates)){
    f = flash_combine(f, flash_init_from_L(column_covariates,fixl))
  }
  # this not yet written

  # need to have initialization from L or F alone, with ability to fix
  # also need to think about how to initialize factors from a given sparsity pattern
  if(!is.null(row_covariates)){
    f = flash_combine(f, flash_init_from_F(row_covariates,fixf))
  }

  return(f)
}


