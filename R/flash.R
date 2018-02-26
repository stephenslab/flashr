#' @title Fit Empirical Bayes Matrix Factorization
#' 
#' @description This is the main interface for fitting EBMF models
#'   based on algorithms from Wang and Stephens.  The default behaviour
#'   is simply to run the greedy algorithm and return the result.  To
#'   follow it by backfitting set \code{backfit = TRUE}.
#' 
#' @param data An n by p matrix or a flash data object created using
#'   \code{flash_set_data}.
#' 
#' @param f_init if supplied, a flash object to which a single new
#'   factor is to be added.
#' 
#' @param var_type type of variance structure to assume for residuals.
#' 
#' @param init_fn function to be used to initialize the factor. This
#'   function should take parameters (Y,K) where Y is an n by p matrix
#'   of data (or a flash data object) and K is a number of factors.  It
#'   should output a list with elements (u,d,v) where u is n by K matrix
#'   v is a p by K matrix and d is a K vector. See \code{udv_si} for an
#'   example.  (If the input data includes missing values then this
#'   function must be able to deal with missing values in its input
#'   matrix.)
#' 
#' @param tol Specify how much objective can change in a single
#'   iteration to be considered not converged.
#' 
#' @param ebnm_fn function to solve the Empirical Bayes Normal Means problem
#' 
#' @param ebnm_param named list containing parameters to be passed to
#'   ebnm_fn when optimizing; defaults set by flash_default_ebnm_param()
#' 
#' @param verbose if TRUE various output progress updates will be printed
#' 
#' @param nullcheck flag whether to check, after running hill-climbing
#'   updates, whether the achieved optimum is better than setting factor
#'   to 0. If this check is performed and fails then the factor will be
#'   set to 0 in the returned fit.
#' 
#' @param seed a random number seed to use before running method - for
#'   reproducibility. Set to NULL if you don't want seed set.  (The seed
#'   can affect initialization when there are missing data; otherwise
#'   the algorithm is deterministic)
#' 
#' @param greedy a TRUE/FALSE flag to say whether to start by running
#'   the greedy algorithm (if FALSE then f_init must be supplied)
#' 
#' @param backfit a TRUE/FALSE flag to say whether to run backfitting.
#' 
#' @seealso flash_add_greedy, flash_backfit
#' 
#' @return A fitted flash object. Use \code{flash_get_ldf} to access
#'   standardized loadings and factors; use \code{flash_get_lf} to
#'   access fitted LF'.
#' 
#' @examples
#' 
#' set.seed(1) # for reproducibility
#' ftrue = matrix(rnorm(200),ncol=2)
#' ltrue = matrix(rnorm(40),ncol=2)
#' ltrue[1:10,1] = 0 # set up some sparsity
#' ltrue[11:20,2] = 0
#' Y = ltrue %*% t(ftrue)+rnorm(2000) # set up a simulated matrix
#' f = flash(Y)
#' ldf = flash_get_ldf(f)
#' ldf$d  #show the weights, analogous to singular values showing importance of each factor
#' plot(ltrue[,1],ldf$l[,2]) # plot true l against estimated l; with this seed it turns out the 2nd loading/factor corresponds to the first column of ltrue
#' plot(ftrue[,1],ldf$f[,2]) # plot true f against estimated f (note estimate is normalized)
#' plot(ltrue %*% t(ftrue), flash_get_lf(f)) #plot true lf' against estimated lf'; the scale of the estimate matches the data
#'
#' # Example to use the more flexible ebnm function in ashr.
#' f2 = flash(Y,ebnm_fn = ebnm_ash)
#' 
#' # example to show how to pass parameters to ashr (may be most useful for research use)
#' f3= flash(Y,ebnm_fn = ebnm_ash, ebnm_param = list(mixcompdist = "normal",method="fdr"))
#' 
#' # example to show how to use a different initialization function
#' f4 = flash(Y,init_fn = function(x,K=1){softImpute::softImpute(x,K,lambda=10)})
#' 
#' @export
#' 
flash = function(data,
                 Kmax = 100,
                 f_init = NULL,
                 var_type = c("by_column", "constant"),
                 init_fn = "udv_si",
                 tol = 1e-2,
                 ebnm_fn = ebnm_pn,
                 ebnm_param = flash_default_ebnm_param(ebnm_fn),
                 verbose = FALSE,
                 nullcheck = TRUE,
                 seed = 123,
                 greedy = TRUE,
                 backfit = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if(!greedy & is.null(f_init)){
    stop("If greedy is false then must provide f_init")
  }
  if(!greedy & !backfit){
    warning("If both greedy and backfit are false then nothing to do!")
  }

  if (greedy) {
    f = flash_add_greedy(
      data,
      Kmax,
      f_init,
      var_type,
      init_fn,
      tol,
      ebnm_fn,
      ebnm_param,
      verbose,
      nullcheck,
      seed
    )
  } else {
    f = f_init
  }

  if (backfit) {
    f = flash_backfit(data, f, kset=NULL, var_type, tol, ebnm_fn,
                      ebnm_param, verbose, nullcheck)
  }
  return(f)
}

#' @title Fit Empirical Bayes Matrix Factorization (greedy algorithm)
#' 
#' @description This implements the greedy algorithm from Wang and
#'   Stephens. It can be used to adds factors to an existing fit, or
#'   start from scratch.  It adds factors iteratively, at each stage
#'   adding a new factor and then optimizing it.  It is "greedy" in that
#'   it does not return to re-optimize previous factors.  The function
#'   stops when an added factor contributes nothing, or Kmax is reached.
#'   Each new factor is intialized by applying the function `init_fn` to
#'   the residuals after removing previously-fitted factors.
#' 
#' @inheritParams flash_r1
#' 
#' @return A fitted flash object.
#' 
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l,f) + matrix(rnorm(1000),nrow=100)
#' f = flash_add_greedy(Y,10)
#' flash_get_ldf(f)$d #gives the weights for each factor (analogue of singular values)
#' # example to show how to use a different initialization function
#' f2 = flash_add_greedy(Y,10,init_fn = function(x,K=1){softImpute::softImpute(x,K,lambda=10)})
#' 
#' @export
#' 
flash_add_greedy = function(data,
                            Kmax = 1,
                            f_init = NULL,
                            var_type = c("by_column", "constant"),
                            init_fn = "udv_si",
                            tol = 1e-2,
                            ebnm_fn = ebnm_pn,
                            ebnm_param = flash_default_ebnm_param(ebnm_fn),
                            verbose = FALSE,
                            nullcheck = TRUE,
                            seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
  var_type = match.arg(var_type)
  f = f_init

  for (k in 1:Kmax) {
    message("fitting factor/loading ", k)
    f = flash_r1(data,
                 f,
                 var_type,
                 init_fn,
                 tol,
                 ebnm_fn,
                 ebnm_param,
                 verbose,
                 nullcheck,
                 seed = NULL)
    
    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, get_k(f)))
      break
  }

  return(f)
}


#' @title Refines a fit of the flash model to data by "backfitting".
#' 
#' @description Iterates through the factors of a flash object,
#'   updating each until convergence.
#' 
#' @param f_init A fitted flash object to be refined.
#' 
#' @param kset The indices of factors to be optimized (NULL indicates
#'   all factors).
#' 
#' @param maxiter A maximum number of iterations to perform (in the
#'   inner loop). To perform just one iteration we suggest setting
#'   \code{maxiter = 1} and \code{nullcheck = FALSE}.
#' 
#' @inheritParams flash_r1
#' 
#' @return A fitted flash object.
#' 
#' @examples
#' 
#' Y = matrix(rnorm(100),nrow=5,ncol=20)
#' fg = flash_add_greedy(Y,10)
#' fb = flash_backfit(Y,fg) # refines fit from greedy by backfitting
#' fsi = flash_init_fn(flash_set_data(Y),"udv_si",4)
#' fb2 = flash_backfit(Y,fsi)
#' flash_get_ldf(fb2)$d
#' 
#' @export
#' 
flash_backfit = function(data,
                         f_init,
                         kset = NULL,
                         var_type = c("by_column", "constant"),
                         tol = 1e-2,
                         ebnm_fn = ebnm_pn,
                         ebnm_param = flash_default_ebnm_param(ebnm_fn),
                         verbose = FALSE,
                         nullcheck = TRUE,
                         maxiter = 1000) {
  f = f_init
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
  if (is.null(kset)) {
    kset = 1:get_k(f)
  }
  var_type = match.arg(var_type)

  if (verbose) {
    message("iteration:1")
  }

  for (k in kset) {
    f = flash_update_single_fl(data, f, k, var_type, ebnm_fn, ebnm_param)
  }

  c = flash_get_objective(data, f)
  if (verbose) {
    message("objective: ", c)
  }

  diff = Inf
  iteration = 2

  while((diff > tol) & (iteration <= maxiter)) {
      
    # There are two steps; first backfit, then null check
    # (if nullcheck removes some factors then the whole process is repeated)
    while((diff > tol) & (iteration <= maxiter)){
      if(verbose){message("iteration:", iteration)}
      for(k in kset){
        f = flash_update_single_fl(data,f,k,var_type,ebnm_fn,ebnm_param)
      }
      cnew = flash_get_objective(data, f)
      diff = cnew - c
      c = cnew
      if(verbose){message("objective: ",c)}
      iteration = iteration + 1
    }

    if (nullcheck) {
      #remove factors that actually hurt objective
      kset = 1:get_k(f)
      f = perform_nullcheck(data, f, kset, var_type, verbose)

      # recompute objective; if it changes then whole process will be repeated
      cnew = flash_get_objective(data, f)
      diff = cnew - c
      c = cnew
      iteration = 1
    }

  }

  return(f)
}

#' @title Fits a rank 1 Empirical Bayes Matrix Factorization model.
#' 
#' @return A fitted flash object. Use \code{flash_get_ldf} to access
#'   standardized loadings and factors; use \code{flash_get_lf} to
#'   access fitted LF'.
#' 
#' @inheritParams flash
#' 
#' @examples
#' 
#' ftrue = rnorm(100)
#' ltrue = rnorm(20)
#' Y = ltrue %*% t(ftrue)+rnorm(2000) # set up a simulated matrix with rank 1 plus noise structure
#' f = flash_r1(Y)
#' ldf = flash_get_ldf(f)
#' plot(ltrue,ldf$l) # plot true l against estimated l (note estimate is normalized);
#' plot(ftrue,ldf$f) # plot true f against estimated f (note estimate is normalized)
#' plot(ltrue %*% t(ftrue), flash_get_lf(f)) #plot true lf' against estimated lf'; the scale of the estimate matches the data
#'
#' # example to use the more flexible ebnm function in ashr; show how to pass parameters to
#' f2 = flash_r1(Y,ebnm_fn = ebnm_ash)
#' # example to show how to pass parameters to ashr
#' f3= flash_r1(Y,ebnm_fn = ebnm_ash, ebnm_param = list(mixcompdist = "normal",method="fdr"))
flash_r1 = function(data,
                    f_init = NULL,
                    var_type = c("by_column", "constant"),
                    init_fn = "udv_si",
                    tol = 1e-2,
                    ebnm_fn = ebnm_pn,
                    ebnm_param = flash_default_ebnm_param(ebnm_fn),
                    verbose = FALSE,
                    nullcheck = TRUE,
                    seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
  var_type = match.arg(var_type)
  f = flash_add_factors_from_data(data,
                                  f_init = f_init,
                                  init_fn = init_fn,
                                  K = 1)
  f = flash_optimize_single_fl(data,
                               f,
                               get_k(f),
                               var_type,
                               nullcheck,
                               tol,
                               ebnm_fn,
                               ebnm_param,
                               verbose)
  return(f)
}
