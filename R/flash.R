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
#' @param Kmax The maximum number of factors to be added to the flash
#'   object. (If \code{nullcheck = TRUE}, the actual number of factors
#'   added might be less than \code{Kmax}.)
#'
#' @param f_init The flash object to which new factors are to be added.
#'   If \code{f_init = NULL}, then a new flash fit object is created.
#'
#' @param var_type The type of variance structure to assume for
#'   residuals.
#'
#' @param init_fn The function used to initialize factors. This
#'   function should take parameters (Y,K) where Y is an n by p matrix
#'   of data (or a flash data object) and K is a number of factors.  It
#'   should output a list with elements (u,d,v) where u is n by K matrix
#'   v is a p by K matrix and d is a K vector. See \code{udv_si} for an
#'   example.  (If the input data includes missing values then this
#'   function must be able to deal with missing values in its input
#'   matrix.)
#'
#' @param tol Specifies how much the objective can change in a single
#'   iteration to be considered not converged.
#'
#' @param ebnm_fn The function used to solve the Empirical Bayes Normal
#'   Means problem.
#'
#' @param ebnm_param A named list containing parameters to be passed to
#'   ebnm_fn when optimizing; defaults are set by
#'   \code{flash_default_ebnm_param()}.
#'
#' @param verbose If TRUE, various progress updates will be printed.
#'
#' @param nullcheck If TRUE, then after running hill-climbing updates,
#'   \code{flash} will check whether the achieved optimum is better than
#'   setting the factor to 0. If the check is performed and fails then
#'   the factor will be set to 0 in the returned fit.
#'
#' @param seed A random number seed to use before running \code{flash}
#'   - for reproducibility. Set to NULL if you don't want the seed set.
#'   (The seed can affect initialization when there are missing data;
#'   otherwise the algorithm is deterministic.)
#'
#' @param greedy If TRUE, factors are added via the greedy algorithm.
#'   If FALSE, then \code{f_init} must be supplied.
#'
#' @param backfit If TRUE, factors are refined via the backfitting
#'   algorithm.
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
#'
#' # Show the weights, analogous to singular values showing importance
#' # of each factor.
#' ldf$d
#'
#' # Plot true l against estimated l; with this seed it turns out the
#' # 2nd loading/factor corresponds to the first column of ltrue.
#' plot(ltrue[,1],ldf$l[,2])
#'
#' # Plot true f against estimated f (note estimate is normalized).
#' plot(ftrue[,1],ldf$f[,2])
#'
#' # Plot true lf' against estimated lf'; the scale of the estimate
#' # matches the data.
#' plot(ltrue %*% t(ftrue), flash_get_lf(f))
#'
#' # Example to use the more flexible ebnm function in ashr.
#' f2 = flash(Y,ebnm_fn = ebnm_ash)
#'
#' # Example to show how to pass parameters to ashr (may be most
#' # useful for research use).
#' f3 = flash(Y,ebnm_fn = ebnm_ash,
#'            ebnm_param = list(mixcompdist = "normal",method="fdr"))
#'
#' # Example to show how to use a different initialization function.
#' library(softImpute)
#' f4 = flash(Y,init_fn = function(x,K=1){softImpute(x,K,lambda=10)})
#'
#' @export
#'
flash = function(data,
                 Kmax = 100,
                 f_init = NULL,
                 var_type = c("by_column", "by_row", "constant",
                              "zero", "kroneker"),
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
    f = flash_add_greedy(data,
                         Kmax,
                         f_init,
                         var_type,
                         init_fn,
                         tol,
                         ebnm_fn,
                         ebnm_param,
                         gl = NULL,
                         fixgl = FALSE,
                         gf = NULL,
                         fixgf = FALSE,
                         verbose,
                         nullcheck,
                         seed
    )
  } else {
    f = f_init
  }

  if (backfit) {
    f = flash_backfit(data,
                      f,
                      kset = NULL,
                      var_type,
                      tol,
                      ebnm_fn,
                      ebnm_param,
                      gl = NULL,
                      fixgl = FALSE,
                      gf = NULL,
                      fixgf = FALSE,
                      verbose,
                      nullcheck)
  }
  return(f)
}

#' @title Fit Empirical Bayes Matrix Factorization (greedy algorithm)
#'
#' @description This implements the greedy algorithm from Wang and
#'   Stephens.  It can be used to add factors to an existing fit, or
#'   start from scratch.  It adds factors iteratively, at each stage
#'   adding a new factor and then optimizing it.  It is "greedy" in that
#'   it does not return to re-optimize previous factors.  The function
#'   stops when an added factor contributes nothing or when \code{Kmax}
#'   is reached. Each new factor is intialized by applying the function
#'   \code{init_fn} to the residuals after removing previously-fitted
#'   factors.
#'
#' @inheritParams flash
#'
#' @param gl Passed into \code{ebnm_fn} as parameter \code{g} (used to
#'   fix or initialize priors on the loadings). This can be a single
#'   prior or a list of length \code{Kmax}, with \code{gl[[k]]}
#'   specifying the prior for the kth loading.
#'
#' @param fixgl Passed into \code{ebnm_fn} as parameter \code{fixg} (used
#'   to fix priors on the loadings). This can be a single boolean which
#'   specifies \code{fixg} for all loadings or a vector of booleans,
#'   with \code{fixg}[k] specifying \code{fixg} for the kth loading.
#'
#' @param gl Passed into \code{ebnm_fn} as parameter \code{g} (used to
#'   fix or initialize priors on the factors). This can be a single
#'   prior or a list of length \code{Kmax}, with \code{gf[[k]]}
#'   specifying the prior for the kth factor.
#'
#' @param fixgl Passed into \code{ebnm_fn} as parameter \code{fixg} (used
#'   to fix priors on the factors). This can be a single boolean which
#'   specifies \code{fixg} for all factors or a vector of booleans,
#'   with \code{fixg}[k] specifying \code{fixg} for the kth factor.
#'
#' @return A fitted flash object.
#'
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l,f) + matrix(rnorm(1000),nrow=100)
#' f = flash_add_greedy(Y,10)
#'
#' # Gives the weights for each factor (analogue of singular values).
#' flash_get_ldf(f)$d
#'
#' # Example to show how to use a different initialization function.
#' library(softImpute)
#' f2 = flash_add_greedy(Y,10,init_fn = function(x,K=1){
#'   softImpute(x,K,lambda=10)
#' })
#'
#' @export
#'
flash_add_greedy = function(data,
                            Kmax = 1,
                            f_init = NULL,
                            var_type = c("by_column", "by_row", "constant",
                                         "zero", "kroneker"),
                            init_fn = "udv_si",
                            tol = 1e-2,
                            ebnm_fn = ebnm_pn,
                            ebnm_param = flash_default_ebnm_param(ebnm_fn),
                            gl = NULL,
                            fixgl = FALSE,
                            gf = NULL,
                            fixgf = FALSE,
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
  if (!is.list(gl[[1]])) {
    gl = rep(list(gl), Kmax)
  }
  if (length(fixgl) == 1) {
    fixgl = rep(fixgl, Kmax)
  }
  if (!is.list(gf[[1]])) {
    gf = rep(list(gf), Kmax)
  }
  if (length(fixgf) == 1) {
    fixgf = rep(fixgf, Kmax)
  }
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
                 gl[[k]],
                 fixgl[k],
                 gf[[k]],
                 fixgf[k],
                 verbose,
                 nullcheck,
                 seed = NULL)

    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, flash_get_k(f)))
      break
  }

  return(f)
}


#' @title Refines a fit of the flash model to data by "backfitting".
#'
#' @description Iterates through the factors of a flash object,
#'   updating each until convergence.
#'
#' @param kset The indices of factors to be optimized (NULL indicates
#'   all factors).
#'
#' @param maxiter A maximum number of iterations to perform (in the
#'   inner loop). To perform just one iteration we suggest setting
#'   \code{maxiter = 1} and \code{nullcheck = FALSE}.
#'
#' @inheritParams flash
#'
#' @param gl If nonnull, then the priors on the loadings will be fixed at
#'   the specified values. This should be a list, with \code{gl[[k]]}
#'   specifying the prior for the kth loading.
#'
#' @param gf If nonnull, then the priors on the factors will be fixed
#'   the specified values. This should be a list, with \code{gf[[k]]}
#'   specifying the prior for the kth factor.
#'
#' @return A fitted flash object.
#'
#' @examples
#'
#' l = rnorm(100) #simulate some rank 1 data
#' f = rnorm(10)
#' Y = outer(l,f) + matrix(rnorm(1000),nrow=100)
#' fg = flash_add_greedy(Y,10)
#' fb = flash_backfit(Y,fg) # refines fit from greedy by backfitting
#' flash_get_ldf(fb)$d
#'
#' @export
#'
flash_backfit = function(data,
                         f_init,
                         kset = NULL,
                         var_type = c("by_column", "by_row", "constant",
                                      "zero", "kroneker"),
                         tol = 1e-2,
                         ebnm_fn = ebnm_pn,
                         ebnm_param = flash_default_ebnm_param(ebnm_fn),
                         gl = NULL,
                         fixgl = FALSE,
                         gf = NULL,
                         fixgf = FALSE,
                         verbose = FALSE,
                         nullcheck = TRUE,
                         maxiter = 1000) {
  f = f_init
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
  if (is.null(kset)) {
    kset = 1:flash_get_k(f)
  }
  var_type = match.arg(var_type)
  if (!is.list(gl[[1]])) {
    gl = rep(list(gl), max(kset))
  }
  if (length(fixgl) == 1) {
    fixgl = rep(fixgl, max(kset))
  }
  if (!is.list(gf[[1]])) {
    gf = rep(list(gf), max(kset))
  }
  if (length(fixgf) == 1) {
    fixgf = rep(fixgf, max(kset))
  }

  if (verbose) {
    message("iteration:1")
  }

  for (k in kset) {
    f = flash_update_single_fl(data,
                               f,
                               k,
                               var_type,
                               ebnm_fn,
                               ebnm_param,
                               gl[[k]],
                               fixgl[k],
                               gf[[k]],
                               fixgf[k])
  }

  c = flash_get_objective(data, f)
  if (verbose) {
    message("objective: ", c)
  }

  diff = Inf
  iteration = 2

  while((diff > tol) & (iteration <= maxiter)) {

    # There are two steps; first backfit, then null check (if nullcheck
    # removes some factors then the whole process is repeated).
    while((diff > tol) & (iteration <= maxiter)){
      if(verbose){message("iteration:", iteration)}
      for(k in kset){
        f = flash_update_single_fl(data,
                                   f,
                                   k,
                                   var_type,
                                   ebnm_fn,
                                   ebnm_param,
                                   gl[[k]],
                                   fixgl[k],
                                   gf[[k]],
                                   fixgf[k])
      }
      cnew = flash_get_objective(data, f)
      diff = cnew - c
      c = cnew
      if(verbose){message("objective: ",c)}
      iteration = iteration + 1
    }

    if (nullcheck) {
      # remove factors that actually hurt objective
      kset = 1:flash_get_k(f)
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

# @title Fits a rank 1 Empirical Bayes Matrix Factorization model.
#
# @return A fitted flash object. Use \code{flash_get_ldf} to access
#   standardized loadings and factors; use \code{flash_get_lf} to
#   access fitted LF'.
#
# @inheritParams flash_add_greedy
#
# @examples
#
# ftrue = rnorm(100)
# ltrue = rnorm(20)
#
# # Set up a simulated matrix with rank 1 plus noise structure.
# Y = ltrue %*% t(ftrue)+rnorm(2000)
# f = flash_r1(Y)
# ldf = flash_get_ldf(f)
#
# # Plot true l against estimated l (note estimate is normalized).
# plot(ltrue,ldf$l)
#
# # Plot true f against estimated f (note estimate is normalized).
# plot(ftrue,ldf$f)
#
# # Plot true lf' against estimated lf'; the scale of the estimate
# # matches the data.
# plot(ltrue %*% t(ftrue), flash_get_lf(f))
#
# # Example to use the more flexible ebnm function in ashr.
# f2 = flash_r1(Y,ebnm_fn = ebnm_ash)
#
# # Example to show how to pass parameters to ashr.
# f3 = flash_r1(Y,ebnm_fn = ebnm_ash,
#               ebnm_param = list(mixcompdist = "normal",method="fdr"))
#
flash_r1 = function(data,
                    f_init = NULL,
                    var_type = c("by_column", "by_row", "constant",
                                 "zero", "kroneker"),
                    init_fn = "udv_si",
                    tol = 1e-2,
                    ebnm_fn = ebnm_pn,
                    ebnm_param = flash_default_ebnm_param(ebnm_fn),
                    gl = NULL,
                    fixgl = FALSE,
                    gf = NULL,
                    fixgf = FALSE,
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
                                  K = 1,
                                  f_init = f_init,
                                  init_fn = init_fn)
  f = flash_optimize_single_fl(data,
                               f,
                               flash_get_k(f),
                               var_type,
                               nullcheck,
                               tol,
                               ebnm_fn,
                               ebnm_param,
                               gl,
                               fixgl,
                               gf,
                               fixgf,
                               verbose)
  return(f)
}
