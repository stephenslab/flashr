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
#'   object.
#'
#' @param f_init The flash object to which new factors are to be added.
#'   If \code{f_init = NULL}, then a new flash fit object is created.
#'
#' @param var_type The type of variance structure to assume for
#'   residuals. Options include:
#'   \describe{
#'     \item{\code{"by_column"}}{Residuals in any given column are
#'       assumed to have the same variance.}
#'     \item{\code{"by_row"}}{Residuals in any given row have the
#'       same variance.}
#'     \item{\code{"constant"}}{All residuals are assumed to have the
#'       same variance.}
#'     \item{\code{"zero"}}{The variance of the residuals is fixed. To
#'       use this variance type, the standard errors must be
#'       specified via parameter \code{S} when using
#'       \code{flash_set_data} to set the flash data object.}
#'     \item{\code{"kroneker"}}{This variance type has not yet been
#'       implemented.}
#'   }
#'
#' @param init_fn The function used to initialize factors. Options
#'   include:
#'   \describe{
#'     \item{\code{"udv_si"}}{Provides a simple wrapper to
#'       \code{\link[softImpute]{softImpute}} to provide a rank-one
#'       initialization. Uses option \code{type = "als"}.}
#'     \item{\code{"udv_si_svd"}}{Uses
#'       \code{\link[softImpute]{softImpute}} with option
#'       \code{type = "svd"}.}
#'     \item{\code{"udv_svd"}}{Provides a simple wrapper to \code{svd}.}
#'     \item{\code{"udv_random"}}{Provides a random initialization of
#'       factors.}
#'   }
#'   A user-specified function can also be used. This function should
#'   take parameters \code{(Y, K)}, where \code{Y} is an n by p matrix of
#'   data (or a flash data object) and \code{K} is the number of factors.
#'   It should output a list with elements \code{(u, d, v)}, where
#'   \code{u} is a n by K matrix, \code{v} is a p by K matrix, and
#'   \code{d} is a K vector. (If the input data includes missing values,
#'   then the function must be able to deal with missing values in its
#'   input matrix.)
#'
#' @param tol Specifies how much the objective can change in a single
#'   iteration to be considered not converged.
#'
#' @param ebnm_fn The function used to solve the Empirical Bayes Normal
#'   Means problem. Either a single character string (giving the name of
#'   of the function) or a list with fields \code{l} and \code{f}
#'   (specifying different functions to be used for loadings and factors)
#'   are acceptable arguments. Options include:
#'   \describe{
#'     \item{\code{"ebnm_ash"}}{A wrapper to the function
#'       \code{\link[ashr]{ash}}.}
#'     \item{\code{"ebnm_pn"}}{A wrapper to function
#'       \code{\link[ebnm]{ebnm_point_normal}} in package \pkg{ebnm}.}
#'     \item{\code{"ebnm_pl"}}{A wrapper to function
#'       \code{\link[ebnm]{ebnm_point_laplace}} in \pkg{ebnm}.}
#'   }
#'
#' @param ebnm_param A named list containing parameters to be passed to
#'   \code{ebnm_fn} when optimizing. A list with fields \code{l} and
#'   \code{f} (each of which is a named list) will separately supply
#'   parameters for loadings and factors. Set to \code{NULL} to use
#'   defaults.
#'
#' @param verbose If \code{TRUE}, various progress updates will be
#'   printed.
#'
#' @param nullcheck If \code{TRUE}, then after running hill-climbing
#'   updates \code{flash} will check whether the achieved optimum is
#'   better than setting the factor to zero. If the check is performed
#'   and fails then the factor will be set to zero in the returned fit.
#'
#' @param seed A random number seed to use before running \code{flash}
#'   - for reproducibility. Set to \code{NULL} if you don't want the
#'   seed set. (The seed can affect initialization when there are
#'   missing data; otherwise the algorithm is deterministic.)
#'
#' @param greedy If \code{TRUE}, factors are added via the greedy
#'   algorithm. If \code{FALSE}, then \code{f_init} must be supplied.
#'
#' @param backfit If \code{TRUE}, factors are refined via the backfitting
#'   algorithm.
#'
#' @seealso \code{\link{flash_add_greedy}}, \code{\link{flash_backfit}}
#'
#' @return A fitted flash object. Use \code{flash_get_ldf} to access
#'   standardized loadings and factors; use
#'   \code{flash_get_fitted_values} to access fitted LF'.
#'
#' @examples
#'
#' set.seed(1) # for reproducibility
#' ftrue = matrix(rnorm(200), ncol=2)
#' ltrue = matrix(rnorm(40), ncol=2)
#' ltrue[1:10, 1] = 0 # set up some sparsity
#' ltrue[11:20, 2] = 0
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
#' plot(ltrue[,1], ldf$l[,2])
#'
#' # Plot true f against estimated f (note estimate is normalized).
#' plot(ftrue[,1], ldf$f[,2])
#'
#' # Plot true lf' against estimated lf'; the scale of the estimate
#' # matches the data.
#' plot(ltrue %*% t(ftrue), flash_get_fitted_values(f))
#'
#' # Example to use the more flexible ebnm function in ashr.
#' f2 = flash(Y, ebnm_fn="ebnm_ash")
#'
#' # Example to show how to pass parameters to ashr (may be most
#' # useful for research use).
#' f3 = flash(Y,
#'            ebnm_fn="ebnm_ash",
#'            ebnm_param=list(mixcompdist="normal", method="fdr"))
#'
#' # Example to show how to separately specify parameters for factors
#' # and loadings.
#' f4 = flash(Y,
#'            ebnm_fn=list(l="ebnm_pn", f="ebnm_ash"),
#'            ebnm_param=list(l=list(),
#'                            f=list(g=ashr::normalmix(1,0,1), fixg=TRUE)))
#'
#' # Example to show how to use a different initialization function.
#' library(softImpute)
#' f5 = flash(Y, init_fn=function(x, K=1){softImpute(x, K, lambda=10)})
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
                 ebnm_fn = "ebnm_pn",
                 ebnm_param = NULL,
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
                         verbose,
                         nullcheck,
                         seed)
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
                      verbose,
                      nullcheck)
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
#'   Each new factor is intialized by applying the function
#'   \code{init_fn} to the residuals after removing previously-fitted
#'   factors.
#'
#' @inheritParams flash
#'
#' @return A fitted flash object.
#'
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l, f) + matrix(rnorm(1000), nrow=100)
#' f = flash_add_greedy(Y,10)
#'
#' # Gives the weights for each factor (analogue of singular values).
#' flash_get_ldf(f)$d
#'
#' # Example to show how to use a different initialization function.
#' library(softImpute)
#' f2 = flash_add_greedy(Y, 10, init_fn=function(x, K=1) {
#'   softImpute(x, K, lambda=10)
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
                            ebnm_fn = "ebnm_pn",
                            ebnm_param = NULL,
                            verbose = FALSE,
                            nullcheck = TRUE,
                            seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data = handle_data(data)
  var_type = match.arg(var_type)
  init_fn = handle_init_fn(init_fn)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, Kmax)

  f = f_init

  for (k in 1:Kmax) {
    message("fitting factor/loading ", k)
    old_f = f
    f = flash_r1(data,
                 f,
                 var_type,
                 init_fn,
                 tol,
                 ebnm_fn$l,
                 ebnm_param$l[[k]],
                 ebnm_fn$f,
                 ebnm_param$f[[k]],
                 verbose,
                 nullcheck,
                 maxiter = 5000,
                 seed = NULL)

    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, flash_get_k(f))) {
      if (flash_get_k(f) > 1) {
        # Remove zero factor as long as it doesn't create an empty object.
        f = old_f
      }
      break
    }
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
#' @param kset The indices of factors to be optimized (\code{NULL}
#'   indicates all factors).
#'
#' @param maxiter A maximum number of iterations to perform (not
#'   including repeated fittings if \code{nullcheck} fails). To perform
#'   just one iteration we suggest setting \code{maxiter = 1} and
#'   \code{nullcheck = FALSE}.
#'
#' @inheritParams flash
#'
#' @param ebnm_param A named list containing parameters to be passed to
#'   \code{ebnm_fn} when optimizing. A list with fields \code{l} and
#'   \code{f} (each of which is a named list) will separately supply
#'   parameters for the loadings and factors. An unnamed list of
#'   \code{length(kset)} named lists will separately supply parameters
#'   for each factor/loading in \code{kset}. Finally, a list with fields
#'   \code{l} and \code{f}, each of which contains an unnamed list of
#'   \code{length(kset)} named lists, will separately supply parameters
#'   for each distinct loading and each distinct factor. Set to
#'   \code{NULL} to use defaults.
#'
#' @return A fitted flash object.
#'
#' @examples
#'
#' LL = matrix(rnorm(200), ncol=2) # simulate some rank 2 data
#' FF = matrix(rnorm(20), nrow=2)
#' Y = LL %*% FF + matrix(rnorm(1000), nrow=100)
#' fg = flash_add_greedy(Y, 10)
#' fb = flash_backfit(Y, fg) # refines fit from greedy by backfitting
#' flash_get_ldf(fb)$d
#'
#' # Example to illustrate different types of arguments to ebnm_param.
#' # 1. Fix a N(0, 1) prior on the loadings.
#' ebnm_param_l = list(g=ashr::normalmix(1,0,1), fixg=TRUE)
#' fg2 = flash_add_greedy(Y, 10, ebnm_fn="ebnm_ash",
#'                        ebnm_param=list(l=ebnm_param_l, f=list()))
#' # 2. Now refit factors, forcing loadings to use prior from greedy fit.
#' ebnm_param_f = lapply(fg2$gf, function(g) {list(g=g, fixg=TRUE)})
#' fb2 = flash_backfit(Y, fg2, kset=1:2, ebnm_fn="ebnm_ash",
#'                     ebnm_param=list(l=list(), f=ebnm_param_f))
#'
#' @export
#'
flash_backfit = function(data,
                         f_init,
                         kset = NULL,
                         var_type = c("by_column", "by_row", "constant",
                                      "zero", "kroneker"),
                         tol = 1e-2,
                         ebnm_fn = "ebnm_pn",
                         ebnm_param = NULL,
                         verbose = FALSE,
                         nullcheck = TRUE,
                         maxiter = 1000) {
  data = handle_data(data)
  kset = handle_kset(kset, f_init)
  var_type = match.arg(var_type)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, length(kset))

  f = f_init

  if (verbose) {
    message("iteration:1")
  }

  for (i in 1:length(kset)) {
    f = flash_update_single_fl(data,
                               f,
                               kset[i],
                               var_type,
                               ebnm_fn$l,
                               ebnm_param$l[[i]],
                               ebnm_fn$f,
                               ebnm_param$f[[i]])
  }

  c = flash_get_objective(data, f)
  if (verbose) {
    message("objective: ", c)
  }

  diff = Inf
  iteration = 2

  while((diff > tol) & (iteration <= maxiter)) {

    # There are two steps; first backfit, then null check
    # (if nullcheck removes some factors then the whole process
    # is repeated)
    while((diff > tol) & (iteration <= maxiter)){
      if(verbose){message("iteration:", iteration)}
      for (i in 1:length(kset)) {
        f = flash_update_single_fl(data,
                                   f,
                                   kset[i],
                                   var_type,
                                   ebnm_fn$l,
                                   ebnm_param$l[[i]],
                                   ebnm_fn$f,
                                   ebnm_param$f[[i]])
      }
      cnew = flash_get_objective(data, f)
      diff = cnew - c
      c = cnew
      if(verbose){message("objective: ",c)}
      iteration = iteration + 1
    }

    if (nullcheck) {
      # Remove factors that actually hurt the objective.
      kset = 1:flash_get_k(f)
      f = perform_nullcheck(data, f, kset, var_type, verbose)

      # Recompute objective; if it changes then the whole process will
      # be repeated.
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
# @return A fitted flash object.
#
# @inheritParams flash
#
flash_r1 = function(data,
                    f_init,
                    var_type,
                    init_fn,
                    tol,
                    ebnm_fn_l,
                    ebnm_param_l,
                    ebnm_fn_f,
                    ebnm_param_f,
                    verbose,
                    maxiter,
                    nullcheck,
                    seed) {
  f = flash_add_factors_from_data(data,
                                  K = 1,
                                  f_init,
                                  init_fn)
  f = flash_optimize_single_fl(data,
                               f,
                               flash_get_k(f),
                               var_type,
                               nullcheck,
                               tol,
                               ebnm_fn_l,
                               ebnm_param_l,
                               ebnm_fn_f,
                               ebnm_param_f,
                               verbose,
                               maxiter)
  return(f)
}
