#' @title Fit Empirical Bayes Matrix Factorization
#'
#' @description This is the main interface for fitting EBMF models
#'   based on algorithms from Wang and Stephens.  The default behaviour
#'   is simply to run the greedy algorithm and return the result.  To
#'   follow it by backfitting set \code{backfit = TRUE}.
#'
#' @param Y An n by p data matrix.
#'
#' @param S An n by p matrix of the standard errors of the observations in
#'   Y. (Can be a scalar if all standard errors are equal.) If
#'   \code{S = NULL}, then the standard errors will be estimated during
#'   fitting.
#'
#' @param greedy_Kmax The maximum number of factors to be greedily added to
#'   the flash object. (This does not include any fixed loadings or
#'   factors.)
#'
#' @param backfit If \code{TRUE}, factors are refined via the backfitting
#'   algorithm. The default behavior is to backfit all factors. Optionally,
#'   a vector may be passed in containing the indices of factors to
#'   backfit.
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
#'       specified via parameter \code{S}.}
#'     \item{\code{"kroneker"}}{This variance type has not yet been
#'       implemented.}
#'   }
#'
#' @param method The type of prior to use when fitting factor/loading
#'   pairs. Options include:
#'   \describe{
#'     \item{\code{"fastest"}}{This method uses point-normal priors, which
#'       can be fit quickly.}
#'     \item{\code{"normalmix"}}{A more flexible class of normal-mixture
#'       ash priors.}
#'     \item{\code{"nonnegative"}}{Nonnegative uniform-mixture ash priors
#'       on factors and loadings.}
#'     \item{\code{"nnfactors"}}{Nonnegative priors on factors and
#'       normal-mixture priors on loadings.}
#'     \item{\code{"nnloadings"}}{Nonnegative priors on loadings and
#'       normal-mixture priors on factors.}
#'   }
#'
#' @param f_init The flash object or flash fit object to which new
#'   factors are to be added. If \code{f_init = NULL}, then a new flash
#'   object is created.
#'
#' @param fixed_loadings An n x K matrix of loadings to be added. The
#'   default behavior is to fix supplied values and to allow missing
#'   values to be re-estimated during fitting. If a different behavior is
#'   desired, a list that includes fields \code{vals} and \code{is_fixed}
#'   can be passed in. Both fields must be n x K matrices.
#'
#' @param fixed_factors A p x K matrix of factors or a list that includes
#'   fields \code{vals} and \code{is_fixed} (both of which must be p x K
#'   matrices).
#'
#' @param verbose If \code{TRUE}, various progress updates will be
#'   printed.
#'
#' @param control A list of additional parameters used during optimization.
#'  The following parameters may be specified:
#'  \describe{
#'    \item{\code{init_fn}}{The function used to initialize factors.
#'      Options include:
#'      \describe{
#'        \item{\code{"udv_si"}}{Provides a simple wrapper to
#'          \code{\link[softImpute]{softImpute}} to provide a rank-one
#'          initialization. Uses option \code{type = "als"}.}
#'        \item{\code{"udv_si_svd"}}{Uses
#'          \code{\link[softImpute]{softImpute}} with option
#'          \code{type = "svd"}.}
#'        \item{\code{"udv_svd"}}{Provides a simple wrapper to \code{svd}.}
#'        \item{\code{"udv_random"}}{Provides a random initialization of
#'          factors.}
#'      }
#'      A user-specified function can also be passed in. This function should
#'      take parameters \code{(Y, K)}, where \code{Y} is an n by p matrix of
#'      data (or a flash data object) and \code{K} is the number of factors.
#'      It should output a list with elements \code{(u, d, v)}, where
#'      \code{u} is a n by K matrix, \code{v} is a p by K matrix, and
#'      \code{d} is a K vector. (If the input data includes missing values,
#'      then the function must be able to deal with missing values in its
#'      input matrix.)
#'    }
#'    \item{\code{ebnm_fn}}{The function used to solve the Empirical Bayes
#'      Normal Means problem. Either a single character string (giving the
#'      name of the function) or a list with fields \code{l} and \code{f}
#'      (specifying different functions to be used for loadings and factors)
#'      are acceptable arguments. Options include:
#'      \describe{
#'        \item{\code{"ebnm_ash"}}{A wrapper to the function
#'          \code{\link[ashr]{ash}}.}
#'        \item{\code{"ebnm_pn"}}{A wrapper to function
#'          \code{\link[ebnm]{ebnm_point_normal}} in package \pkg{ebnm}.}
#'        \item{\code{"ebnm_pl"}}{A wrapper to function
#'          \code{\link[ebnm]{ebnm_point_laplace}} in \pkg{ebnm}.}
#'      }
#'    }
#'    \item{\code{ebnm_param}}{A named list containing parameters to be
#'      passed to \code{ebnm_fn} when optimizing. A list with fields
#'      \code{l} and \code{f} (each of which is a named list) will
#'      separately supply parameters for loadings and factors. If parameter
#'      \code{warmstart} is used, the current value of \code{g} (if
#'      available) will be passed to \code{ebnm_fn}. (So, \code{ebnm_fn}
#'      should accept a parameter named \code{g}, not one named
#'      \code{warmstart}.) Set \code{ebnm_param} to \code{NULL} to use
#'      defaults.}
#'    \item{\code{stopping_rule}}{TODO.}
#'    \item{\code{tol}}{Specifies how much the objective can change in a
#'      single iteration to be considered not converged.}
#'    \item{\code{r1opt_maxiter}}{TODO.}
#'    \item{\code{backfit_maxiter}}{TODO.}
#'    \item{\code{nullcheck}}{If \code{TRUE}, then after running
#'      hill-climbing updates \code{flash} will check whether the achieved
#'      optimum is better than setting the factor to zero. If the check is
#'      performed and fails then the factor will be set to zero in the
#'      returned fit.}
#'    \item{\code{verbose_output}}{TODO.}
#'    \item{\code{seed}}{A random number seed to use before running
#'      \code{flash} - for reproducibility. Set to \code{NULL} if you
#'      don't want the seed set. (The seed can affect initialization when
#'      there is missing data; otherwise the algorithm is deterministic.)}
#'  }
#'
#' @return A flash object.
#'
#' @examples
#'
#' set.seed(1) # for reproducibility
#' ftrue = matrix(rnorm(200), ncol=2)
#' ltrue = matrix(rnorm(40), ncol=2)
#' ltrue[1:10, 1] = 0 # set up some sparsity
#' ltrue[11:20, 2] = 0
#' Y = ltrue %*% t(ftrue) + rnorm(2000) # set up a simulated matrix
#' f = flash(Y, greedy_Kmax=1)
#' ldf = f$ldf
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
#' plot(ltrue %*% t(ftrue), f$fitted_values)
#'
#' # Example to use the more flexible ebnm function in ashr.
#' f2 = flash(Y, method="normalmix")
#'
#' # Example to show how to pass parameters to ashr (may be most
#' # useful for research use).
#' f3 = flash(Y, control=list(ebnm_fn="ebnm_ash",
#'                            ebnm_param=list(mixcompdist="normal",
#'                                            method="fdr")))
#'
#' # Example to show how to separately specify parameters for factors
#' # and loadings.
#' f4 = flash(Y, control=list(ebnm_fn=list(l="ebnm_pn", f="ebnm_ash"),
#'                            ebnm_param=list(l=list(),
#'                                            f=list(g=ashr::normalmix(1,0,1),
#'                                                   fixg=TRUE))))
#'
#' # Example to show how to use a different initialization function.
#' library(softImpute)
#' init_fn = function(x, K=1) {softImpute(x, K, lambda=10)}
#' f5 = flash(Y, control=list(init_fn=init_fn))
#'
#' @export
#'
flash = function(Y,
                 S = NULL,
                 greedy_Kmax = 0,
                 backfit = FALSE,
                 var_type = c("by_column",
                              "by_row",
                              "constant",
                              "zero",
                              "kroneker"),
                 method = c("fastest",
                            "normalmix",
                            "nonnegative",
                            "nnfactors",
                            "nnloadings"),
                 f_init = NULL,
                 fixed_loadings = NULL,
                 fixed_factors = NULL,
                 verbose = TRUE,
                 control = list()) {
  if (is.null(fixed_loadings)
      && is.null(fixed_factors)
      && greedy_Kmax < 1
      && (identical(backfit, FALSE) || is.null(f_init))) {
    stop(paste("Nothing to do. Set greedy_Kmax > 0 to add factors from",
               "data. Set fixed_loadings or fixed_factors to add",
               "specified factors. Set backfit to backfit a new or",
               "existing flash object."))
  }

  data = handle_data(flash_set_data(Y, S), f_init)
  var_type = handle_var_type(match.arg(var_type), data)
  fl = handle_f(f_init, init_null_f = TRUE)
  LL_init = handle_fixed(fixed_loadings, flash_get_n(fl))
  FF_init = handle_fixed(fixed_factors, flash_get_p(fl))
  Kmax = LL_init$K + FF_init$K + greedy_Kmax
  backfit_kset = handle_backfit(backfit)

  params = get_control_defaults(match.arg(method))
  if (any(!is.element(names(control), names(params)))) {
    stop("Control argument contains unknown parameter names.")
  }
  params = modifyList(params, control, keep.null = TRUE)

  init_fn = handle_init_fn(params$init_fn)
  ebnm_fn = handle_ebnm_fn(params$ebnm_fn)

  # Do not allow different ebnm_params for each factor/loading if more than
  #   one type of fit (among add_fixed_loadings, add_fixed_factors,
  #   add_greedy, and backfit) is being performed (otherwise the indexing
  #   gets too confusing):
  allow_lists_of_lists = ((LL_init$K > 0) +
                            (FF_init$K > 0) +
                            (greedy_Kmax > 0) +
                            (!identical(backfit, FALSE)) == 1)

  if (allow_lists_of_lists) {
    if (identical(backfit, TRUE)) {
      n_expected = flash_get_k(fl) # backfitting the entire flash object
    } else {
      n_expected = max(Kmax, length(backfit_kset))
    }
  } else {
    n_expected = flash_get_k(fl) + Kmax
  }
  ebnm_param = handle_ebnm_param(params$ebnm_param, ebnm_fn, n_expected,
                                 allow_lists_of_lists)

  # TODO: handle stopping rule, verbose_output, tol, Kmax, maxiter
  stopping_rule = params$stopping_rule
  tol = params$tol
  r1opt_maxiter = params$r1opt_maxiter
  backfit_maxiter = params$backfit_maxiter
  nullcheck = params$nullcheck
  seed = params$seed

  if (!verbose) {
    params$verbose_output = ""
  }
  verbose_output = unlist(strsplit(params$verbose_output, split = NULL))

  history = list()

  if (LL_init$K > 0) {
    verbose_fixed_loadings_announce(LL_init$K)

    res = add_fixed_loadings(data = data,
                             f_init = fl,
                             LL = LL_init$vals,
                             fixl = LL_init$is_fixed,
                             var_type = var_type,
                             ebnm_fn = ebnm_fn,
                             ebnm_param = ebnm_param,
                             stopping_rule = stopping_rule,
                             tol = tol,
                             verbose_output = verbose_output,
                             nullcheck = nullcheck,
                             maxiter = r1opt_maxiter)
    fl = res$f
    history = c(history, res$history)
  }

  if (FF_init$K > 0) {
    verbose_fixed_factors_announce(FF_init$K)

    res = add_fixed_factors(data = data,
                            f_init = fl,
                            FF = FF_init$vals,
                            fixf = FF_init$is_fixed,
                            var_type = var_type,
                            ebnm_fn = ebnm_fn,
                            ebnm_param = ebnm_param,
                            stopping_rule = stopping_rule,
                            tol = tol,
                            verbose_output = verbose_output,
                            nullcheck = nullcheck,
                            maxiter = r1opt_maxiter)
    fl = res$f
    history = c(history, res$history)
  }

  if (greedy_Kmax > 0) {
    res = add_greedy(data = data,
                     Kmax = greedy_Kmax,
                     f_init = fl,
                     var_type = var_type,
                     init_fn = init_fn,
                     ebnm_fn = ebnm_fn,
                     ebnm_param = ebnm_param,
                     stopping_rule = stopping_rule,
                     tol = tol,
                     verbose_output = verbose_output,
                     nullcheck = nullcheck,
                     maxiter = r1opt_maxiter,
                     seed = seed)
    fl = res$f
    history = c(history, res$history)
  }

  if (!identical(backfit, FALSE) && backfit_maxiter > 0) {
    backfit_kset = handle_kset(backfit_kset, fl)
    res = backfit(data = data,
                  f = fl,
                  kset = backfit_kset,
                  var_type = var_type,
                  ebnm_fn = ebnm_fn,
                  ebnm_param = ebnm_param,
                  stopping_rule = stopping_rule,
                  tol = tol,
                  verbose_output = verbose_output,
                  nullcheck = nullcheck,
                  maxiter = backfit_maxiter)
    fl = res$f
    history = c(history, res$history)
  }

  # If a factor is added without doing any optimization, then the
  #   objective will be not valid.
  if (Kmax > 0 && r1opt_maxiter < 1 &&
      (identical(backfit, FALSE) || backfit_maxiter < 1)) {
    compute_obj = FALSE
  } else {
    compute_obj = TRUE
  }

  flash_object = construct_flash_object(data = data,
                                        fit = fl,
                                        history = history,
                                        f_init = f_init,
                                        compute_obj = compute_obj)

  return(flash_object)
}
