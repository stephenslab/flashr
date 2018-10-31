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
#' @return A flash object.
#'
#' @examples
#' l = rnorm(100)
#' f = rnorm(10)
#' Y = outer(l, f) + matrix(rnorm(1000), nrow=100)
#' f = flash_add_greedy(Y,10)
#'
#' # Gives the weights for each factor (analogue of singular values).
#' f$ldf$d
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
                            Kmax = 100,
                            f_init = NULL,
                            var_type = c("by_column",
                                         "by_row",
                                         "constant",
                                         "zero",
                                         "kroneker"),
                            init_fn = "udv_si",
                            tol = 1e-2,
                            ebnm_fn = "ebnm_pn",
                            ebnm_param = NULL,
                            verbose = TRUE,
                            nullcheck = TRUE,
                            seed = 123) {

  if (verbose) {
    verbose_output = "odn" # objective, obj diff, nullcheck
  } else {
    verbose_output = ""
  }

  flash_object = flash_greedy_workhorse(data,
                                        Kmax,
                                        f_init,
                                        var_type,
                                        init_fn,
                                        tol,
                                        ebnm_fn,
                                        ebnm_param,
                                        verbose_output,
                                        nullcheck,
                                        seed)

  return(flash_object)
}

# The "workhorse" function has some additional parameters that are normally
#   hidden to the user.
#
flash_greedy_workhorse = function(data,
                                  Kmax = 100,
                                  f_init = NULL,
                                  var_type = c("by_column",
                                               "by_row",
                                               "constant",
                                               "zero",
                                               "kroneker"),
                                  init_fn = "udv_si",
                                  tol = 1e-2,
                                  ebnm_fn = "ebnm_pn",
                                  ebnm_param = NULL,
                                  verbose_output = "odn",
                                  nullcheck = TRUE,
                                  seed = 123,
                                  maxiter = 5000,
                                  stopping_rule = c("objective",
                                                    "loadings",
                                                    "factors",
                                                    "all_params")) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  f = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f)
  var_type = handle_var_type(match.arg(var_type), data)
  init_fn = handle_init_fn(init_fn)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, Kmax)
  verbose_output = unlist(strsplit(verbose_output, split=NULL))
  stopping_rule = match.arg(stopping_rule)

  res = add_greedy(data,
                   Kmax,
                   f,
                   var_type,
                   init_fn,
                   tol,
                   ebnm_fn,
                   ebnm_param,
                   verbose_output,
                   nullcheck,
                   seed,
                   maxiter,
                   stopping_rule)

  flash_object = construct_flash_object(data = data,
                                        fit = res$f,
                                        history = res$history,
                                        f_init = f_init)

  return(flash_object)
}

# Private function without parameter checks. Returns fit and history rather
#   than full flash object.
#
add_greedy = function(data,
                      Kmax,
                      f_init,
                      var_type,
                      init_fn,
                      tol,
                      ebnm_fn,
                      ebnm_param,
                      verbose_output,
                      nullcheck,
                      seed,
                      maxiter,
                      stopping_rule) {
  f = f_init

  prev_K = flash_get_k(f_init)
  history = list()

  if (length(verbose_output) > 0) {
    verbose_greedy_announce()
  }

  for (k in 1:Kmax) {
    if (length(verbose_output) > 0) {
      verbose_next_fl(prev_K + k, stopping_rule, tol)
    }

    old_f = f
    res = flash_r1(data,
                   f,
                   var_type,
                   init_fn,
                   tol,
                   ebnm_fn$l,
                   ebnm_param$l[[k]],
                   ebnm_fn$f,
                   ebnm_param$f[[k]],
                   verbose_output,
                   nullcheck,
                   maxiter,
                   stopping_rule)

    f = res$f
    next_history = res$history

    # Test whether the factor/loading combination is effectively zero.
    if (is_tiny_fl(f, flash_get_k(f))) {
      # Remove zero factor as long as it doesn't create an empty object.
      if (flash_get_k(f) > 1) {
        f = old_f
        next_history$zeroed_out = prev_K + k
      }
      history = c(history, list(next_history))
      break
    }

    history = c(history, list(next_history))
  }

  return(list(f = f, history = history))
}


# @title Fits a rank 1 Empirical Bayes Matrix Factorization model.
#
# @inherit flash_add_greedy
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
                    verbose_output,
                    nullcheck,
                    maxiter,
                    stopping_rule) {

  f = add_factors_from_data(data, K = 1, f_init, init_fn)

  # TODO: deal with case maxiter = 0

  opt_res = flash_optimize_single_fl(data,
                                     f,
                                     flash_get_k(f),
                                     var_type,
                                     tol,
                                     ebnm_fn_l,
                                     ebnm_param_l,
                                     ebnm_fn_f,
                                     ebnm_param_f,
                                     verbose_output,
                                     maxiter,
                                     stopping_rule)

  f = opt_res$f

  if (nullcheck) {
    null_res = perform_nullcheck(data,
                                 f,
                                 flash_get_k(f),
                                 var_type,
                                 verbose = ("n" %in% verbose_output))
    f = null_res$f
    # zeroed_out field is handled in flash_greedy_workhorse
  }

  return(list(f = f, history = opt_res$history))
}

