# TODO: fix up public function

#' @title Add a set of fixed loadings to a flash fit object
#'
#' @inheritParams flash
#'
#' @param LL The loadings, an n by K matrix. Missing values will be
#'   initialized by the mean of the relevant column (but will generally
#'   be re-estimated when refitting the model).
#'
#' @param fixl An n by K matrix of \code{TRUE}/\code{FALSE} values
#'   indicating which elements of \code{LL} should be considered fixed
#'   and not changed during updates.  The default is to fix all
#'   non-missing values, so missing values will be updated when the
#'   flash object is updated.
#'
#' @param ... Additional parameters to be passed to \code{flash_backfit}.
#'   Note that \code{nullcheck} defaults to \code{FALSE} here.
#'
#' @return A flash object, with loadings initialized from \code{LL},
#'   and corresponding factors initialized to zero.
#'
#' @export
#'
flash_add_fixed_loadings = function(data,
                                    LL,
                                    f_init = NULL,
                                    fixl = NULL,
                                    var_type,
                                    ebnm_fn,
                                    ebnm_param,
                                    stopping_rule,
                                    tol,
                                    verbose_output,
                                    nullcheck,
                                    maxiter) {
  f = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f)
  LL = handle_LL(LL, expected_nrow = flash_get_n(f))
  fixl = handle_fix(fixl, LL, default_val = TRUE)
  var_type = handle_var_type(match.arg(var_type), data)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, n)
  verbose_output = unlist(strsplit(verbose_output, split = NULL))
  stopping_rule = match.arg(stopping_rule)

  f = add_fixed_loadings(data,
                         LL,
                         f,
                         fixl,
                         var_type,
                         ebnm_fn,
                         ebnm_param,
                         stopping_rule,
                         tol,
                         verbose_output,
                         nullcheck,
                         maxiter)

  flash_object = construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init,
                                        compute_obj = backfit)

  return(flash_object)
}

# Private function without parameter checks. Returns fit and history rather
#   than full flash object.
#
add_fixed_loadings = function(data,
                              LL,
                              f_init,
                              fixl,
                              var_type,
                              ebnm_fn,
                              ebnm_param,
                              stopping_rule,
                              tol,
                              verbose_output,
                              nullcheck,
                              maxiter,
                              init_tol = 1e-3) {
  f = f_init

  prev_K = flash_get_k(f_init)
  history = list()

  for (k in 1:ncol(LL)) {
    if (length(verbose_output) > 0) {
      verbose_next_fl(prev_K + k, stopping_rule, tol)
    }

    # Use matrix of residuals with NAs set to zero:
    R = flash_get_R(data, f)

    ll = LL[, k, drop = FALSE]

    if (all(fixl)) {
      ff = crossprod(R, ll) / as.numeric(crossprod(ll))
    } else {
      # First impute the missing data:
      missing = is.na(ll)
      if (all(ll[!missing] == 0)) { # sparse case
        ll[missing] = 1
      } else { # not sparse
        ll[missing] = mean(ll[!missing])
      }

      # Then estimate ll and ff via alternating least squares:
      ff = crossprod(R, ll) / as.numeric(crossprod(ll))
      oldvals = c(ll[!fixl], ff)
      rel_chg = Inf
      while (rel_chg > init_tol) {
        ll[!fixl] = tcrossprod(R[!fixl, ], t(ff)) / as.numeric(crossprod(ff))
        ff = crossprod(R, ll) / as.numeric(crossprod(ll))
        newvals = c(ll[!fixl], ff)
        rel_chg = max(abs(1 - newvals / oldvals))
        oldvals = newvals
      }
    }

    f2 = flash_init_lf(ll, ff, fixl = fixl)
    f = flash_combine(f, f2)

    if (maxiter > 0) {
      opt_res = flash_optimize_single_fl(data,
                                         f,
                                         prev_K + k,
                                         var_type,
                                         tol,
                                         ebnm_fn$l,
                                         ebnm_param$l[[prev_K + k]],
                                         ebnm_fn$f,
                                         ebnm_param$f[[prev_K + k]],
                                         verbose_output,
                                         maxiter,
                                         stopping_rule)

      f = opt_res$f
      opt_res$history$type = "add-fixed"
      history = c(history, list(opt_res$history))
    }
  }

  return(list(f = f, history = history))
}


#' @title Add a set of fixed factors to a flash fit object
#'
#' @inheritParams flash_add_fixed_loadings
#'
#' @param FF The factors, a p vector or p by K matrix. Missing values
#'   will be initialized by the mean of the relevant column (but will
#'   generally be re-estimated when refitting the model).
#'
#' @param fixf A p by K matrix of of \code{TRUE}/\code{FALSE} values
#'   indicating which elements of \code{FF} should be considered fixed
#'   and not changed during updates.  The default is to fix all
#'   non-missing values, so missing values will be updated when the
#'   flash object is updated.
#'
#' @return A flash object, with factors initialized from \code{FF},
#'   and corresponding loadings initialized to zero.
#'
#' @export
#'
flash_add_fixed_factors = function(data,
                                   FF,
                                   f_init = NULL,
                                   fixf = NULL,
                                   init_fn = "udv_si",
                                   backfit = TRUE,
                                   ...) {
  f = handle_f(f_init)
  data = handle_data(data, f)
  # FF, fixf, and init_fn are handled by flash_add_fixed_loadings

  flash_object = flash_add_fixed_loadings(flash_transpose_data(data),
                                          FF,
                                          flash_transpose(f),
                                          fixf,
                                          init_fn,
                                          backfit,
                                          ...)
  f = flash_transpose(get_flash_fit(flash_object))
  history = get_flash_fit_history(flash_object)

  flash_object = construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init,
                                        compute_obj = backfit)

  return(flash_object)
}

# Private function without parameter checks. Returns fit and history rather
#   than full flash object.
#
add_fixed_factors = function(data,
                             FF,
                             f_init,
                             fixf,
                             var_type,
                             ebnm_fn,
                             ebnm_param,
                             stopping_rule,
                             tol,
                             verbose_output,
                             nullcheck,
                             maxiter,
                             init_tol = 1e-3) {
  if (stopping_rule == "loadings") {
    stopping_rule = "factors"
  } else if (stopping_rule == "factors") {
    stopping_rule = "loadings"
  }

  res = add_fixed_loadings(flash_transpose_data(data),
                           FF,
                           flash_transpose(f_init),
                           fixf,
                           var_type,
                           list(l = ebnm_fn$f, f = ebnm_fn$l),
                           list(l = ebnm_param$f, f = ebnm_param$l),
                           stopping_rule,
                           tol,
                           verbose_output,
                           nullcheck,
                           maxiter,
                           init_tol)
  res$f = flash_transpose(res$f)

  return(res)
}


# @title Add factor/loading pairs to a flash object
#
# @description Adds specified factor/loading pairs to a flash object.
#
# @inheritParams flash
#
# @param LL The loadings, an n by K matrix.
#
# @param FF The factors, a p by K matrix.
#
# @param fixl An n by K matrix of \code{TRUE}/\code{FALSE} values
#   indicating which elements of \code{LL} should be considered fixed
#   and not changed during updates. Useful for including a mean factor
#   for example.
#
# @param fixf A p by K matrix of \code{TRUE}/\code{FALSE} values; same
#   as \code{fixl} but for factors \code{FF}.
#
# @return A flash fit object, with additional factors initialized
#   using \code{LL} and \code{FF}.
#
flash_add_lf = function(data,
                        LL,
                        FF,
                        f_init = NULL,
                        fixl = NULL,
                        fixf = NULL) {
  f_init = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f_init)
  LL = handle_LL(LL, expected_nrow = flash_get_n(f_init))
  FF = handle_LL(FF, expected_nrow = flash_get_p(f_init))
  # fixl and fixf are handled by flash_init_lf

  f2 = flash_init_lf(LL, FF, fixl, fixf)
  f = flash_combine(f_init, f2)

  return(f)
}
