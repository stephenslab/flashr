#' @title Refines a fit of the flash model to data by "backfitting".
#'
#' @description Iterates through the factors of a flash object,
#'   updating each until convergence.
#'
#' @inheritParams flash
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
                         var_type = c("by_column",
                                      "by_row",
                                      "constant",
                                      "zero",
                                      "kroneker"),
                         tol = 1e-2,
                         ebnm_fn = "ebnm_pn",
                         ebnm_param = NULL,
                         verbose = TRUE,
                         nullcheck = TRUE,
                         maxiter = 1000) {

  if (verbose) {
    verbose_output = "odn" # objective, obj diff, nullcheck
  } else {
    verbose_output = ""
  }

  this_call = match.call()

  f = flash_backfit_workhorse(data,
                              f_init,
                              kset,
                              var_type,
                              tol,
                              ebnm_fn,
                              ebnm_param,
                              verbose_output,
                              nullcheck,
                              maxiter,
                              this_call)

  return(f)
}

# The "workhorse" function has some additional parameters that are
#   normally hidden to the user.
#
flash_backfit_workhorse = function(data,
                                   f,
                                   kset = NULL,
                                   var_type = c("by_column",
                                                "by_row",
                                                "constant",
                                                "zero",
                                                "kroneker"),
                                   tol = 1e-2,
                                   ebnm_fn = "ebnm_pn",
                                   ebnm_param = NULL,
                                   verbose_output = "odn",
                                   nullcheck = TRUE,
                                   maxiter = 1000,
                                   this_call = match.call(),
                                   stopping_rule = c("objective",
                                                     "loadings",
                                                     "factors",
                                                     "all_params")) {

  flash_object = handle_f(f)
  f = get_flash_fit(flash_object)

  data = handle_data(data, f)
  kset = handle_kset(kset, f)
  var_type = handle_var_type(match.arg(var_type), data)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, length(kset))
  verbose_output = unlist(strsplit(verbose_output, split=NULL))
  stopping_rule = match.arg(stopping_rule)
  # if (is.null(this_call)) {
  #   this_call = match.call()
  # }

  if (length(verbose_output) > 0) {
    verbose_backfit_announce(length(kset), stopping_rule, tol)
  }

  if (is_max_chg_needed(stopping_rule, verbose_output)) {
    res = normalize_lf(f$EL[, kset], f$EF[, kset])
    old_EL = res$EL
    old_EF = res$EF
  }

  history = list()

  # There are two steps: first backfit (inner loop), then nullcheck
  #   (outer loop). If nullcheck removes any factors then the whole
  #   process is repeated.
  continue_outer_loop = TRUE
  while (continue_outer_loop) {
    timing <- system.time({

      if (length(verbose_output) > 0) {
        verbose_obj_table_header(verbose_output)
      }

      iter = 0
      diff = Inf
      diff_track = rep(NA, maxiter)
      obj_track = rep(NA, maxiter)

      while ((iter < maxiter) && (diff > tol)) {
        iter = iter + 1

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

        if (is_obj_needed(stopping_rule, verbose_output)) {
          obj_track[iter] = flash_get_objective(data, f)
          if (iter > 1) {
            obj_diff = obj_track[iter] - obj_track[iter - 1]
          } else {
            obj_diff = Inf
          }
        }

        if (is_max_chg_needed(stopping_rule, verbose_output)) {
          res = normalize_lf(f$EL, f$EF)
          max_chg_l = calc_max_chg(res$EL, old_EL)
          max_chg_f = calc_max_chg(res$EF, old_EF)

          old_EL = res$EL
          old_EF = res$EF
        }

        diff = calc_diff(stopping_rule, obj_diff, max_chg_l, max_chg_f)
        diff_track[iter] = diff

        if (length(verbose_output) > 0) {
          verbose_obj_table_entry(verbose_output,
                                  iter,
                                  obj_track[iter],
                                  obj_diff,
                                  max_chg_l,
                                  max_chg_f,
                                  f$gl[kset],
                                  f$gf[kset])
        }

        next_history = list(call = this_call,
                            niter = iter,
                            obj.track = obj_track[1:iter],
                            diff.track = diff_track[1:iter])
        if (!is_obj_needed(stopping_rule, verbose_output)) {
          next_history$obj_track = NULL
        }
      }

      continue_outer_loop = FALSE
      if (nullcheck) {
        # Remove factors that hurt the objective.
        kset = 1:flash_get_k(f)
        old_f = f
        res = perform_nullcheck(data, f, kset, var_type,
                                verbose = ("n" %in% verbose_output))

        # If nullcheck removes anything, start over.
        if (length(res$zeroed_out) > 0) {
          continue_outer_loop = TRUE
          f = res$f
          next_history$zeroed.out = res$zeroed_out
        }
      }
    })

    next_history$elapsed.time = timing[["elapsed"]]
    history = c(history, list(next_history))
  }

  flash_object = update_flash_object(flash_object,
                                     fit = f,
                                     obj = flash_get_objective(data, f),
                                     history = history)

  return(flash_object)
}