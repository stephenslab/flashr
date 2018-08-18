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

  f = flash_backfit_workhorse(data,
                              f_init,
                              kset,
                              var_type,
                              tol,
                              ebnm_fn,
                              ebnm_param,
                              verbose_output,
                              nullcheck,
                              maxiter)

  return(f)
}

# This function has some additional parameters that are normally
#   hidden to the user.
#
flash_backfit_workhorse = function(data,
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
                                   verbose_output = "odn",
                                   nullcheck = TRUE,
                                   maxiter = 1000,
                                   stopping_rule = c("objective",
                                                     "loadings",
                                                     "factors",
                                                     "all_params")) {

  f = handle_f(f_init)
  data = handle_data(data, f)
  kset = handle_kset(kset, f)
  var_type = handle_var_type(match.arg(var_type), data)
  ebnm_fn = handle_ebnm_fn(ebnm_fn)
  ebnm_param = handle_ebnm_param(ebnm_param, ebnm_fn, length(kset))
  stopping_rule = match.arg(stopping_rule)

  if (!identical(verbose_output, "")) {
    verbose_output = unlist(strsplit(verbose_output, split=NULL))
    verbose_backfit_announce(length(kset), stopping_rule, tol)
    verbose_obj_table_header(verbose_output)
  }

  obj = NULL
  obj_diff = Inf
  old_obj = -Inf
  max_chg_l = max_chg_f = Inf

  if (stopping_rule != "objective"
      || "L" %in% verbose_output || "F" %in% verbose_output) {
    norms = sqrt(colSums(f$EL[, kset, drop = FALSE]^2))
    old_EL = as.vector(sweep(f$EL[, kset, drop = FALSE], 2, norms, `/`))
    old_EF = as.vector(sweep(f$EF[, kset, drop = FALSE], 2, norms, `*`))
  }

  # There are two steps: first backfit (inner loop), then nullcheck
  #   (outer loop). If nullcheck removes any factors then the whole
  #   process is repeated.
  continue_outer_loop = TRUE
  while (continue_outer_loop) {

    iter = 0
    while ((iter < maxiter) &&
           !is_converged(stopping_rule, tol, obj_diff, max_chg_l, max_chg_f)) {

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

      if (stopping_rule == "objective"
          || "o" %in% verbose_output || "d" %in% verbose_output) {
        obj = flash_get_objective(data, f)
        obj_diff = obj - old_obj
        old_obj = obj
      }

      if (stopping_rule != "objective"
          || "L" %in% verbose_output || "F" %in% verbose_output) {
        norms = sqrt(colSums(f$EL[, kset, drop = FALSE]^2))
        EL = as.vector(sweep(f$EL[, kset, drop = FALSE], 2, norms, `/`))
        EF = as.vector(sweep(f$EF[, kset, drop = FALSE], 2, norms, `*`))
        max_chg_l = calc_max_chg(EL, old_EL)
        max_chg_f = calc_max_chg(EF, old_EF)

        old_EL = EL
        old_EF = EF
      }

      if (!identical(verbose_output, "")) {
        verbose_obj_table_entry(verbose_output, iter, obj, obj_diff,
                                max_chg_l, max_chg_f, f$gl[kset], f$gf[kset])
      }
    }

    continue_outer_loop = FALSE
    if (nullcheck) {
      # Remove factors that hurt the objective.
      kset = 1:flash_get_k(f)
      old_f = f
      f = perform_nullcheck(data, f, kset, var_type,
                            verbose = ("n" %in% verbose_output))

      # If nullcheck removes anything, start over.
      continue_outer_loop = !identical(old_f, f)
    }
  }

  return(f)
}
