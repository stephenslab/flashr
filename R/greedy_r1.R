# @title Optimize a flash factor-loading combination.
#
# @details Iteratively updates factor and loading k of f (as well as
#   residual precision) to convergence of objective (used in the greedy
#   algorithm for example).
#
# @param data A flash data object.
#
# @param f A flash object.
#
# @param k The index of the factor/loading to optimize.
#
# @param var_type Type of variance structure to assume for residuals.
#
# @param tol A tolerance for the optimization.
#
# @param ebnm_fn_l Function to solve the Empirical Bayes Normal Means
#   problem (used for loadings).
#
# @param ebnm_param_l Parameters to be passed to ebnm_fn_l when
#   optimizing.
#
# @param ebnm_fn_f Function to solve the Empirical Bayes Normal Means
#   problem (used for factors).
#
# @param ebnm_param_f Parameters to be passed to ebnm_fn_f when
#   optimizing.
#
# @param verbose If TRUE, various output progress updates will be printed.
#
# @param maxiter The maximum number of iterations to perform.
#
# @param calc_obj If TRUE, convergence will be determined using the
#   objective function. Otherwise, it will be determined by looking at
#   percentage changes in EL and EF.
#
# @return An updated flash object.
#
flash_optimize_single_fl = function(data,
                                    f,
                                    k,
                                    var_type,
                                    tol,
                                    ebnm_fn_l,
                                    ebnm_param_l,
                                    ebnm_fn_f,
                                    ebnm_param_f,
                                    verbose_output,
                                    maxiter,
                                    stopping_rule) {

  if (!identical(verbose_output, "")) {
    verbose_obj_table_header(verbose_output)
  }

  obj = NULL
  obj_diff = Inf
  old_obj = -Inf
  max_chg_l = max_chg_f = Inf

  if (stopping_rule != "objective"
      || "L" %in% verbose_output || "F" %in% verbose_output) {
    res = normalize_lf(f$EL[, k], f$EF[, k])
    old_EL = res$EL
    old_EF = res$EF
  }

  R2 = flash_get_R2(data, f)

  # Expected residuals and squared residuals with factor k excluded:
  Rk = flash_get_Rk(data, f, k)
  R2k = (R2 + 2 * outer(f$EL[, k], f$EF[, k]) * Rk
         - outer(f$EL2[, k], f$EF2[, k]))

  iter = 0
  while ((iter < maxiter) &&
         !is_converged(stopping_rule, tol, obj_diff, max_chg_l, max_chg_f)) {

    iter = iter + 1
    f = flash_update_single_fl(data,
                               f,
                               k,
                               var_type,
                               ebnm_fn_l,
                               ebnm_param_l,
                               ebnm_fn_f,
                               ebnm_param_f,
                               Rk,
                               R2)

    R2 = (R2k - 2 * outer(f$EL[, k], f$EF[, k]) * Rk
          + outer(f$EL2[, k], f$EF2[, k]))

    if (stopping_rule == "objective"
        || "o" %in% verbose_output || "d" %in% verbose_output) {
      obj = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
               e_loglik_from_R2_and_tau(R2, f$tau, data$missing))
      obj_diff = obj - old_obj
      old_obj = obj
    }

    if (stopping_rule != "objective"
        || "L" %in% verbose_output || "F" %in% verbose_output) {
      res = normalize_lf(f$EL[, k], f$EF[, k])
      max_chg_l = calc_max_chg(res$EL, old_EL)
      max_chg_f = calc_max_chg(res$EF, old_EF)

      old_EL = res$EL
      old_EF = res$EF
    }

    if (!identical(verbose_output, "")) {
      verbose_obj_table_entry(verbose_output, iter, obj, obj_diff,
                              max_chg_l, max_chg_f, f$gl[k], f$gf[k])
    }
  }

  return(f)
}


# Compute the expected log-likelihood (at non-missing locations) based
#   on expected squared residuals and tau.
e_loglik_from_R2_and_tau = function(R2, tau, missing) {
  -0.5 * sum(log((2 * pi)/tau[!missing]) + tau[!missing] * R2[!missing])
}
