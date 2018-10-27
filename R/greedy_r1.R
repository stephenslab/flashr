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
  if (length(verbose_output) > 0) {
    verbose_obj_table_header(verbose_output)
  }

  if (is_max_chg_needed(stopping_rule, verbose_output)) {
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
  diff = Inf
  diff_track = rep(NA, maxiter)
  obj_track = rep(NA, maxiter)

  while ((iter < maxiter) && (diff > tol)) {
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

    if (is_obj_needed(stopping_rule, verbose_output)) {
      obj_track[iter] = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
                           e_loglik_from_R2_and_tau(R2, f$tau))
      obj_diff = calc_obj_diff(obj_track, iter)
    }

    if (is_max_chg_needed(stopping_rule, verbose_output)) {
      res = normalize_lf(f$EL[, k], f$EF[, k])
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
                              f$gl[k], f
                              $gf[k])
    }
  }


  history = list(type = "greedy",
                 kset = k,
                 niter = iter,
                 obj_track = obj_track[1:iter],
                 diff_track = diff_track[1:iter])

  if (!is_obj_needed(stopping_rule, verbose_output)) {
    history$obj_track = NULL
  }

  return(list(f = f, history = history))
}
