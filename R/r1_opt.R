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
# @param nullcheck Flag whether to check, after running hill-climbing
#   updates, whether the achieved optimum is better than setting factor
#   to 0. If this check is performed and fails then the factor will be
#   set to 0 in the returned fit.
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
                                    verbose,
                                    maxiter,
                                    calc_obj = TRUE) {

  if (verbose) {
    if (calc_obj) {
      message("  Iteration          Objective")
    } else {
      message("  Iteration         Difference")
    }
  }

  diff = Inf
  old_obj = -Inf

  R2 = flash_get_R2(data, f)

  # Expected residuals and squared residuals with factor k excluded:
  Rk = flash_get_Rk(data, f, k)
  R2k = (R2 + 2 * outer(f$EL[, k], f$EF[, k]) * Rk
         - outer(f$EL2[, k], f$EF2[, k]))

  iter = 0
  while ((diff > tol) & (iter < maxiter)) {
    iter = iter + 1

    old_vals = c(f$EL[, k], f$EF[, k])

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

    if (calc_obj) {
      # Check convergence by increase in objective function.
      obj = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
               e_loglik_from_R2_and_tau(R2, f$tau, data$missing))
      diff = obj - old_obj
      old_obj = obj

      if (diff < 0) {
        display_obj_decr_warning()
      }
      if (verbose) {
        message(sprintf("%11d", iter), sprintf("%19.3f", obj))
      }
    } else {
      # Check convergence by percentage changes in EL and EF.
      #   Normalize EL and EF so that EF has unit norm. Note that this
      #   messes up stored log-likelihoods etc... so not recommended.
      warning("Renormalization step not fully tested; be careful!")

      norm = sqrt(sum(f$EF[, k]^2))
      f$EF[, k] = f$EF[, k] / norm
      f$EF2[, k] = f$EF2[, k] / (norm^2)
      f$EL[, k] = f$EL[, k] * norm
      f$EL2[, k] = f$EL2[, k] * (norm^2)

      all_diff = abs(c(f$EL[, k], f$EF[, k])/old_vals - 1)
      if (all(is.nan(all_diff))) {
        # All old and new entries of EL and EF are zero.
        diff = 0
      } else {
        # Ignore entries where both old and new values are zero.
        diff = max(all_diff[!is.nan(all_diff)])
      }

      if (verbose) {
        message(sprintf("%11d", iter),
                sprintf("%19.3f", diff))
      }
    }
  }

  return(f)
}


# Compute the expected log-likelihood (at non-missing locations) based
#   on expected squared residuals and tau.
e_loglik_from_R2_and_tau = function(R2, tau, missing) {
  -0.5 * sum(log((2 * pi)/tau[!missing]) + tau[!missing] * R2[!missing])
}


# Warning to be displayed whenever the objective decreases.
display_obj_decrease_warning = function() {
  warning(paste("An iteration decreased the objective.",
                "This happens occasionally, perhaps due to",
                "numeric reasons. You could ignore this",
                "warning, but you might like to check out",
                "https://github.com/stephenslab/flashr/issues/26",
                "for more details."))
}
