# @title Update a single flash factor/loading combination (and precision).
#
# @param data A flash data object.
#
# @param f A flash fit object.
#
# @param k Index of factor/loading pair to update.
#
# @param var_type Variance structure to assume for residuals.
#
# @param ebnm_fn_l Function to solve EBNM problem (loadings updates).
#
# @param ebnm_param_l Parameters to be passed to ebnm_fn_l.
#
# @param ebnm_fn_f Function to solve EBNM problem (factor updates).
#
# @param ebnm_param_f Parameters to be passed to ebnm_fn_f.
#
# @return An updated flash object.
#
flash_update_single_fl = function(data,
                                  f,
                                  k,
                                  var_type,
                                  ebnm_fn_l,
                                  ebnm_param_l,
                                  ebnm_fn_f,
                                  ebnm_param_f) {
  f = flash_update_precision(data, f, var_type)
  f = flash_update_single_loading(data, f, k, ebnm_fn_l, ebnm_param_l)
  f = flash_update_single_factor(data, f, k, ebnm_fn_f, ebnm_param_f)

  return(f)
}


# @title Update a flash loading
#
# @inheritParams flash_update_single_fl
#
# @return An updated flash object.
#
flash_update_single_loading = function(data, f, k, ebnm_fn, ebnm_param) {
  subset = which(!f$fixl[, k])
  res = calc_update_vals(data, f, k, subset, ebnm_fn, ebnm_param,
                         loadings = TRUE)

  if (!is.null(res)) {
    f$EL[subset, k] = res$EX
    f$EL2[subset, k] = res$EX2
    f$gl[[k]] = res$g
    f$ebnm_fn_l[[k]] = ebnm_fn
    f$ebnm_param_l[[k]] = ebnm_param
    f$KL_l[[k]] = res$KL
  }

  return(f)
}


# @title Update a flash factor
#
# @inheritParams flash_update_single_fl
#
# @return An updated flash object.
#
flash_update_single_factor = function(data, f, k, ebnm_fn, ebnm_param) {
  subset = which(!f$fixf[, k])
  res = calc_update_vals(data, f, k, subset, ebnm_fn, ebnm_param,
                         loadings = FALSE)

  if (!is.null(res)) {
    f$EF[subset, k] = res$EX
    f$EF2[subset, k] = res$EX2
    f$gf[[k]] = res$g
    f$ebnm_fn_f[[k]] = ebnm_fn
    f$ebnm_param_f[[k]] = ebnm_param
    f$KL_f[[k]] = res$KL
  }

  return(f)
}


# @title Calculate updated values for factor/loading updates
#
# @inheritParams flash_update_single_loading
#
# @param subset The subset of factor or loadings entries that are not
#   considered fixed (and can thus be updated).
#
# @param loadings Should be TRUE for loadings updates and FALSE for
#   factor updates
#
# @param R Optionally, a matrix of residuals can be passed in. This
#   should speed up calculations when updates are done in parallel.
#
# @return A list with elements EX, EX2, g, and KL (these are updated
#   values of either EL, EL2, gl, and KL_l or EF, EF2, gf, and KL_f).
#   If no update should be performed, returns NULL.
#
calc_update_vals = function(data, f, k, subset,
                            ebnm_fn, ebnm_param, loadings, R = NULL) {
  # Do not update if all elements are fixed:
  if (length(subset) == 0) {
    return(NULL)
  }

  if (loadings) {
    ebnm_args = calc_ebnm_l_args(data, f, k, subset, R)
  } else {
    ebnm_args = calc_ebnm_f_args(data, f, k, subset, R)
  }
  if (is.null(ebnm_args)) {
    return(NULL)
  }
  a = do.call(ebnm_fn, list(ebnm_args$x, ebnm_args$s, ebnm_param))

  return(list(EX = a$postmean,
              EX2 = a$postmean2,
              g = a$fitted_g,
              KL = a$penloglik - NM_posterior_e_loglik(ebnm_args$x,
                                                       ebnm_args$s,
                                                       a$postmean,
                                                       a$postmean2)))
}


# @title Call calc_ebnm_args for loadings updates
#
# @inheritParams calc_update_vals
#
# @return A list with elements x and s (vectors of observations and
#   standard errors to be passed into ebnm_fn).
#
calc_ebnm_l_args = function(data, f, k, subset, R = NULL) {
  calc_ebnm_args(subset, flash_get_Rk(data, f, k, R),
                 data$missing, f$tau, f$EF, f$EF2, k)
}


# @title Call calc_ebnm_args for factor updates
#
# @inherit calc_ebnm_l_args
#
calc_ebnm_f_args = function(data, f, k, subset, R = NULL) {
  calc_ebnm_args(subset, t(flash_get_Rk(data, f, k, R)),
                 t(data$missing), t(f$tau), f$EL, f$EL2, k)
}


# @title Calculates arguments to be passed into ebnm_fn
#
# @inheritParams calc_update_vals
#
# @param Rk The matrix of residuals excluding factor/loading pair k.
#
# @param missing A matrix indicating which observations are missing.
#
# @param tau Either f$tau or t(f$tau), depending on whether loadings
#   or factors are being updated.
#
# @param EX Either f$EL or f$EF.
#
# @param EX2 Either f$EL2 or f$EF2.
#
# @return A list with elements x and s (vectors of observations and
#   standard errors to be passed into ebnm_fn).
#
calc_ebnm_args = function(subset, Rk, missing, tau, EX, EX2, k) {
  tau = tau[subset, , drop = FALSE]
  missing = missing[subset, , drop = FALSE]
  tau[missing] = 0
  Rk = Rk[subset, , drop = FALSE]

  s2 = 1/(tau %*% EX2[, k])
  if (sum(is.finite(s2)) == 0) {
    return(NULL)
  } # return NULL if there are no finite values

  x = ((Rk * tau) %*% EX[, k]) * s2
  # if a value of s2 becomes numerically negative, set it to a
  # small positive number
  s = sqrt(pmax(s2, .Machine$double.eps))

  return(list(x = x, s = s))
}


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
                                    nullcheck,
                                    tol,
                                    ebnm_fn_l,
                                    ebnm_param_l,
                                    ebnm_fn_f,
                                    ebnm_param_f,
                                    verbose,
                                    maxiter) {
  f_subset = which(!f$fixf[, k])
  l_subset = which(!f$fixl[, k])
  KLobj = (sum(unlist(f$KL_l)) + sum(unlist(f$KL_f))
           - f$KL_l[[k]] - f$KL_f[[k]])

  res = r1_opt(flash_get_Rk(data, f, k),
               flash_get_R2k(data, f, k),
               f$EL[, k],
               f$EF[, k],
               f$EL2[, k],
               f$EF2[, k],
               l_subset,
               f_subset,
               ebnm_fn_l,
               ebnm_param_l,
               ebnm_fn_f,
               ebnm_param_f,
               var_type,
               tol,
               calc_F = TRUE,
               data$missing,
               verbose,
               maxiter,
               KLobj,
               data$S)

  f = update_f_from_r1_opt_results(f, k, res)

  if (nullcheck) {
    f = perform_nullcheck(data, f, k, var_type, verbose)
  }

  return(f)
}


# @title Zeros out factors when that improves the objective.
#
# @description Sometimes zeroing out a factor can improve the
#   objective. This function iterates over factors with indices in
#   kset and checks whether zeroing it out will improve the objective;
#   if so then that factor is set to 0 (and precision is updated).
#   Returns the final flash fit object obtained when this iterative
#   process stops (ie a complete pass is performed with no factor being
#   zeroed).
#
# @param data A flash data object.
#
# @param f A flash object.
#
# @param kset The indices of the factor/loading to check.
#
# @param var_type Type of variance structure to assume for residuals.
#
# @param verbose If TRUE, various output progress updates will be
#   printed.
#
# @return A flash object.
#
perform_nullcheck = function(data, f, kset, var_type, verbose) {
  if (verbose) {
    message("Performing nullcheck...")
  }

  f_changed = TRUE  # We are going to iterate until f does not change.
  while (f_changed) {

    f_changed = FALSE
    for (k in kset) {

      f0 = flash_zero_out_factor(f, k)
      f0 = flash_update_precision(data, f0, var_type)
      F0 = flash_get_objective(data, f0)
      F1 = flash_get_objective(data, f)

      if (F0 > F1) {
        if (verbose) {
          message("  Deleting factor ", k,
                  " increases objective by ",
                  signif(F0 - F1, digits=3),
                  ". Factor zeroed out.")
        }
        f = f0
        f_changed = TRUE
      } else if (F1 > F0) {
        if (verbose) {
          message ("  Deleting factor ", k,
                   " decreases objective by ",
                   signif(F1 - F0, digits=3),
                   ". Factor retained.")
        }
      }

    }
  }
  if (verbose) {
    message("  Nullcheck complete. Objective: ",
            round(flash_get_objective(data, f), digits=3))
  }
  return(f)
}
