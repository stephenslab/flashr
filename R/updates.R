# @title Update a flash loading
#
# @description Updates loading k of f to increase the objective F.
#   Updates only the loading, once (not the factor).
#
# @inheritParams flash_update_single_fl
#
# @return an updated flash object
#
flash_update_single_loading = function(data,
                                       f,
                                       k,
                                       ebnm_fn = ebnm_pn,
                                       ebnm_param = flash_default_ebnm_param(ebnm_fn),
                                       gl = NULL,
                                       fixgl = FALSE,
                                       return_sampler = F) {
    if (!is.null(gl)) {
      ebnm_param = modifyList(ebnm_param, list(fixg=fixgl, g=gl))
    }

    subset = which(!f$fixl[, k])  # check which elements are not fixed
    if (length(subset) > 0) {
        # and only do the update if some elements are not fixed

        tau = f$tau[subset, , drop = FALSE]

        if (data$anyNA)
            {
                tau = tau * (!data$missing[subset, ])
            }  # set missing values to have precision 0

        s2 = 1/(tau %*% f$EF2[, k])
        if (sum(is.finite(s2)) > 0) {
            # check some finite values before proceeding
            Rk = flash_get_Rk(data, f, k)[subset, ]  #residuals excluding factor k
            x = ((Rk * tau) %*% f$EF[,k]) * s2
            # if a value of s2 becomes numerically negative, set it to a small positive number
            s = sqrt(pmax(s2, .Machine$double.eps))
            a = ebnm_fn(x, s, ebnm_param, return_sampler)
            if (return_sampler) {
                if (is.null(a$post_sampler)) {
                    stop("No sampler implemented for that ebnm function.")
                }
                return(sampler(f$fixl[, k], a$post_sampler, f$EL[f$fixl[, k], k]))
            }

            f$EL[subset, k] = a$postmean
            f$EL2[subset, k] = a$postmean2
            f$gl[[k]] = a$fitted_g
            f$ebnm_param_l[[k]] = ebnm_param
            f$KL_l[[k]] = a$penloglik - NM_posterior_e_loglik(x, s, a$postmean, a$postmean2)
            f$penloglik_l[[k]] = a$penloglik
        } else if (return_sampler) { # if all else fails, sample values at their expectation
            return(sampler(rep(TRUE, length(f$EL[, k])), NULL, f$EL[, k]))
        }
    }
    return(f)
}


# @title  Update a flash factor
#
# @description Updates factor k of f to increase the objective F.
#   Updates only the factor, once (not the loading).
#
# @inheritParams flash_update_single_fl
#
# @return an updated flash object
#
flash_update_single_factor = function(data,
                                      f,
                                      k,
                                      ebnm_fn = ebnm_pn,
                                      ebnm_param = flash_default_ebnm_param(ebnm_fn),
                                      gf = NULL,
                                      fixgf = FALSE,
                                      return_sampler = F) {
    if (!is.null(gf)) {
      ebnm_param = modifyList(ebnm_param, list(fixg=fixgf, g=gf))
    }

    subset = which(!f$fixf[, k])  # check which elements are not fixed
    if (length(subset) > 0) {
        # and only do the update if some elements are not fixed

        tau = f$tau[, subset, drop = FALSE]
        if (data$anyNA)
            {
                tau = tau * (!data$missing[, subset])
            }  # set missing values to have precision 0

        s2 = 1/(t(tau) %*% f$EL2[, k])
        if (sum(is.finite(s2)) > 0) {
            # check some finite values before proceeding
            Rk = flash_get_Rk(data, f, k)[, subset]  #residuals excluding factor k
            x = (t(Rk * tau) %*% f$EL[, k]) * s2
            # if a value of s2 becomes numerically negative, set it to a small positive number
            s = sqrt(pmax(s2, .Machine$double.eps))
            a = ebnm_fn(x, s, ebnm_param, return_sampler)
            if (return_sampler) {
                if (is.null(a$post_sampler)) {
                    stop("No sampler implemented for that ebnm function.")
                }
                return(sampler(f$fixf[, k], a$post_sampler, f$EF[f$fixl[, k], k]))
            }

            f$EF[subset, k] = a$postmean
            f$EF2[subset, k] = a$postmean2
            f$gf[[k]] = a$fitted_g
            f$ebnm_param_f[[k]] = ebnm_param
            f$KL_f[[k]] = a$penloglik - NM_posterior_e_loglik(x, s, a$postmean, a$postmean2)
            f$penloglik_f[[k]] = a$penloglik
        }
    } else if (return_sampler) { # if all else fails, sample values at their expectation
        return(sampler(rep(TRUE, length(f$EF[, k])), NULL, f$EF[, k]))
    }
    return(f)
}

# @title Update a single flash factor-loading combination (and precision).
#
# @param data a flash data object
# @param f a flash fit object
# @param k the index of the loading to update
# @param ebnm_fn function to solve the Empirical Bayes normal means problem
# @param ebnm_param parameters to be passed to ebnm_fn when optimizing
# @param gl if nonnull, fixes the prior on the loading
# @param gf if nonnull, fixes the prior on the factor
#
# @return an updated flash object
#
flash_update_single_fl = function(data,
                                  f,
                                  k,
                                  var_type,
                                  ebnm_fn = ebnm_pn,
                                  ebnm_param = flash_default_ebnm_param(ebnm_fn),
                                  gl = NULL,
                                  fixgl = FALSE,
                                  gf = NULL,
                                  fixgf = FALSE) {
    f = flash_update_precision(data, f, var_type)
    f = flash_update_single_factor(data, f, k, ebnm_fn, ebnm_param, gf, fixgf)
    f = flash_update_single_loading(data, f, k, ebnm_fn, ebnm_param, gl, fixgl)
    return(f)
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
# @param ebnm_fn Function to solve the Empirical Bayes normal means
#   problem.
#
# @param gl Passed into ebnm_fn as parameter g.
#
# @param fixgl Passed into ebnm_fn as parameter fixg.
#
# @param gf Passed into ebnm_fn as parameter g.
#
# @param fixgf Passed into ebnm_fn as parameter fixg.
#
# @param ebnm_param Parameters to be passed to ebnm_fn when optimizing.
#
# @param verbose If TRUE, various output progress updates will be printed.
#
# @return An updated flash object.
#
flash_optimize_single_fl = function(data,
                                    f,
                                    k,
                                    var_type,
                                    nullcheck = TRUE,
                                    tol = 0.01,
                                    ebnm_fn = ebnm_pn,
                                    ebnm_param = flash_default_ebnm_param(ebnm_fn),
                                    gl = NULL,
                                    fixgl = FALSE,
                                    gf = NULL,
                                    fixgf = FALSE,
                                    verbose = FALSE) {
    f_subset = which(!f$fixf[, k])
    l_subset = which(!f$fixl[, k])
    KLobj = sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) - f$KL_l[[k]] - f$KL_f[[k]]

    res = r1_opt(flash_get_Rk(data, f, k),
                 flash_get_R2k(data, f, k),
                 f$EL[, k],
                 f$EF[, k],
                 f$EL2[, k],
                 f$EF2[, k],
                 l_subset,
                 f_subset,
                 ebnm_fn,
                 ebnm_param,
                 gl,
                 fixgl,
                 gf,
                 fixgf,
                 var_type,
                 tol,
                 calc_F = TRUE,
                 missing = data$missing,
                 verbose = verbose,
                 KLobj = KLobj,
                 S = data$S)

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

    f_changed = TRUE  #we are going to iterate until f does not change
    while (f_changed) {

        f_changed = FALSE
        for (k in kset) {

            f0 = flash_zero_out_factor(data, f, k)
            f0 = flash_update_precision(data, f0, var_type)
            F0 = flash_get_objective(data, f0)
            F1 = flash_get_objective(data, f)

            if (verbose) {
                message("performing nullcheck")
                message("objective from deleting factor:", F0)
                message("objective from keeping factor:", F1)
            }

            if (F0 > F1) {
                if (verbose) {
                  message("factor zeroed out")
                }
                f = f0
                f_changed = TRUE
            }

        }
    }
    if (verbose) {
        message("nullcheck complete, objective:", flash_get_objective(data, f))
    }
    return(f)
}


