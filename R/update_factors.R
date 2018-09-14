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
# @param Rk Optionally, a matrix of residuals (excluding factor k) may
#   be passed in (for performance reasons).
#
# @param R2 A matrix of squared residuals may also be passed in.
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
                                  ebnm_param_f,
                                  Rk = NULL,
                                  R2 = NULL) {
  # Update precision:
  if (is.null(R2)) {
    R2 = flash_get_R2(data, f)
  }
  f$tau = compute_precision(R2, data$missing, var_type, data$S)

  if (is.null(Rk)) {
    Rk = flash_get_Rk(data, f, k)
  }
  f = flash_update_single_loading(data,
                                  f,
                                  k,
                                  ebnm_fn_l,
                                  ebnm_param_l,
                                  Rk,
                                  calc_obj = TRUE)
  f = flash_update_single_factor(data,
                                 f,
                                 k,
                                 ebnm_fn_f,
                                 ebnm_param_f,
                                 Rk,
                                 calc_obj = TRUE)

  return(f)
}


# @title Update a flash loading
#
# @inheritParams flash_update_single_fl
#
# @param calc_obj Specifies whether to calculate KL divergences.
#
# @return An updated flash object.
#
flash_update_single_loading = function(data,
                                       f,
                                       k,
                                       ebnm_fn,
                                       ebnm_param,
                                       Rk,
                                       calc_obj = TRUE) {
  subset = which(!f$fixl[, k])
  res = calc_update_vals(data,
                         f,
                         k,
                         subset,
                         ebnm_fn,
                         ebnm_param,
                         loadings = TRUE,
                         Rk,
                         calc_obj)

  if (!is.null(res)) {
    f$EL[subset, k] = res$EX
    f$EL2[subset, k] = res$EX2
    f$gl[[k]] = res$g
    f$ebnm_fn_l[[k]] = ebnm_fn
    f$ebnm_param_l[[k]] = ebnm_param
    if (calc_obj) {
      f$KL_l[[k]] = res$KL
    }
  }

  return(f)
}


# @title Update a flash factor
#
# @inherit flash_update_single_loading
#
flash_update_single_factor = function(data,
                                      f,
                                      k,
                                      ebnm_fn,
                                      ebnm_param,
                                      Rk,
                                      calc_obj = TRUE) {
  subset = which(!f$fixf[, k])
  res = calc_update_vals(data,
                         f,
                         k,
                         subset,
                         ebnm_fn,
                         ebnm_param,
                         loadings = FALSE,
                         Rk,
                         calc_obj)

  if (!is.null(res)) {
    f$EF[subset, k] = res$EX
    f$EF2[subset, k] = res$EX2
    f$gf[[k]] = res$g
    f$ebnm_fn_f[[k]] = ebnm_fn
    f$ebnm_param_f[[k]] = ebnm_param
    if (calc_obj) {
      f$KL_f[[k]] = res$KL
    }
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
# @return A list with elements EX, EX2, g, and KL (these are updated
#   values of either EL, EL2, gl, and KL_l or EF, EF2, gf, and KL_f).
#   If no update should be performed, returns NULL.
#
calc_update_vals = function(data,
                            f,
                            k,
                            subset,
                            ebnm_fn,
                            ebnm_param,
                            loadings,
                            Rk,
                            calc_obj = TRUE) {
  # Do not update if all elements are fixed:
  if (length(subset) == 0) {
    return(NULL)
  }

  if (loadings) {
    ebnm_args = calc_ebnm_l_args(data, f, k, subset, Rk)
  } else {
    ebnm_args = calc_ebnm_f_args(data, f, k, subset, Rk)
  }
  if (is.null(ebnm_args)) {
    return(NULL)
  }

  if (!is.null(ebnm_param$warmstart)) {
    if (ebnm_param$warmstart) {
      if (loadings && length(f$gl) >= k) {
        ebnm_param$g = f$gl[[k]]
      } else if (!loadings && length(f$gf) >= k) {
        ebnm_param$g = f$gf[[k]]
      }
    }
    ebnm_param$warmstart = NULL
  }

  a = do.call(ebnm_fn, list(ebnm_args$x, ebnm_args$s, ebnm_param))

  res = list(EX = a$postmean,
             EX2 = a$postmean2,
             g = a$fitted_g)
  if (calc_obj) {
    KL = a$penloglik - NM_posterior_e_loglik(ebnm_args$x, ebnm_args$s,
                                             a$postmean, a$postmean2)
    res = c(res, KL = KL)
  }

  return(res)
}


# @title Calculate EBNM arguments for loadings updates
#
# @inheritParams calc_update_vals
#
# @return A list with elements x and s (vectors of observations and
#   standard errors to be passed into ebnm_fn).
#
calc_ebnm_l_args = function(data, f, k, subset, Rk) {
  tau = f$tau[subset, , drop = FALSE]
  missing = data$missing[subset, , drop = FALSE]
  tau[missing] = 0
  Rk = Rk[subset, , drop = FALSE]

  s2 = 1/(tau %*% f$EF2[, k])
  if (sum(is.finite(s2)) == 0) {
    return(NULL)
  }

  x = ((Rk * tau) %*% f$EF[, k]) * s2

  # Avoid NaNs when s2 is infinite (in which case the value of x
  #   doesn't matter).
  x[is.infinite(s2)] = 0

  # If a value of s2 becomes numerically negative, set it to a
  #   small positive number.
  s = sqrt(pmax(s2, .Machine$double.eps))

  return(list(x = x, s = s))
}


# @title Calculate EBNM arguments for factor updates
#
# @inherit calc_ebnm_l_args
#
calc_ebnm_f_args = function(data, f, k, subset, Rk) {
  tau = f$tau[, subset, drop = FALSE]
  missing = data$missing[, subset, drop = FALSE]
  tau[missing] = 0
  Rk = Rk[, subset, drop = FALSE]

  s2 = 1/(t(f$EL2[, k]) %*% tau)
  if (sum(is.finite(s2)) == 0) {
    return(NULL)
  }

  x = (t(f$EL[, k]) %*% (Rk * tau)) * s2
  x[is.infinite(s2)] = 0
  s = sqrt(pmax(s2, .Machine$double.eps))

  return(list(x = x, s = s))
}
