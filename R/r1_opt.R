# @title Optimize a single loading and factor ('rank 1' model).
#
# @description This function iteratively optimizes the loading,
#   factor and residual precision, from residuals and their expected
#   squared values Currently the tolerance is on the changes in l and f
#   (not on the objective function)
#
# @param R An n times p matrix of data (expected residuals).
#
# @param R2 An n times p matrix of expected squared residuals.
#
# @param l_init The initial value of loading used for iterative
#   scheme (n vector).
#
# @param f_init The initial value of the factor for iterative scheme
#   (p vector).
#
# @param l2_init initial value of l2 (optional)
#
# @param f2_init initial value of f2 (optional)
#
# @param l_subset Vector of indices of l to update (default is all).
#
# @param f_subset Vector of indices of f to update (default is all).
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
# @param var_type The type of variance structure to assume.
#
# @param tol A tolerance on changes in l and f to diagnose convergence.
#
# @param calc_F Whether to compute the objective function (useful for
#   testing purposes).
#
# @param missing An n times matrix of TRUE/FALSE indicating which
#   elements of R and R2 should be considered missing (note neither R
#   nor R2 must have missing values; eg set them to 0).
#
# @param verbose If TRUE, then trace of objective function is printed.
#
# @param maxiter An upper bound on the number of iterations before
#   terminating.
#
# @param KLobj The value of the KL part of the objective for the
#   other factors not being optimized (optional, but allows objective
#   to be computed accurately).
#
# @param S Standard errors from flash data object (only used when
#   var_type = "zero")
#
# @return An updated flash object.
#
r1_opt = function(R,
                  R2,
                  l_init,
                  f_init,
                  l2_init,
                  f2_init,
                  l_subset,
                  f_subset,
                  ebnm_fn_l,
                  ebnm_param_l,
                  ebnm_fn_f,
                  ebnm_param_f,
                  var_type,
                  tol,
                  calc_F,
                  missing,
                  verbose,
                  maxiter,
                  KLobj,
                  S) {
  l = l_init
  f = f_init
  l2 = l2_init
  f2 = f2_init

  # Default initialization of l2 and f2.
  if (is.null(l2)) {
    l2 = l^2
  }
  if (is.null(f2)) {
    f2 = f^2
  }

  gl = NULL
  gf = NULL
  penloglik_l = NULL
  penloglik_f = NULL

  if (calc_F) {
    F_obj = -Inf   # Variable to store value of objective function.
    KL_f = 0
    KL_l = 0
    if (verbose) {
      message("  Iteration          Objective")
    }
  } else {
    F_obj = NULL
    KL_f = NULL
    KL_l = NULL
    if (verbose) {
      message("  Iteration         Difference")
    }
  }

  diff = Inf

  # Expected squared residuals with l and f included.
  R2new = R2 - 2 * outer(l, f) * R + outer(l2, f2)
  iter = 0

  while ((diff > tol) & (iter < maxiter)) {
    iter = iter + 1
    l_old = l
    f_old = f

    tau = compute_precision(R2new, missing, var_type, S)

    if (length(f_subset) > 0) {
      s2 = 1/(t(l2) %*% tau[, f_subset, drop = FALSE])
      if (any(is.finite(s2))) {

        # Check some finite values before proceeding.
        x = (t(l) %*% (R[, f_subset, drop = FALSE] *
                         tau[, f_subset, drop = FALSE])) * s2

        # If a value of s2 is numerically negative, set it to
        # a small positive number.
        s = sqrt(pmax(s2, .Machine$double.eps))
        ebnm_f = do.call(ebnm_fn_f, list(x, s, ebnm_param_f))
        f[f_subset] = ebnm_f$postmean
        f2[f_subset] = ebnm_f$postmean2
        gf = ebnm_f$fitted_g
        penloglik_f = ebnm_f$penloglik

        if (calc_F) {
          KL_f = ebnm_f$penloglik -
            NM_posterior_e_loglik(x, s, ebnm_f$postmean,
                                  ebnm_f$postmean2)
        }
      }
    }

    if (length(l_subset) > 0) {
      s2 = 1/(tau[l_subset, , drop = FALSE] %*% f2)
      if (any(is.finite(s2))) {

        # Check some finite values before proceeding.
        x = ((R[l_subset, , drop = FALSE] *
                tau[l_subset, , drop = FALSE]) %*% f) * s2

        # If a value of s2 is numerically negative, set it to
        # a small positive number.
        s = sqrt(pmax(s2, .Machine$double.eps))
        ebnm_l = do.call(ebnm_fn_l, list(x, s, ebnm_param_l))
        l[l_subset] = ebnm_l$postmean
        l2[l_subset] = ebnm_l$postmean2
        gl = ebnm_l$fitted_g
        penloglik_l = ebnm_l$penloglik

        if (calc_F) {
          KL_l = ebnm_l$penloglik -
            NM_posterior_e_loglik(x, s, ebnm_l$postmean,
                                  ebnm_l$postmean2)
        }
      }
    }

    R2new = R2 - 2 * outer(l, f) * R + outer(l2, f2)

    if (calc_F) {
      Fnew = KLobj + KL_l + KL_f +
        e_loglik_from_R2_and_tau(R2new, tau, missing)
      if (verbose) {
        message(sprintf("%11d", iter),
                sprintf("%19.3f", Fnew))
      }
      diff = Fnew - F_obj

      if (diff < 0) {
        warning(paste("An iteration decreased the objective.",
                      "This happens occasionally, perhaps due to",
                      "numeric reasons. You could ignore this",
                      "warning, but you might like to check out",
                      "https://github.com/stephenslab/flashr/issues/26",
                      "for more details."))
      }
      F_obj = Fnew
    } else {

      # Check convergence by percentage changes in l and f
      # normalize l and f so that f has unit norm note that this
      # messes up stored log-likelihoods etc... so not
      # recommended.
      warning("renormalization step not fully tested; be careful!")
      norm = sqrt(sum(f^2))
      f = f/norm
      f2 = f2/(norm^2)
      l = l * norm
      l2 = l2 * (norm^2)

      all_diff = abs(c(l, f)/c(l_old, f_old) - 1)
      if (all(is.nan(all_diff))) {

        # All old and new entries of l and f are zero.
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

  return(list(l = l, f = f, l2 = l2, f2 = f2, tau = tau, F_obj = F_obj,
              KL_l = KL_l, KL_f = KL_f, gl = gl, gf = gf,
              ebnm_fn_l = ebnm_fn_l, ebnm_fn_f = ebnm_fn_f,
              ebnm_param_l = ebnm_param_l, ebnm_param_f = ebnm_param_f))
}

# Put the results into f.
update_f_from_r1_opt_results = function(f, k, res) {
  f$EL[, k] = res$l
  f$EF[, k] = res$f
  f$EL2[, k] = res$l2
  f$EF2[, k] = res$f2
  f$tau = res$tau

  f$ebnm_fn_f[[k]] = res$ebnm_fn_f
  f$ebnm_fn_l[[k]] = res$ebnm_fn_l
  f$ebnm_param_f[[k]] = res$ebnm_param_f
  f$ebnm_param_l[[k]] = res$ebnm_param_l

  if (!is.null(res$gf)) {
    f$gf[[k]] = res$gf
  }
  if (!is.null(res$gl)) {
    f$gl[[k]] = res$gl
  }

  if (!is.null(res$KL_f)) {
    f$KL_f[[k]] = res$KL_f
  }
  if (!is.null(res$KL_l)) {
    f$KL_l[[k]] = res$KL_l
  }

  return(f)
}

# Compute the expected log-likelihood (at non-missing locations) based
# on expected squared residuals and tau.
e_loglik_from_R2_and_tau = function(R2, tau, missing) {
  -0.5 * sum(log((2 * pi)/tau[!missing]) + tau[!missing] * R2[!missing])
}
