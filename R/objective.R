# @title Get value of objective function
#
# @description Get value of objective function from data and flash fit
#   object.
#
# @inheritParams flash
#
# @param f A flash fit object.
#
flash_get_objective = function(data, f) {
  return(sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) +
           e_loglik_from_R2_and_tau(flash_get_R2(data, f), f$tau))
}

# Compute the expected log-likelihood (at non-missing locations) based
#   on expected squared residuals and tau.
#
e_loglik_from_R2_and_tau = function(R2, tau) {
  # We assume that R2 already has NAs where data is missing.

  # If tau is stored as a scalar or as a n x 1 or 1 x p matrix (i.e., when
  #   var_type is constant, by_row, or by_column), we will rely on
  #   broadcasting. This requires conversion to a vector. Further, with a
  #   1 x p matrix, we need to transpose R2 to broadcast correctly:
  if (is.matrix(tau) && nrow(tau) == 1) {
    R2 = t(R2)
  }
  if (is.matrix(tau) && (nrow(tau) == 1 || ncol(tau) == 1)) {
    tau = as.vector(tau)
  }

  return(-0.5 * sum(log(2 * pi) - log(tau) + tau * R2, na.rm = TRUE))
}


# @title Expected log likelihood for normal means model.
#
# @description The likelihood is for x | theta ~ N(theta, s^2);
#   The expectation is taken over the posterior on theta.
#
# @param x Observations in normal means.
#
# @param s Standard errors of x.
#
# @param Et Posterior mean of theta.
#
# @param Et2 Posterior second moment of theta.
#
NM_posterior_e_loglik = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}
