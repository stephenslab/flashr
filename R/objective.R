#' @title Get value of objective function
#'
#' @description Get value of objective function from data and flash fit
#'   object.
#'
#' @inheritParams flash
#'
#' @param f A flash fit object.
#'
#' @export
#'
flash_get_objective = function(data, f) {
  f = handle_f(f)
  data = handle_data(data, f)

  return(sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) + e_loglik(data, f))
}


# @title Get expected log likelihood under current fit.
#
# @inheritParams flash_get_objective
#
e_loglik = function(data, f) {
  return(-0.5 * sum((log((2*pi)/f$tau[!data$missing]) +
                     f$tau[!data$missing] *
                     flash_get_R2(data, f)[!data$missing])))
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
