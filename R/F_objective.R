#' @title Get objective function from flash data and flash fit objects.
#' 
#' @param data A flash data object.
#' 
#' @param f A flash fit object.
#' 
#' @export
#' 
flash_get_objective = function(data, f) {
  return(sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) + e_loglik(data, f))
}

# @title Get expected log likelihood under current fit.
# @param data A flash data object.
# @param f A flash fit object.
e_loglik = function(data, f) {
  return(-0.5 * sum((log((2 * pi)/f$tau[!data$missing]) +
                     f$tau[!data$missing] * get_R2(data, f)[!data$missing])))
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
# @param Em Posterior mean of theta.
# 
# @param Em2 Posterior second moment of theta.
# 
NM_posterior_e_loglik = function(x, s, Et, Et2) {
  return(-0.5 * sum(log(2 * pi * s^2) + (1/s^2) * (Et2 - 2 * x * Et + x^2)))
}




