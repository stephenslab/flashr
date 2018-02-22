#' @title Get objective function from flash data and flash fit objects
#' @param data a flash data object
#' @param f a flash fit object
#' @export
flash_get_objective = function(data, f) {
    return(sum(unlist(f$KL_l)) + sum(unlist(f$KL_f)) + e_loglik(data, f))
}


#' @title Get expected log likelihood under current fit
#' @param data a flash data object
#' @param f a flash fit object
e_loglik = function(data, f) {
    -0.5 * sum((log((2 * pi)/f$tau[!data$missing]) + f$tau[!data$missing] * get_R2(data, f)[!data$missing]))
}


#' @title expected log likelihood for normal means model
#' @details The likelihood is for x | theta ~ N(theta, s^2);
#' The expectation is taken over the posterior on theta
#' @param x observations in normal means
#' @param s standard errors of x
#' @param Em posterior mean of theta
#' @param Em2 posterior second moment of theta
NM_posterior_e_loglik = function(x, s, Et, Et2) {
    -0.5 * sum(log(2 * pi * s^2) + (1/s^2) * (Et2 - 2 * x * Et + x^2))
}




