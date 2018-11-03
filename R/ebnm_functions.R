# Provides functions to solve the Empirical Bayes Normal Means problem.
# Each function must take arguments x, s, ebnm_param, and output, and
# must return a list with elements postmean, postmean2, fitted_g,
# and penloglik.
#
# If sampling from the posterior is desired, then the function must also
# be able to return a sampling function when output = "post_sampler".
# This sampling function should take a single argument nsamp and return
# a matrix with nsamp rows and length(x) columns. (In other words, the
# (i,j)-entry of the matrix should correspond to the ith sample from the
# posterior for the jth element in the vector of observations.)


# @title EBNM using ash
#
# @description A wrapper to the ash function for flash.
#
# @param x A vector of observations.
#
# @param s A vector of standard errors.
#
# @param ash_param A list of parameters to be passed into ash.
#
# @param output If output = "post_sampler", then the return value is a
#   function that samples from the posterior.
#
#' @importFrom ashr ash
#'
ebnm_ash = function(x, s, ash_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ash_param$output = "post_sampler"
  } else {
    ash_param$output = "flash_data"
  }

  a = try(do.call(ash,
                  c(list(betahat = as.vector(x), sebetahat = as.vector(s)),
                    ash_param)),
          silent = TRUE)

  if (is(a, "try-error") && ash_param$optmethod != "mixIP") {
      ash_param$method = "mixIP"
      a = do.call(ash,
                  c(list(betahat = as.vector(x), sebetahat = as.vector(s)),
                    ash_param))
  }

  if (is(a, "try-error")) {
    stop(paste("Error occurred while calling ebnm_ash:", a))
  }

  if (identical(output, "post_sampler")) {
    out = a$post_sampler
  } else if (is.null(a$flash_data$postmean)) {
    stop(paste("ashr is not outputting flashr data in the right format.",
               "Maybe ashr needs updating to latest version?"))
  } else {
    out = a$flash_data
  }
  return(out)
}


# @title EBNM using point-normal prior, from ebnm package
#
# @description A wrapper to the function
#   \code{\link[ebnm]{ebnm_point_normal}}.
#
# @inheritParams ebnm_ash
#
# @param ebnm_param A list of parameters to be passed to
#   \code{ebnm_point_normal}.
#
# Do not include an @importFrom field here until ebnm is on CRAN.
#
ebnm_pn = function(x, s, ebnm_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    ebnm_param$output = "post_sampler"
  }

  res = do.call(ebnm::ebnm_point_normal,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  if (identical(output, "post_sampler")) {
    out = res$post_sampler
  } else {
    out = list(postmean = res$result$PosteriorMean,
               postmean2 = res$result$PosteriorMean2,
               fitted_g = res$fitted_g,
               penloglik = res$loglik)
  }

  return(out)
}


# @title EBNM using point-laplace prior, from ebnm package
#
# @description A wrapper to the function
# \code{\link[ebnm]{ebnm_point_laplace}}.
#
# @inheritParams ebnm_pn
#
# @param ebnm_param A list of parameters to be passed to the function
#   \code{ebnm_point_laplace}.
#
# @param output Sampling from the posterior has not yet been implemented
#   for point-laplace priors.
#
# Do not include an @importFrom field here until ebnm is on CRAN.
#
ebnm_pl = function(x, s, ebnm_param, output = NULL) {
  if (identical(output, "post_sampler")) {
    stop("Sampling from the posterior not yet implemented for ebnm_pl.")
  }

  res = do.call(ebnm::ebnm_point_laplace,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  return(list(postmean = res$result$PosteriorMean,
              postmean2 = res$result$PosteriorMean2,
              fitted_g = res$fitted_g,
              penloglik = res$loglik))
}
