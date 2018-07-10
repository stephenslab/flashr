# Provides functions to solve the Empirical Bayes Normal Means problem
# function must take the arguments x,s and output a list with elements
# postmean, postmean2, fitted_g and penloglik.

# @title ebnm_ash
#
# @description A wrapper to the ash function for flash.
#
# @param x A vector of observations.
#
# @param s A vector of standard errors.
#
# @param ash_param A list of parameters to be passed into ash.
#
#' @importFrom ashr ash
#'
ebnm_ash = function(x, s, ash_param) {
  a = do.call(ash,
              c(list(betahat = as.vector(x), sebetahat = as.vector(s)),
                ash_param))

  if (identical(ash_param$output, "post_sampler")) {
    return(a)
  }

  if (is.null(a$flash_data$postmean)) {
    stop(paste("ashr is not outputting flashr data in the right format.",
               "Maybe ashr needs updating to latest version?"))
  }
  return(a$flash_data)
}


# @title ebnm_pn
#
# @description A wrapper to the function
#   \code{\link[ebnm]{ebnm_point_normal}}.
#
# @param x A vector of observations.
#
# @param s A vector of standard errors.
#
# @param ebnm_param A list of parameters to be passed to the function
#   \code{ebnm_point_normal}.
#
ebnm_pn = function(x, s, ebnm_param) {
  res = do.call(ebnm::ebnm_point_normal,
                c(list(x = as.vector(x), s = as.vector(s)),
                  ebnm_param))

  if (identical(ebnm_param$output, "post_sampler")) {
    return(res)
  }

  return(list(postmean  = res$result$PosteriorMean,
              postmean2 = res$result$PosteriorMean2,
              fitted_g  = res$fitted_g,
              penloglik = res$loglik))
}


# @title EBNM using point-laplace prior, from ebnm package.
#
# @description A wrapper to the function
# \code{\link[ebnm]{ebnm_point_laplace}}.
#
# @param x A vector of observations.
#
# @param s A vector of standard errors.
#
# @param ebnm_param A list of parameters to be passed to the function
#   \code{ebnm_point_laplace}.
#
# @importFrom ebnm ebnm_point_laplace
#
ebnm_pl = function(x, s, ebnm_param) {
    res = do.call(ebnm::ebnm_point_laplace,
                  c(list(x = as.vector(x), s = as.vector(s)),
                    ebnm_param))
    return(list(postmean  = res$result$PosteriorMean,
                postmean2 = res$result$PosteriorMean2,
                fitted_g  = res$fitted_g,
                penloglik = res$loglik))
}
