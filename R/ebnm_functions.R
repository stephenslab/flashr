# Provides functions to solve the Empirical Bayes Normal Means problem
# function must take the arguments x,s and output a list with elements
# postmean, postmean2, fitted_g and penloglik.

#' @title A wrapper to the ash function for flash.
#'
#' @param x Describe parameter here.
#'
#' @param s Describe parameter here.
#'
#' @param ash_param Describe parameter here.
#'
#' @importFrom utils modifyList
#' @importFrom ashr ash
#'
#' @export
#'
ebnm_ash = function(x, s, ash_param, return_sampler = FALSE) {
    if (return_sampler) {
      ash_param = modifyList(ash_param, list(output = "post_sampler"))
    }
    a = do.call(ash, c(list(betahat = as.vector(x), sebetahat = as.vector(s)), ash_param))
    if (!return_sampler) {
      if (is.null(a$flash_data$postmean)) {
        stop("ashr is not outputting flashr data in the right format. Maybe ashr needs updating to latest version?")
      }
      a = a$flash_data
    }
    return(a)
}

#' @title EBNM using point-laplace prior, from ebnm package.
#'
#' @description A wrapper to the function
#' \code{\link[ebnm]{ebnm_point_laplace}}.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard errors.
#'
#' @param ebnm_param A list of parameters to be passed to the function
#'   ebnm_point_laplace.
#'
#' @importFrom ebnm ebnm_point_laplace
#'
#' @export
#'
ebnm_pl = function(x, s, ebnm_param, return_sampler = FALSE) {
    res = do.call(ebnm_point_laplace, c(list(x = as.vector(x),
                  s = as.vector(s)), ebnm_param))
    return(list(postmean  = res$result$PosteriorMean,
                postmean2 = res$result$PosteriorMean2,
                fitted_g  = res$fitted_g,
                penloglik = res$loglik))
}

#' @title ebnm_pn
#'
#' @description A wrapper to the function
#'   \code{\link[ebnm]{ebnm_point_normal}}.
#'
#' @param x A vector of observations.
#'
#' @param s A vector of standard errors.
#'
#' @param ebnm_param A list of parameters to be passed to the function
#'   \code{ebnm_point_normal}.
#'
#' @importFrom ebnm ebnm_point_normal
#'
#' @export
#'
ebnm_pn = function(x, s, ebnm_param, return_sampler = FALSE) {
  if (return_sampler) {
    ash_param = modifyList(ebnm_param, list(output = "post_sampler"))
  }
  res = do.call(ebnm_point_normal, c(list(x = as.vector(x),s = as.vector(s)), ebnm_param))
  if (!return_sampler) {
    res = list(postmean  = res$result$PosteriorMean,
               postmean2 = res$result$PosteriorMean2,
               fitted_g  = res$fitted_g,
               penloglik = res$loglik)
  }
  return(res)
}


