#' provides functions to solve the Empirical Bayes Normal Means problem
#' function must take the arguments x,s and output a list with elements
#' postmean, postmean2, fitted_g, penloglik
#' @details a wrapper to the ash function for flash
#' @export
ebnm_ash = function(x,s,ash_param){
  ash_param = modifyList(ash_param, list(outputlevel=5))
  a = do.call(ashr::ash,c(list(betahat=as.vector(x),sebetahat=as.vector(s)),ash_param))
  if(is.null(a$flash_data$postmean)){stop("ashr is not outputting flashr data. Maybe your version of ashr needs updating to latest version?")}
  return(a$flash_data)
}

# EBNM using point-laplace prior, from ebnm package
#' @title ebnm_pl
#' @details a wrapper to the function ebnm::point_laplace
#' @param x a vector of observations
#' @param s a vector of standard errors
#' @param ebnm_param a list of parameters to be passed to the function ebnm_point_laplace
#' @export
ebnm_pl = function(x,s,ebnm_param){
  res = do.call(ebnm::ebnm_point_laplace,c(list(x=as.vector(x),s=as.vector(s)),ebnm_param) )
  return(list(postmean = res$result$PosteriorMean, postmean2 = res$result$PosteriorMean2, fitted_g = res$fitted_g, penloglik = res$loglik))
}

#' @title ebnm_pn
#' @details a wrapper to the function ebnm::point_normal
#' @param x a vector of observations
#' @param s a vector of standard errors
#' @param ebnm_param a list of parameters to be passed to the function ebnm_point_normal
#' @export
ebnm_pn = function(x,s,ebnm_param){
  res = do.call(ebnm::ebnm_point_normal,c(list(x=as.vector(x),s=as.vector(s)),ebnm_param) )
  return(list(postmean = res$result$PosteriorMean, postmean2 = res$result$PosteriorMean2, fitted_g = res$fitted_g, penloglik = res$loglik))
}


