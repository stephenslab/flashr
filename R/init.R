# This file contains functions related to initializing flash fit
# object.

# @title Initialize a flash fit object from the results of a factor
# analysis.
#
# @param LL The loadings, an n by K matrix.
#
# @param FF The factors, a p by K matrix.
#
# @param fixl an n by K matrix of TRUE/FALSE values indicating which
# elements of LL should be considered fixed and not changed during
# updates. Useful for including a mean factor for example.
#
# @param fixf as p by K matrix of TRUE/FALSE values; same as fixl but
# for factors FF.
#
# @return A flash fit object, with factors initialized using L and F.
#
#' @importFrom assertthat assert_that
#'
flash_init_lf = function(LL, FF, fixl = NULL, fixf = NULL) {
  assert_that(ncol(LL) == ncol(FF))

  fixl = handle_fix(fixl, LL, default_val = FALSE)
  fixf = handle_fix(fixf, FF, default_val = FALSE)

  f = list(EL = LL,
           EF = FF,
           EL2 = LL^2,
           EF2 = FF^2,
           fixl = fixl,
           fixf = fixf)

  f$gl = list()
  f$gf = list()
  f$ebnm_fn_l = list()
  f$ebnm_fn_f = list()
  f$ebnm_param_l = list()
  f$ebnm_param_f = list()
  f$KL_l = as.list(rep(0, flash_get_k(f)))
  f$KL_f = as.list(rep(0, flash_get_k(f)))
  f = c(f, list(tau = NULL))
  class(f) = "flash"

  return(f)
}


# @title Initialize an empty flash fit object.
#
# @return An empty flash fit object.
#
flash_init_null = function() {
  f = list(EL = NULL,
           EF = NULL,
           EL2 = NULL,
           EF2 = NULL,
           fixl = NULL,
           fixf = NULL,
           gl = NULL,
           gf = NULL,
           ebnm_fn_l = NULL,
           ebnm_fn_f = NULL,
           ebnm_param_l = NULL,
           ebnm_param_f = NULL,
           KL_l = NULL,
           KL_f = NULL,
           tau = NULL)
  class(f) = "flash"

  return(f)
}


# @title Initialize a flash fit object by applying a function to data.
#
# @param data a flash data object.
#
# @param init_fn An initialization function, which takes as input an
#   (n by p matrix, or flash data object) and K, a number of factors,
#   and and outputs a list with elements (u,d,v).
#
# @return A flash fit object.
#
flash_init_fn = function(data, init_fn, K=1) {
  s = do.call(init_fn, list(get_Yorig(data), K))
  f = flash_init_udv(s, K)

  l_names = rownames(get_Yorig(data))
  f_names = colnames(get_Yorig(data))
  rownames(f$EL) = rownames(f$EL2) = rownames(f$fixl) = l_names
  rownames(f$EF) = rownames(f$EF2) = rownames(f$fixf) = f_names

  return(f)
}


# @title Initialize a flash fit object from a list with elements (u,d,v).
#
# @param s List with elements (u,v,d).
#
# @param K The number of factors to use (factors \code{1:K} are used).
#
# @return A flash fit object ready for optimization.
#
flash_init_udv = function(s, K = 1) {
  s$u = as.matrix(s$u)
  s$v = as.matrix(s$v)

  # Deals with case these are vectors (K = 1).
  if (ncol(s$u) > K)
  {
    s$u = s$u[, 1:K, drop = FALSE]
  }
  if (ncol(s$v) > K) {
    s$v = s$v[, 1:K, drop = FALSE]
  }
  if (length(s$d) > K) {
    s$d = s$d[1:K]
  }

  f = flash_init_lf(t(s$d * t(s$u)), s$v)

  return(f)
}


# @title udv_si
#
# @description Provides a simple wrapper to \code{softImpute} to
#   provide a rank 1 initialization. Uses \code{type = "als"} option.
#
# @param Y An n by p matrix.
#
# @param K Number of factors to use.
#
# @return A list with components (u,d,v).
#
#' @importFrom softImpute softImpute
#'
udv_si = function(Y, K = 1) {
  suppressWarnings(
    res <- softImpute(Y, rank.max = K, type = "als", lambda = 0)
  )
  return(res)
}


# @title udv_si_svd
#
# @description provides a simple wrapper to \code{softImpute} to
#   provide a rank 1 initialization. Uses \code{type = "svd"} option.
#
# @inherit udv_si
#
#' @importFrom softImpute softImpute
#'
udv_si_svd = function(Y, K = 1) {
  suppressWarnings(
    res <- softImpute(Y, rank.max = K, type = "svd", lambda = 0)
  )
  return(res)
}


# @title udv_svd
#
# @description Provides a simple wrapper to svd.
#
# @inherit udv_si
#
udv_svd = function (Y, K = 1) {
  svd(Y, K, K)
}


# @title udv_random
#
# @description Provides a random initialization of factors.
#
# @inherit udv_si
#
# @return A list with components (u,d,v), with elements of u and v
#   i.i.d. N(0,1).
#
#' @importFrom stats rnorm
#'
udv_random = function (Y, K = 1) {
  n = nrow(Y)
  p = ncol(Y)
  return(list(u = matrix(rnorm(n * K), ncol = K),
              d = 1,
              v = matrix(rnorm(p * K), ncol = K)))
}
