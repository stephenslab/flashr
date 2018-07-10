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
    if (is.null(fixl)) {
        fixl = matrix(FALSE, ncol = ncol(LL), nrow = nrow(LL))
    }
    if (is.null(fixf)) {
        fixf = matrix(FALSE, ncol = ncol(FF), nrow = nrow(FF))
    }

    f = list(EL = LL, EF = FF, EL2 = LL^2, EF2 = FF^2, fixl = fixl,
             fixf = fixf)

    f$gl = list()
    f$gf = list()
    f$ebnm_fn_l = list()
    f$ebnm_fn_f = list()
    f$ebnm_param_l = list()
    f$ebnm_param_f = list()
    f$KL_l = as.list(rep(0, flash_get_k(f)))
    f$KL_f = as.list(rep(0, flash_get_k(f))) # KL divergences for each l and f.
    f = c(f, list(tau = NULL))
    class(f) = "flash"
    return(f)
}

# @title Initialize an empty flash fit object.
#
# @return An empty flash fit object.
#
flash_init_null = function() {
    f = list(EL = NULL, EF = NULL, fixl = NULL, fixf = NULL, EL2 = NULL, EF2 = NULL, gl = NULL, gf = NULL, ebnm_param_l = NULL,
        ebnm_param_f = NULL, KL_l = NULL, KL_f = NULL, tau = NULL)
    class(f) = "flash"
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
  suppressWarnings(res <- softImpute(Y, rank.max = K,
                                     type = "als", lambda = 0))
  return(res)
}

# @title udv_si_svd
#
# @description provides a simple wrapper to \code{softImpute} to
#   provide a rank 1 initialization. Uses \code{type = "svd"} option.
#
# @param Y an n by p matrix
#
# @param K number of factors to use
#
# @return A list with components (u,d,v)
#
#' @importFrom softImpute softImpute
#'
udv_si_svd = function(Y, K = 1) {
  suppressWarnings(res <- softImpute(Y,rank.max = K,type = "svd",lambda = 0))
  return(res)
}

# @title udv_svd
#
# @description Provides a simple wrapper to svd.
#
# @param Y An n by p matrix.
#
# @param K Number of factors to use.
#
# @return A list with components (u,d,v).
#
udv_svd = function (Y, K = 1) {
  svd(Y, K, K)
}

# @title udv_random
#
# @description Provides a random initialization of factors.
#
# @param Y An n by p matrix.
#
# @param K Number of factors.
#
# @return A list with components (u,d,v), with elements of u and v
#   i.i.d. N(0,1).
#
#' @importFrom stats rnorm
#'
udv_random = function (Y, K = 1) {
  n = nrow(Y)
  p = ncol(Y)
  return(list(u = matrix(rnorm(n * K), ncol = K), d = 1,
              v = matrix(rnorm(p * K), ncol = K)))
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
