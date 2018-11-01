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

# Initialize a nonnegative factor/loading pair via NNLM.
#
udv_nn <- function(Y, K = 1) {
  WH <- NNLM::nnmf(Y, K, verbose = 0)
  return(list(u = WH$W, d = rep(1, K), v = t(WH$H)))
}

# Initialize a nonnegative factor (with no constraints on loadings).
#
udv_nnfactors <- function(Y, K = 1) {
  return(udv_partialnn(Y, K, nn = "factors"))
}

# Initialize a nonnegative loading vector (with no constraints on factors).
#
udv_nnloadings <- function(Y, K = 1) {
  return(udv_partialnn(Y, K, nn = "loadings"))
}

udv_partialnn <- function(Y, K = 1, nn = c("factors", "loadings")) {
  if (K > 1) {
    stop(paste("K > 1 not yet implemented for nonnegative initialization",
               "functions"))
  }
  nn <- match.arg(nn)

  udv <- flashr:::udv_si(Y, K)
  vector_to_check <- switch(nn, loadings = udv$u, factors = udv$v)
  if (is_neg_part_larger(vector_to_check)) {
    udv$u <- -udv$u
    udv$v <- -udv$v
  }

  if (identical(nn, "loadings")) {
    udv$u <- get_pos_part(udv$u)
  } else if (identical(nn, "factors")) {
    udv$v <- get_pos_part(udv$v)
  }

  return(udv)
}

is_neg_part_larger <- function(x) {
  pos_part <- x[x > 0]
  neg_part <- x[x < 0]
  return(sum(neg_part^2) > sum(pos_part^2))
}

get_pos_part <- function(x) {
  x[x < 0] <- 0
  return(x)
}
