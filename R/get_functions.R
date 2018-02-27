# Functions for extracting useful information about the object.

#' @title Return the estimated LF' matrix.
#' 
#' @param f A flash fit object.
#' 
#' @return The estimated value of LF', an n by p matrix.
#' 
#' @export
#' 
flash_get_lf = function(f) {
    if (is.null(f$EL)) {
        return(NULL)
    }
    return(f$EL %*% t(f$EF))
}

#' @title flash_get_ldf
#' 
#' @param f A flash fit object.
#' 
#' @param k Indices of loadings/factors to be returned.
#' 
#' @param drop_zero_factors Flag whether to remove any factor/loadings
#'   that are zero.
#' 
#' @details Returns standardized loadings, factors, and weights from a
#' flash object.
#' 
#' @return A list with the following elements:
#' \itemize{
#' \item{l} A matrix whose columns contain the standardized loadings
#'   (ie norm 1).
#' \item{d} A vector of weights (analogous to the singular values in
#'   an svd).
#' \item{f} A matrix whose columns contain the standardized factors
#'   (ie norm 1).}
#' These are analogous to the u, d and v returned by \code{svd}, but
#' columns of l and f are not orthogonal.
#' 
#' @export
#' 
flash_get_ldf = function(f, k = NULL, drop_zero_factors =TRUE) {
  if (is.null(k)) {
    k = 1:flash_get_k(f)
  }
  ll = f$EL[, k, drop=FALSE]
  ff = f$EF[, k, drop=FALSE]
  d= sqrt(colSums(ll^2) * colSums(ff^2))

  ll = scale(ll, scale=sqrt(colSums(ll^2)), center=FALSE)
  ff = scale(ff, scale=sqrt(colSums(ff^2)), center=FALSE)

  if(drop_zero_factors){
    ll = ll[,d!=0,drop=FALSE]
    ff = ff[,d!=0,drop=FALSE]
    d = d[d!=0,drop=FALSE]
  }

  list(d=d,
       l = ll,
       f = ff)
}

# @title Get the residuals from a flash data and fit object,
#   excluding factor k.
flash_get_Rk = function(data, f, k) {
    if (flash_get_k(f) < k) {
        stop("factor k does not exist")
    }
    return(data$Y - f$EL[, -k, drop = FALSE] %*% t(f$EF[, -k, drop = FALSE]))
}

# @title Get the residuals from data and a flash fit object.
flash_get_R = function(data, f) {
    if (is.null(f$EL))
        {
            return(data$Y)
        }  # if f is null, return Y
 else {
        return(data$Y - flash_get_lf(f))
    }
}

# @title Get the residuals from data and a flash fit object, with
#   missing data as in original.
flash_get_R_withmissing = function(data, f) {
    if (is.null(f$EL))
        {
            return(get_Yorig(data))
        }  # if f is null, return Y
 else {
        return(get_Yorig(data) - flash_get_lf(f))
    }
}

# @title Get the expected squared residuals from a flash data and fit
# object.
flash_get_R2 = function(data, f) {
    if (is.null(f$EL)) {
        return(data$Y^2)
    } else {
        LF = f$EL %*% t(f$EF)
        return((data$Y - LF)^2 + f$EL2 %*% t(f$EF2) - f$EL^2 %*% t(f$EF^2))
    }
}

# @title Get the expected squared residuals from a flash data and fit
#   object excluding factor k.
flash_get_R2k = function(data, f, k) {
    if (is.null(f$EL)) {
        return(data$Y^2)
    } else {
        return(flash_get_Rk(data, f, k)^2 +
               f$EL2[, -k, drop = FALSE] %*% t(f$EF2[, -k, drop = FALSE]) -
               f$EL[, -k, drop = FALSE]^2 %*%
               t(f$EF[, -k, drop = FALSE]^2))
    }
}

# @title is_tiny_fl
# @details Checks whether kth factor/loading combination is tiny.
is_tiny_fl = function(f, k, tol = 1e-08) {
    return(sum(f$EL[, k]^2) * sum(f$EF[, k]^2) < tol)
}

flash_get_l = function(f) {
  f$EL
}

flash_get_f = function(f) {
  f$EF
}

#' @title Get number of factors in a fit object.
#' 
#' @param f A flash fit object.
#' 
#' @details Returns the number of factors in a flash fit.
#' 
#' @export
#' 
flash_get_k = function(f) {
    k = ncol(f$EL)
    if (is.null(k)) {
        return(0)
    } else {
        return(k)
    }
}

flash_get_n = function(f) {
    nrow(f$EL)
}

flash_get_p = function(f) {
  nrow(f$EF)
}

#' @title flash_get_pve
#' 
#' @description Returns the factor contributions ('proportion of
#' variance explained') by each factor/loading combination in flash
#' fit f. Because the factors are not required to be orthogonal this
#' should be interpreted loosely: eg PVE could total more than 1.
#'
#' @params f Description of input argument goes here.
#' 
#' @export
#' 
flash_get_pve = function(f) {
    s = (flash_get_ldf(f)$d)^2
    s/(sum(s) + sum(1/f$tau))
    # wei's version:
    #
    #  sapply(seq(1,K),function(x){
    #    sum(f$EL2[,x] %*% t(f$EF2[,x])) }))/sum(Y^2)
    #
}

flash_get_conv_criteria = function(data, f) {
    flash_get_objective(data, f)
}
