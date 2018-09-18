# Functions for extracting useful information about the object.

# @title Get fitted values.
#
# @description Return the estimated LF' matrix.
#
# @inheritParams flash_get_ldf
#
# @return The estimated value of LF', an n by p matrix.
#
flash_get_fitted_values = function(f) {
  f = handle_f(f, allow_null = FALSE)
  if (is.null(f$EL)) {
    return(NULL)
  }
  return(f$EL %*% t(f$EF))
}


# @title Get LDF from a flash object
#
# @description Returns standardized loadings, factors, and weights from
#   a flash object.
#
# @param f A flash fit object.
#
# @param kset Indices of loadings/factors to be returned.
#
# @param drop_zero_factors If \code{TRUE}, then any factor/loadings
#   that are zero will be removed.
#
# @return A list with the following elements. These are analogous to
#   the \code{u}, \code{d} and \code{v} returned by \code{svd}, but
#   the columns of \code{l} and \code{f} are not orthogonal.
#
#   \item{\code{l}}{A matrix whose columns contain the standardized
#     loadings (i.e., with norm 1).}
#
#   \item{\code{d}}{A vector of weights (analogous to the singular
#     values in an SVD).}
#
#   \item{\code{f}}{A matrix whose columns contain the standardized
#     factors (i.e., with norm 1).}
#
flash_get_ldf = function(f, kset = NULL, drop_zero_factors = TRUE) {
  f = handle_f(f, allow_null = FALSE)
  kset = handle_kset(kset, f)

  ll = f$EL[, kset, drop=FALSE]
  ff = f$EF[, kset, drop=FALSE]
  d = sqrt(colSums(ll^2) * colSums(ff^2))

  ll = scale(ll, scale=sqrt(colSums(ll^2)), center=FALSE)
  ff = scale(ff, scale=sqrt(colSums(ff^2)), center=FALSE)

  if(drop_zero_factors) {
    ll = ll[, d!=0, drop=FALSE]
    ff = ff[, d!=0, drop=FALSE]
    d = d[d!=0, drop=FALSE]
  }

  return(list(d = d,
              l = ll,
              f = ff))
}


# @title Get number of factors in a flash object
#
# @description Returns the number of factors in a flash object.
#   Factors that have been zeroed out are not counted.
#
# @inheritParams flash_get_ldf
#
flash_get_nfactors = function(f) {
  f = handle_f(f, allow_null = FALSE)
  ldf = flash_get_ldf(f)
  return(length(ldf$d))
}


# @title Get PVE from a flash object
#
# @description Returns the factor contributions (proportion of
#   variance explained) for each factor/loading combination in flash
#   fit \code{f}. Because the factors are not required to be orthogonal,
#   this should be interpreted loosely: e.g., PVE could total more
#   than 1.
#
# @inheritParams flash_get_ldf
#
flash_get_pve = function(data, f, drop_zero_factors = TRUE) {
  f = handle_f(f, allow_null = FALSE)

  s = (flash_get_ldf(f, drop_zero_factors=drop_zero_factors)$d)^2
  if (is.null(f$tau)) {
    var_from_tau = 0
  } else if (is.matrix(f$tau)) {
    tau = f$tau[f$tau != 0]
    var_from_tau = sum(1/tau)
  } else { # tau is a scalar
    var_from_tau = sum(!data$missing) / f$tau
  }

  return(s/(sum(s) + var_from_tau))
}


# @title Get the residuals from a flash data and fit object,
#   excluding factor k.
flash_get_Rk = function(data, f, k, R = NULL) {
  if (is.null(R)) {
    return(
      data$Y - f$EL[, -k, drop = FALSE] %*% t(f$EF[, -k, drop = FALSE])
    )
  } else {
    return(R + outer(f$EL[, k], f$EF[, k]))
  }
}

# @title Get the residuals from data and a flash fit object.
flash_get_R = function(data, f) {
  # If f is null, return Y.
  if (is.null(f$EL)) {
    return(data$Y)
  } else {
    return(data$Y - f$EL %*% t(f$EF))
  }
}

# @title Get the residuals from data and a flash fit object, with
#   missing data as in original.
flash_get_R_withmissing = function(data, f) {
  # If f is null, return Y.
  if (is.null(f$EL)) {
    return(get_Yorig(data))
  } else {
    return(get_Yorig(data) - f$EL %*% t(f$EF))
  }
}

# @title Get the expected squared residuals from a flash data and fit
#   object.
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
    return(
      flash_get_Rk(data, f, k)^2 +
        f$EL2[, -k, drop = FALSE] %*% t(f$EF2[, -k, drop = FALSE]) -
        f$EL[, -k, drop = FALSE]^2 %*%
        t(f$EF[, -k, drop = FALSE]^2)
    )
  }
}

# @title is_tiny_fl
# @details Checks whether kth factor/loading combination is tiny.
is_tiny_fl = function(f, k, tol = 1e-08) {
  return(sum(f$EL[, k]^2) * sum(f$EF[, k]^2) < tol)
}

# @title Get number of factors in a fit object.
# @details Returns the number of factors in a flash fit.
flash_get_k = function(f) {
  k = ncol(f$EL)
  if (is.null(k)) {
    return(0)
  } else {
    return(k)
  }
}

flash_get_l = function(f) {
  return(f$EL)
}

flash_get_f = function(f) {
  return(f$EF)
}

flash_get_n = function(f) {
  return(nrow(f$EL))
}

flash_get_p = function(f) {
  return(nrow(f$EF))
}
