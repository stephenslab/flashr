#' @title Use a flash fit to fill in missing entries
#'
#' @description Fills in missing entries of Y by using the relevant
#'   entries of the estimated LDF' from the flash fit.
#'
#' @inheritParams flash
#'
#' @param f A fitted flash object obtained from running \code{flash} on
#'   \code{data}.
#'
#' @return A matrix with non-missing entries the same as Y, and
#'   missing entries imputed from the flash fit.
#'
#' @export
#'
flash_fill = function(data, f) {
  f = handle_f(f)
  data = handle_data(data, f, output = "matrix")

  data[is.na(data)] = flash_get_fitted_values(f)[is.na(data)]

  return(data)
}


#' @title Zero out factor from flash object
#'
#' @description The factor and loadings of the kth factor of \code{f}
#'   are made to be zero (except for elements of the factor/loading that
#'   are designated to be fixed). This effectively reduces the rank by 1,
#'   but the zero factor/loading is retained in \code{f} so that
#'   the number and indexing of factor/loading matrices in \code{f}
#'   remain the same.
#'
#' @param f A flash fit object.
#'
#' @param k The index of the factor/loading pair to zero out.
#'
#' @return The flash object with the kth factor/loading zeroed out.
#'
#' @export
#'
flash_zero_out_factor = function(f, k) {
  f = handle_f(f, allow_null = FALSE)
  k = handle_k(k, f)

  f$EL[!f$fixl[, k], k] = 0
  f$EL2[!f$fixl[, k], k] = 0
  f$EF[!f$fixf[, k], k] = 0
  f$EF2[!f$fixf[, k], k] = 0
  f$gl[[k]] = list(NULL)
  f$gf[[k]] = list(NULL)
  f$KL_l[[k]] = 0
  f$KL_f[[k]] = 0

  return(f)
}


# @title Transpose a flash fit object
#
# @param f A flash fit object.
#
# @return A new flash fit object, with the factors and loadings of the
#   original flash fit object interchanged.
#
flash_transpose = function(f) {
  if (is.null(f)) {
    return(NULL)
  }

  tmp = names(f)
  tmp[c(which(tmp == "EL"), which(tmp == "EF"))] = c("EF", "EL")
  tmp[c(which(tmp == "EL2"), which(tmp == "EF2"))] = c("EF2", "EL2")
  tmp[c(which(tmp == "fixl"), which(tmp == "fixf"))] = c("fixf", "fixl")
  tmp[c(which(tmp == "gl"), which(tmp == "gf"))] = c("gf", "gl")
  tmp[c(which(tmp == "KL_l"), which(tmp == "KL_f"))] = c("KL_f", "KL_l")
  tmp[c(which(tmp == "ebnm_fn_l"),
        which(tmp == "ebnm_fn_f"))] = c("ebnm_fn_f", "ebnm_fn_l")
  tmp[c(which(tmp == "ebnm_param_l"),
        which(tmp == "ebnm_param_f"))] = c("ebnm_param_f", "ebnm_param_l")
  names(f) = tmp

  if (is.matrix(f$tau)) {
    f$tau = t(f$tau)
  }

  return(f)
}


# @title Transpose a flash data object
#
# @param f The flash data object.
#
# @return A new flash data object, with the matrices of the original
#   flash data object transposed.
#
flash_transpose_data = function(data) {
  if (is.matrix(data$Yorig)) {
    data$Yorig = t(data$Yorig)
  }
  if (is.matrix(data$missing)) {
    data$missing = t(data$missing)
  }
  if (is.matrix(data$Y)) {
    data$Y = t(data$Y)
  }
  if (is.matrix(data$S)) {
    data$S = t(data$S)
  }

  return(data)
}


# @title Combine two flash fit objects
#
# @param f1 The first flash fit object.
#
# @param f2 The second flash fit object.
#
# @return A flash fit object whose factors are concatenations of f1
#   and f2. If both precision matrices (tau) are nonnull, then the
#   combined fit inherits tau from f2.
#
flash_combine = function(f1, f2) {
  if (is.null(f2$tau)) {
    tau = f1$tau
  } else {
    tau = f2$tau
  }

  f = list(EL = cbind(f1$EL, f2$EL),
           EF = cbind(f1$EF, f2$EF),
           EL2 = cbind(f1$EL2, f2$EL2),
           EF2 = cbind(f1$EF2, f2$EF2),
           fixl = cbind(f1$fixl, f2$fixl),
           fixf = cbind(f1$fixf, f2$fixf),
           gl = c(f1$gl, f2$gl),
           gf = c(f1$gf, f2$gf),
           ebnm_fn_l = c(f1$ebnm_fn_l, f2$ebnm_fn_l),
           ebnm_fn_f = c(f1$ebnm_fn_f, f2$ebnm_fn_f),
           ebnm_param_l = c(f1$ebnm_param_l, f2$ebnm_param_l),
           ebnm_param_f = c(f1$ebnm_param_f, f2$ebnm_param_f),
           KL_l = c(f1$KL_l, f2$KL_l),
           KL_f = c(f1$KL_f, f2$KL_f),
           tau = tau)
  class(f) = "flash"

  return(f)
}


# @title Subset a flash object with respect to its loadings
#
# @param f A flash fit object.
#
# @param subset The subset of loading elements to be retained.
#
# @return A subsetted flash fit object.
#
flash_subset_l = function(f, subset) {
  subf = f
  subf$EL = subf$EL[subset, , drop = F]
  subf$EL2 = subf$EL2[subset, , drop = F]
  subf$fixl = subf$fixl[subset, , drop = F]
  subf$tau = subf$tau[subset, , drop = F]
  subf$KL_l = list(NULL)
  subf$KL_f = list(NULL)

  return(subf)
}


# @title Subset a flash object with respect to its factors
#
# @param f A flash fit object.
#
# @param Subset the subset of factor elements to be retained.
#
# @return A subsetted flash fit object.
#
flash_subset_f = function(f, subset) {
  subf = f
  subf$EF = subf$EF[subset, , drop = F]
  subf$EF2 = subf$EF2[subset, , drop = F]
  subf$fixf = subf$fixf[subset, , drop = F]
  subf$tau = subf$tau[, subset, drop = F]
  subf$KL_l = list(NULL)
  subf$KL_f = list(NULL)

  return(subf)
}


# @title Subset a flash data object
#
# @param f A flash fit object.
#
# @param row_subset The subset of rows to be retained.
#
# @param col_subset The subset of columns to be retained.
#
# @return A subsetted flash data object.
#
flash_subset_data = function(data, row_subset = NULL, col_subset = NULL) {
  if (is.null(row_subset)) {
    row_subset = 1:nrow(data$Y)
  }
  if (is.null(col_subset)) {
    col_subset = 1:ncol(data$Y)
  }

  subdata = data
  subdata$Yorig = subdata$Yorig[row_subset, col_subset, drop = F]
  subdata$anyNA = anyNA(subdata$Yorig)
  subdata$missing = subdata$missing[row_subset, col_subset, drop = F]
  subdata$Y = subdata$Y[row_subset, col_subset, drop = F]
  class(subdata) = "flash_data"

  return(subdata)
}
