#' @title Use a flash fit to fill in missing entries
#' @details Fills in missing entries of Y by using the relevant entries of the estimated LDF' from the flash fit.
#' @param Y an n by p matrix with missing entries used to fit flash
#' @param f the flash fit object obtained from running flash on Y
#' @return a matrix with non-missing entries the same as $Y$, and missing entries imputed from the flash fit
#' @export
flash_fill = function(Y, f){
  missing = is.na(Y)
  if(!is.matrix(Y)){stop("for flash_fill Y must be a matrix")}
  if(dim(Y)[1]!=get_n(f)){stop("dimensions of Y must match flash fit")}
  if(dim(Y)[2]!=get_p(f)){stop("dimensions of Y must match flash fit")}
  Y[is.na(Y)] = flash_get_lf(f)[is.na(Y)]
  return(Y)
}




#' @title transpose a flash fit object
#' @param f the flash fit object
#' @return A new flash fit object, with the factors and loadings of the original flash
#' fit object interchanged.
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
    tmp[c(which(tmp == "ebnm_param_l"), which(tmp == "ebnm_param_f"))] = c("ebnm_param_f", "ebnm_param_l")
    tmp[c(which(tmp == "penloglik_l"), which(tmp == "penloglik_f"))] = c("penloglik_f", "penloglik_l")
    names(f) = tmp
    if (is.matrix(f$tau)) {
        f$tau = t(f$tau)
    }
    return(f)
}

#' @title transpose a flash data object
#' @param f the flash data object
#' @return A new flash data object, with the matrices of the original flash data
#' object transposed.
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
    return(data)
}

#' @title combine two flash fit objects
#' @param f1 first flash fit object
#' @param f2 second flash fit object
#' @return A flash fit object whose factors are concatenations of f1 and f2.
#' The precision (tau) of the combined fit is inherited from f2.
flash_combine = function(f1, f2) {
    list(EL = cbind(f1$EL, f2$EL), EF = cbind(f1$EF, f2$EF), EL2 = cbind(f1$EL2, f2$EL2), EF2 = cbind(f1$EF2, f2$EF2),
        fixl = cbind(f1$fixl, f2$fixl), fixf = cbind(f1$fixf, f2$fixf), gl = c(f1$gl, f2$gl), gf = c(f1$gf, f2$gf),
        ebnm_param_l = c(f1$ebnm_param_l, f2$ebnm_param_l), ebnm_param_f = c(f1$ebnm_param_f, f2$ebnm_param_f), KL_l = c(f1$KL_l,
            f2$KL_l), KL_f = c(f1$KL_f, f2$KL_f), penloglik_l = c(f1$penloglik_l, f2$penloglik_l), penloglik_f = c(f1$penloglik_f,
            f2$penloglik_f), tau = f2$tau)
}

#' @title Subset a flash object with respect to its loadings
#' @param f a flash fit object
#' @param subset the subset of loading elements to be retained
#' @return a subsetted flash fit object
flash_subset_l = function(f, subset) {
    subf = f
    subf$EL = subf$EL[subset, , drop = F]
    subf$EL2 = subf$EL2[subset, , drop = F]
    subf$fixl = subf$fixl[subset, , drop = F]
    subf$tau = subf$tau[subset, , drop = F]
    subf$KL_l = NULL
    subf$KL_f = NULL
    return(subf)
}

#' @title Subset a flash object with respect to its factors
#' @param f a flash fit object
#' @param subset the subset of factor elements to be retained
#' @return a subsetted flash fit object
flash_subset_f = function(f, subset) {
    subf = f
    subf$EF = subf$EF[subset, , drop = F]
    subf$EF2 = subf$EF2[subset, , drop = F]
    subf$fixf = subf$fixf[subset, , drop = F]
    subf$tau = subf$tau[, subset, drop = F]
    subf$KL_l = NULL
    subf$KL_f = NULL
    return(subf)
}

#' @title Subset a flash data object
#' @param f a flash fit object
#' @param row_subset the subset of rows to be retained
#' @param col_subset the subset of columns to be retained
#' @return a subsetted flash data object
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
    return(subdata)
}

#' @title zero out a factor from f
#' @param data a flash data object
#' @param f a flash fit object
#' @param k index of factor/loading to zero out
#' @details The factor and loadings of the kth factor of f are made to be zero
#' (except for elements of the factor/loading that are designated to be fixed).
#' This effectively reduces the rank by 1, although the zero factor/loading is kept in f
#' so the number and indexing of factor/loading matrices in f remains the same.
#' @export
flash_zero_out_factor = function(data, f, k = 1) {
    f$EL[!f$fixl[, k], k] = 0
    f$EL2[!f$fixl[, k], k] = 0
    f$EF[!f$fixf[, k], k] = 0
    f$EF2[!f$fixf[, k], k] = 0
    f$gl[[k]] = ashr::normalmix(1, 0, 0)
    f$gf[[k]] = ashr::normalmix(1, 0, 0)
    f$KL_l[[k]] = 0
    f$KL_f[[k]] = 0  #KL divergences for each l and f
    f = flash_update_precision(data, f)
    return(f)
}
