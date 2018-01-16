#' @title add factors or loadings to f
#' @details The precision parameter in f is updated after adding
#' @param data a flash data object
#' @param f a flash fit object
#' @param LL the loadings, an n by K matrix
#' @param FF the factors, a p by K matrix
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Useful for including a mean factor for example.
#' @param fixf a p by K matrix of TRUE/FALSE values; same as fixl but for factors FF.
#' @return a flash fit object, with additional factors initialized using LL and FF
#' @export
flash_add_lf = function(data,LL,FF,f_init=NULL,fixl=NULL,fixf=NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  f2 = flash_init_lf(LL,FF,fixl,fixf)
  f = flash_combine(f_init,f2)
  return(flash_update_precision(data,f))
}

#' @title add factors to a flash fit object based on data
#' @param data a flash data object
#' @param K number of factors to add
#' @param f_init an existing flash fit object to add to
#' @param init_fn the function to use to initialize new factors (typically some kind of svd-like function)
#' @details Computes the current residuals from data and f_init and adds K new factors based
#' on init_fn applied to these residuals. (If f_init is NULL then the residuals are the data)
#' @export
flash_add_factors_from_data = function(data,K,f_init=NULL,init_fn="udv_si"){
  if(is.null(f_init)){f_init = flash_init_null()}
  R = get_R_withmissing(data,f_init)
  f2 = flash_init_fn(flash_set_data(R),init_fn,K)
  f = flash_combine(f_init,f2)
  return(flash_update_precision(data,f))
}



#' @title  Add a set of fixed loadings to a flash fit object
#' @param data a flash data object
#' @param LL the loadings, an n by K matrix. Missing values will be initialized by the mean of the relevant column (but will generally be
#' re-estimated when refitting the model).
#' @param f_init a flash fit object to which loadings are to be added (if NULL then a new fit object is created)
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Default is to fix all non-missing values, so missing values will be updated when f is updated.
#' @return a flash fit object, with loadings initialized from LL, and corresponding factors initialized to 0.
#' @export
flash_add_fixed_l = function(data, LL, f_init=NULL, fixl = NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixl)){fixl = !is.na(LL)}
  LL = fill_missing_with_column_mean(LL)
  FF = matrix(0,nrow=ncol(data$Y),ncol=ncol(LL))

  f_new = flash_init_lf(LL,FF,fixl=fixl)
  f = flash_combine(f_init,f_new)

  # maybe in future we want to give a fit option? But then would
  # need to pass in var_type? possibly not needed.
  # if(fit){
  #   k1 = get_k(f_init)
  #   k2 = get_k(f)
  #   f = flash_backfit(data,f,kset=((k1+1):k2),var_type=xx)
  # }
  return(f)
}

#' @title  Add a set of fixed factors to a flash fit object
#' @param data a flash data object
#' @param FF the factors, a p by K matrix. Missing values will be initialized by the mean of the relevant column (but will generally be
#' re-estimated when refitting the model).
#' @param f_init a flash fit object to which factors are to be added (if NULL then a new fit object is created)
#' @param fixf a p by K matrix of TRUE/FALSE values indicating which elements of FF should be considered fixed and not changed during updates.
#' Default is to fix all non-missing values, so missing values will be updated when f is updated.
#' @return a flash fit object, with factors initialized from FF, and corresponding loadings initialized to 0.
#' @export
flash_add_fixed_f = function(data, FF, f_init=NULL, fixf = NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixf)){fixf = !is.na(FF)}
  FF = fill_missing_with_column_mean(FF)
  LL = matrix(0,nrow=nrow(data$Y),ncol=ncol(FF))

  f_new = flash_init_lf(LL,FF,fixf=fixf)
  f = flash_combine(f_init,f_new)

  # maybe in future we want to give a fit option? But then would
  # need to pass in var_type? possibly not needed.
  # if(fit){
  #   k1 = get_k(f_init)
  #   k2 = get_k(f)
  #   f = flash_backfit(data,f,kset=((k1+1):k2),var_type=xx)
  # }
  return(f)
}


NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))

fill_missing_with_column_mean = function(X){
  apply(X, 2, NA2mean)
}


#' @title Add a set of sparse loadings to a flash fit object
#' @param data an n by p matrix or a flash data object.
#' @param nonnull_LL an n by K matrix of TRUE/FALSE values. Each column corresponds to
#' a single sparse loading. Within each column, a value of FALSE indicates that the
#' respective element should be initialized to zero. Other elements are treated as
#' potentially non-null, and are initialized using init_fn.
#' @param f_init a flash fit object to which factors are to be added (if NULL then a
#' new fit object is created).
#' @param init_fn the function to use to initialize new factors (typically some kind
#' of svd-like function). Since only a subset of the loading elements need to be
#' initialized, init_fn is applied to a subsetted flash object.
#' @param fix_nulls indicates whether nulls should be regarded as fixed during
#' subsequent updates. Default is TRUE.
#' @return a flash fit object, with additional sparse factors initialized using nonnull_LL
#' @examples
#' sparse_ll = c(rep(1, 2), rep(0, 8))
#' Y = matrix(rnorm(200),nrow=10,ncol=20) + outer(sparse_ll, rnorm(20))
#' f = flash_add_sparse_l(Y, matrix(c(rep(TRUE, 2), rep(FALSE, 8)), ncol=1))
#' f = flash_backfit(Y, f, var_type="constant")
#' @export
flash_add_sparse_l = function(data, nonnull_LL, f_init=NULL, init_fn="udv_si", fix_nulls=TRUE) {
  if(is.matrix(data)){data = flash_set_data(data)}
  f = f_init
  if(is.null(f)){f = flash_init_null()}

  assertthat::assert_that(sum(colSums(nonnull_LL == TRUE) < 2) == 0,
                          msg = "Each loadings vector must contain at least two non-null elements")

  n = nrow(nonnull_LL)
  K = ncol(nonnull_LL)
  fixl = NULL

  for (i in 1:K) {
    subset = nonnull_LL[,i]
    subf = flash_subset_l(f, subset)
    subdata = flash_subset_l_data(data, subset)

    subf = flash_add_factors_from_data(subdata, 1, subf, init_fn)
    next_ll = matrix(0, nrow=n, ncol=1)
    next_ll[subset,] = subf$EL[,ncol(subf$EL)]
    next_ff = matrix(subf$EF[,ncol(subf$EF)], ncol=1)

    if (fix_nulls == TRUE) {
      fixl = matrix(!subset, ncol=1)
    }
    new_f = flash_init_lf(next_ll, next_ff, fixl=fixl)
    f = flash_combine(f, new_f)
  }
  return(f)
}


#' @title Add a set of sparse factors to a flash fit object
#' @param data an n by p matrix or a flash data object.
#' @param nonnull_FF a p by K matrix of TRUE/FALSE values. Each column corresponds to
#' a single sparse factor. Within each column, a value of FALSE indicates that the
#' respective element should be initialized to zero. Other elements are treated as
#' potentially non-null, and are initialized using init_fn.
#' @param f_init a flash fit object to which factors are to be added (if NULL then a
#' new fit object is created).
#' @param init_fn the function to use to initialize new factors (typically some kind
#' of svd-like function). Since only a subset of the factor elements need to be
#' initialized, init_fn is applied to a subsetted flash object.
#' @param fix_nulls indicates whether nulls should be regarded as fixed during
#' subsequent updates. Default is TRUE.
#' @return a flash fit object, with additional sparse factors initialized using nonnull_FF
#' @examples
#' sparse_ff = c(rep(5, 2), rep(0, 18))
#' Y = matrix(rnorm(100),nrow=5,ncol=20) + outer(rnorm(5), sparse_ff)
#' f = flash_add_sparse_f(Y, matrix(c(rep(TRUE, 2), rep(FALSE, 18)), ncol=1))
#' f = flash_backfit(Y, f, var_type="constant")
#' @export
flash_add_sparse_f = function(data, nonnull_FF, f_init=NULL, init_fn="udv_si", fix_nulls=TRUE) {
  if(is.matrix(data)){data = flash_set_data(data)}
  f = f_init
  if(is.null(f)){f = flash_init_null()}

  assertthat::assert_that(sum(colSums(nonnull_FF == TRUE) < 2) == 0,
                          msg = "Each loadings vector must contain at least two non-null elements")

  p = nrow(nonnull_FF)
  K = ncol(nonnull_FF)
  fixf = NULL

  for (i in 1:K) {
    subset = nonnull_FF[,i]
    subf = flash_subset_f(f, subset)
    subdata = flash_subset_f_data(data, subset)

    subf = flash_add_factors_from_data(subdata, 1, subf, init_fn)
    next_ff = matrix(0, nrow=p, ncol=1)
    next_ff[subset,] = subf$EF[,ncol(subf$EF)]
    next_ll = matrix(subf$EL[,ncol(subf$EL)], ncol=1)

    if (fix_nulls == TRUE) {
      fixf = matrix(!subset, ncol=1)
    }
    new_f = flash_init_lf(next_ll, next_ff, fixf=fixf)
    f = flash_combine(f, new_f)
  }
  return(f)
}


#' @title Subset a flash object with respect to its loadings
#' @param f a flash fit object
#' @param subset the subset of loading elements to be retained
#' @return a subsetted flash fit object
flash_subset_l = function(f, subset){
  subf = f
  subf$EL = subf$EL[subset,,drop=F]
  subf$EL2 = subf$EL2[subset,,drop=F]
  subf$fixl = subf$fixl[subset,,drop=F]
  subf$tau = subf$tau[subset,,drop=F]
  subf$KL_l = NULL
  subf$KL_f = NULL
  return(subf)
}


#' @title Subset a flash object with respect to its factors
#' @param f a flash fit object
#' @param subset the subset of factor elements to be retained
#' @return a subsetted flash fit object
flash_subset_f = function(f, subset){
  subf = f
  subf$EF = subf$EF[subset,,drop=F]
  subf$EF2 = subf$EF2[subset,,drop=F]
  subf$fixf = subf$fixf[subset,,drop=F]
  subf$tau = subf$tau[,subset,drop=F]
  subf$KL_l = NULL
  subf$KL_f = NULL
  return(subf)
}


#' @title Subset a flash data object with respect to the rows of Y
#' @param f a flash fit object
#' @param subset the subset of rows to be retained
#' @return a subsetted flash data object
flash_subset_l_data = function(data, subset){
  subdata = data
  subdata$Yorig = subdata$Yorig[subset,,drop=F]
  subdata$anyNA = anyNA(subdata$Yorig)
  subdata$missing = subdata$missing[subset,,drop=F]
  subdata$Y = subdata$Y[subset,,drop=F]
  return(subdata)
}


#' @title Subset a flash data object with respect to the columns of Y
#' @param f a flash fit object
#' @param subset the subset of columns to be retained
#' @return a subsetted flash data object
flash_subset_f_data = function(data, subset){
  subdata = data
  subdata$Yorig = subdata$Yorig[,subset,drop=F]
  subdata$anyNA = anyNA(subdata$Yorig)
  subdata$missing = subdata$missing[,subset,drop=F]
  subdata$Y = subdata$Y[,subset,drop=F]
  return(subdata)
}


#' @title  Initialize a flash fit object by applying a function to data
#' @param data a flash data object
#' @param init_fn an initialization function, which takes as input an (n by p matrix, or flash data object)
#' and K, a number of factors, and and outputs a list with elements (u,d,v)
#' @return a flash fit object
flash_init_fn = function(data,init_fn,K=1){
  s = do.call(init_fn,list(get_Yorig(data),K))
  f = flash_init_udv(s,K)
  f = flash_update_precision(data,f)
  return(f)
}
