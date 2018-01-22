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
