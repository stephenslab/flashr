# Contains functions related to initializing flash fit object

#' @title  Initialize a flash fit object from the results of a factor analysis
#' @param LL the loadings, an n by K matrix
#' @param FF the factors, a p by K matrix
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Useful for including a mean factor for example.
#' @param fixf a p by K matrix of TRUE/FALSE values; same as fixl but for factors FF.
#' @return a flash fit object, with factors initialized using L and F
flash_init_LF = function(LL,FF, fixl = NULL, fixf = NULL){
  assertthat::assert_that(ncol(LL)==ncol(FF))
  if(is.null(fixl)){fixl = matrix(FALSE,ncol=ncol(LL),nrow=nrow(LL))}
  if(is.null(fixf)){fixf = matrix(FALSE,ncol=ncol(FF),nrow=nrow(FF))}

  f = list(EL = LL, EF = FF, EL2 = LL^2, EF2= FF^2, fixl = fixl, fixf=fixf)

  f$gl = list()
  f$gf = list()
  f$ash_param_l = list()
  f$ash_param_f = list()
  f$KL_l = as.list(rep(0,get_k(f)))
  f$KL_f = as.list(rep(0,get_k(f))) #KL divergences for each l and f
  f$tau = NULL
  return(f)
}

#' @title  Add a set of fixed loadings to a flash fit object
#' @param data a flash data object
#' @param LL the loadings, an n by K matrix. Missing values will be initialized by the mean of the relevant column (but will generally be
#' re-estimated when refitting the model).
#' @param f_init a flash fit object to which loadings are to be added (if NULL then a new fit object is created)
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Default is to fix all non-missing values, so missing values will be updated when f is updated.
#' @return a flash fit object, with loadings initialized from LL, and corresponding factors initialized to 0.
flash_add_fixed_L = function(data, LL, f_init=NULL, fixl = NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixl)){fixl = !is.na(LL)}
  LL = fill_missing_with_column_mean(LL)
  FF = matrix(0,nrow=ncol(data$Y),ncol=ncol(LL))

  f_new = flash_init_LF(LL,FF,fixl=fixl)
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
flash_add_fixed_F = function(data, FF, f_init=NULL, fixf = NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixf)){fixf = !is.na(FF)}
  FF = fill_missing_with_column_mean(FF)
  LL = matrix(0,nrow=nrow(data$Y),ncol=ncol(FF))

  f_new = flash_init_LF(LL,FF,fixf=fixf)
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

#' @title  Initialize an empty flash fit object
#' @return an empty flash fit object
flash_init_null = function(){
  f = list(EL = NULL, EF = NULL, fixl = NULL,
           fixf=NULL, EL2=NULL, EF2=NULL,
           gl =NULL, gf = NULL, ash_param_l = NULL,
           ash_param_f = NULL, KL_l = NULL, KL_f = NULL, tau=NULL)
  return(f)
}


#' @title  Initialize a flash fit object by applying a function to data
#' @param data a flash data object
#' @param init_fn an initialization function, which takes as input an (n by p matrix, or flash data object)
#' and K, a number of factors, and and outputs a list with elements (u,d,v)
#' @return a flash fit object
#' @export
flash_init_fn = function(data,init_fn,K=1){
  s = do.call(init_fn,list(get_Yorig(data),K))
  f = flash_init_udv(s,K)
  f = flash_update_precision(data,f)
  return(f)
}

#' @title udv_si
#' @details provides a simple wrapper to \code{softImpute} to provide a rank 1 initialization
#' @param Y an n by p matrix
#' @param K number of factors to use
#' @return a list with components (u,d,v)
#' @export
udv_si = function(Y,K=1){
  softImpute::softImpute(Y, rank.max = K, type = "svd",lambda = 0)
}

#' @title udv_svd
#' @details provides a simple wrapper to svd
#' @param Y an n by p matrix
#' @param K number of factors to use
#' @return a list with components (u,d,v)
#' @export
udv_svd = function(Y,K=1){
  svd(Y,K,K)
}

#' @title udv_random
#' @details provides a random initialization of factors
#' @param Y an n by p matrix
#' @param K number of factors
#' @return a list with components (u,d,v), with elements of u and v iid N(0,1)
#' @export
udv_random = function(Y,K=1){
  n = nrow(Y)
  p = ncol(Y)
  list(u=matrix(rnorm(n*K),ncol=K),d=1,v=matrix(rnorm(p*K),ncol=K))
}


#' @title  Initialize a flash fit object from a list with elements (u,d,v)
#' @param s list with elements (u,v,d)
#' @param K the number of factors to use (factors 1:K are used)
#' @return a flash fit object ready for optimization
flash_init_udv = function(s,K=1){
  s$u = as.matrix(s$u)
  s$v = as.matrix(s$v)
  if(ncol(s$u)>K){s$u = s$u[,1:K,drop=FALSE]} # deals with case these are vectors (K=1)
  if(ncol(s$v)>K){s$v = s$v[,1:K,drop=FALSE]}
  if(length(s$d)>K){s$d = s$d[1:K]}

  f = flash_init_LF(t(s$d * t(s$u)) , s$v)
  return(f)
}



#' @title combine two flash fit objects
#' @param f1 first flash fit object
#' @param f2 second flash fit object
#' @return a flash fit object whose factors are concatenations of f1 and f2
#' The precision (tau) of the combined fit is inherited from f2
flash_combine = function(f1,f2){
  list(
    EL = cbind(f1$EL,f2$EL),
    EF = cbind(f1$EF,f2$EF),
    EL2 = cbind(f1$EL2,f2$EL2),
    EF2 = cbind(f1$EF2,f2$EF2),
    fixl = cbind(f1$fixl, f2$fixl),
    fixf = cbind(f1$fixf, f2$fixf),
    gl = c(f1$gl,f2$gl),
    gf = c(f1$gf,f2$gf),
    ash_param_l = c(f1$ash_param_l, f2$ash_param_l),
    ash_param_f = c(f1$ash_param_f, f2$ash_param_f),
    KL_l = c(f1$KL_l,f2$KL_l),
    KL_f = c(f1$KL_f,f2$KL_f),
    tau = f2$tau
  )
}



#' @title add factors to f based on residuals
#' @param data a flash data object
#' @param f a flash fit object
#' @param init_fn the function to use to initialize new factors
#' @param K number of factors
#' @details Computes the current residuals from data and f and adds K new factors based
#' on init_fn applied to these residuals
flash_add_factors_from_residuals = function(data,f,init_fn,K=1){
  R = get_R(data,f)
  f2 = flash_init_fn(set_flash_data(R),init_fn,K)
  f = flash_combine(f,f2)
  return(flash_update_precision(data,f))
}

#' @title zero out a factor from f
#' @param data a flash data object
#' @param f a flash fit object
#' @param k index of factor/loading to zero out
#' @details The factor and loadings of the kth factor of f are made to be zero.
#' (Except for elements of the factor/loading that are designated to be fixed)
#' This effectively reduces the rank by 1, although the zero factor/loading is kept in f
#' so the number and indexing of factor/loading matrices in f remains the same
#' @export
flash_zero_out_factor = function(data,f,k=1){
  f$EL[!f$fixl[,k],k] = 0
  f$EL2[!f$fixl[,k],k] = 0
  f$EF[!f$fixf[,k],k] = 0
  f$EF2[!f$fixf[,k],k] = 0
  f$gl[[k]] = ashr::normalmix(1,0,0)
  f$gf[[k]] = ashr::normalmix(1,0,0)
  f$KL_l[[k]] = 0
  f$KL_f[[k]] = 0 #KL divergences for each l and f
  f=flash_update_precision(data,f)
  return(f)
}

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
flash_add_LF = function(data,f,LL,FF,fixl=NULL,fixf=NULL){
  f2 = flash_init_LF(LL,FF,fixl,fixf)
  f = flash_combine(f,f2)
  return(flash_update_precision(data,f))
}
