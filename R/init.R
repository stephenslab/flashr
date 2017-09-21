# Contains functions related to initializing flash fit object

#' @title  Initialize a flash fit object with K factors
#' @param data a flash data object
#' @param K a number of factors to initialize
#' @param method indicated how to initialize: can be "svd" or "random"
#' @return a flash fit object, with loadings and factors initialized
#' @export
flash_init = function(data,K=1,method=c("softImpute","svd","random")){
  method = match.arg(method)
  if(method=="softImpute")
    f=flash_init_softImpute(data,K)
  else if(method=="svd")
    f=flash_init_svd(data,K)
  else if (method=="random")
    f= flash_init_random(data,K)
  else stop("illegal method")
  return(f)
}


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

  f = list(EL = LL, EF = FF, fixl = fixl, fixf=fixf)

  f$EL2 = f$EL^2
  f$EF2 = f$EF^2
  f$gl = list()
  f$gf = list()
  f$ash_param_l = list()
  f$ash_param_f = list()
  f$KL_l = as.list(rep(0,get_k(f)))
  f$KL_f = as.list(rep(0,get_k(f))) #KL divergences for each l and f
  f$tau = NULL
  return(f)
}

#' @title  Initialize a flash fit object from a list with elements (u,d,v)
#' @param data a flash data object
#' @param r1_fn an initialization function, which takes as input an (n by p matrix, or flash data object)
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
#' The precision (tau) of the combined fit is inherited from f1
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
    tau = f1$tau
  )
}

#' @title add a factor to f based on residuals
#' @param data a flash data object
#' @param f a flash fit object
#' @param init_fn the function to use to initialize new factors
#' @param K number of factors
#' @details Computes the current residuals from data and f and adds K new factors based
#' on init_fn applied to these residuals
flash_add_factor = function(data,f,init_fn,K=1){
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
