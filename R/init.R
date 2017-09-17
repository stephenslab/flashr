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
#' @return a flash fit object, with factors initialized using L and F
flash_init_LF = function(LL,FF){
  assertthat::assert_that(ncol(LL)==ncol(FF))
  f = list(EL = LL, EF = FF)
  f$EL2 = f$EL^2
  f$EF2 = f$EF^2
  f$gl = list()
  f$gf = list()
  f$ash_param_l = list()
  f$ash_param_f = list()
  f$comp_post_l = list() # component posteriors
  f$comp_post_f = list() # component posteriors
  f$KL_l = list()
  f$KL_f = list() #KL divergences for each l and f
  f$tau = NULL
  return(f)
}

#' @title  Initialize a flash fit object from data using softimpute (without penalty - lambda=0)
#' @param data a flash data object
#' @param K number of factors to use
#' @return a flash fit object, initialized using first K components of softimpute
flash_init_softImpute = function(data,K=1){
  # if missing values, need to use the original data to initialize
  if(data$anyNA){Y.si = softImpute::softImpute(data$Yorig, rank.max = K, type = "als",lambda = 0)}
  else{Y.si = softImpute::softImpute(data$Y, rank.max = K, type = "als",lambda = 0)}
  Y.si$u = as.matrix(Y.si$u)
  Y.si$v = as.matrix(Y.si$v)
  LL = t(Y.si$d * t(Y.si$u))
  f = flash_init_LF(LL,Y.si$v)
  f=flash_update_precision(data,f)
  return(f)
}

#' @title  Initialize a flash fit object from data using SVD
#' @param data a flash data object
#' @param K number of factors to use
#' @return a flash fit object, initialized using first K components of SVD
flash_init_svd = function(data,K=1){
  if(data$anyNA){stop("svd initialization can't be used with missing data")}
  Y.svd = svd(data$Y,nu=K,nv=K)
  LL = t(Y.svd$d[1:K,drop=FALSE] * t(Y.svd$u))
  f = flash_init_LF(LL,Y.svd$v)
  f=flash_update_precision(data,f)
  return(f)
}


#' @title  Initialize a flash fit object using random N(0,1) factors
#' @param data a flash data object
#' @param K a number of factors to initialize
#' @return a flash fit object, with loadings and factors initialized randomly iid N(0,1)
flash_init_random = function(data,K=1){
  n = nrow(data$Y)
  p = ncol(data$Y)
  LL = matrix(rnorm(n*K),ncol=K)
  FF = matrix(rnorm(p*K),ncol=K)
  f = flash_init_LF(LL,FF)
  f=flash_update_precision(data,f)
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
    gl = c(f1$gl,f2$gl),
    gf = c(f1$gf,f2$gf),
    ash_param_l = c(f1$ash_param_l, f2$ash_param_l),
    ash_param_f = c(f1$ash_param_f, f2$ash_param_f),
    comp_post_l = c(f1$comp_post_l,f2$comp_post_l),
    comp_post_f = c(f1$comp_post_f,f2$comp_post_f),
    KL_l = c(f1$KL_l,f2$KL_l),
    KL_f = c(f1$KL_f,f2$KL_f),
    tau = f1$tau
  )
}

#' @title add a factor to f based on residuals
#' @param data a flash data object
#' @param f a flash fit object
#' @details Computes the current residuals from data and f and adds K new factors based
#' on a simple initialization scheme applied to these residuals
flash_add_factor = function(data,f,K=1,init_method=c("softImpute","svd","random")){
  init_method = match.arg(init_method)
  R = get_R(data,f)
  f2 = flash_init(set_flash_data(R),K,init_method)
  f = flash_combine(f,f2)
  return(flash_update_precision(data,f))
}

#' @title zero out a factor from f
#' @param data a flash data object
#' @param f a flash fit object
#' @details the specified factor is made to be 0, effectively reducing the rank by 1
#' @export
flash_zero_out_factor = function(data,f,k=1){
  f$EL[,k] = 0
  f$EL2[,k] = 0
  f$EF[,k] = 0
  f$EF2[,k] = 0
  f$gl[[k]] = ashr::normalmix(1,0,0)
  f$gf[[k]] = ashr::normalmix(1,0,0)
  f$KL_l[[k]] = 0
  f$KL_f[[k]] = 0 #KL divergences for each l and f
  f=flash_update_precision(data,f)
  return(f)
}
