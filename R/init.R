# Contains functions related to initializing flash fit object

#' @title  Initialize a flash fit object with K factors
#' @param data a flash data object
#' @param K a number of factors to initialize
#' @param method indicated how to initialize: can be "svd" or "random"
#' @return a flash fit object, with loadings and factors initialized
#' @export
flash_init = function(data,K=1,method=c("svd","random")){
  method = match.arg(method)
  if(method=="svd")
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
  f$tau = NULL
  return(f)
}

#' @title  Initialize a flash fit object from data using SVD
#' @param data a flash data object
#' @param K number of factors to use
#' @return a flash fit object, initialized using first K components of SVD
flash_init_svd = function(data,K=1){
  Y.svd = svd(data$Y,nu=K,nv=K)
  LL = t(Y.svd$d[1:K,drop=FALSE] * t(Y.svd$u))
  return(flash_init_LF(LL,Y.svd$v))
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
  return(flash_init_LF(LL,FF))
}



#' @title combine two flash fit objects
#' @param f1 first flash fit object
#' @param f2 second flash fit object
#' @return a flash fit object whose factors are concatenations of f1 and f2
#' The precision (tau) of the combined fit is undefined (set to NULL)
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
    tau = NULL
  )
}

#' @title add a factor to f based on residuals
#' @param data a flash data object
#' @param f a flash fit object
#' @details Computes the current residuals from data and f and adds K new factors based
#' on a simple initialization scheme applied to these residuals
flash_add_factor = function(data,f,K=1,init_method=c("svd","random")){
  init_method = match.arg(init_method)
  R = get_R(data,f)
  f2 = flash_init(set_flash_data(R),K,init_method)
  flash_combine(f,f2)
}
