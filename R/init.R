# Contains functions related to initializing flash

#' @title  Initialize a flash object from the results of a factor analysis
#' @param LL the loadings, an n by K matrix
#' @param FF the factors, a p by K matrix
#' @return a flash object, with factors initialized using L and F
flash_init_LF = function(LL,FF){
  assertthat::assert_that(ncol(LL)==ncol(FF))
  f = list(EL = LL, EF = FF)
  f$EL2 = f$EL^2
  f$EF2 = f$EF^2
  f$gl = list()
  f$gf = list()
  f$tau = NULL
  return(f)
}

#' @title  Initialize a flash object from data using SVD
#' @param Y an n by p data matrix
#' @param K number of factors to use
#' @return a flash object, initialized using first K components of SVD
flash_init_svd = function(Y,K=1){
  Y.svd = svd(Y,nu=K,nv=K)
  LL = t(Y.svd$d[1:K,drop=FALSE] * t(Y.svd$u))
  return(flash_init_LF(LL,Y.svd$v))
}


#' @title  Initialize a flash object using random N(0,1) factors
#' @param Y and n by p data matrix (used to obtain n and p)
#' @param K a number of factors to initialize
#' @return a flash object, with loadings and factors initialized randomly iid N(0,1)
flash_init_random = function(Y,K=1){
  n = nrow(Y)
  p = ncol(Y)
  LL = matrix(rnorm(n*K),ncol=K)
  FF = matrix(rnorm(p*K),ncol=K)
  return(flash_init_LF(LL,FF))
}

#' @title  Initialize a flash object with K factors
#' @param Y an n times p data matrix
#' @param K a number of factors to initialize
#' @param method indicated how to initialize: can be "svd" or "random"
#' @return a flash object, with loadings and factors initialized
#' @export
flash_init = function(Y,K=1,method=c("svd","random")){
  method = match.arg(method)
  if(method=="svd")
    f=flash_init_svd(Y,K)
  else if (method=="random")
    f= flash_init_random(Y,K)
  else stop("illegal method")
  return(f)
}

#' @title combine two flash objects
#' @param f1 first flash object
#' @param f2 second flash object
#' @return a flash object whose factors are concatenations of f1 and f2
flash_combine = function(f1,f2){
  list(
    EL = cbind(f1$EL,f2$EL),
    EF = cbind(f1$EF,f2$EF),
    EL2 = cbind(f1$EL2,f2$EL2),
    EF2 = cbind(f1$EF2,f2$EF2),
    gl = c(f1$gl,f2$gl),
    gf = c(f1$gf,f2$gf),
    tau = NULL
  )
}

#' @title add a factor to f based on residuals
#' @param Y an n by p data matrix
#' @param f a flash object
#' @details Computes the current residuals from Y and f and adds K new factors based
#' on a simple initialization scheme applied to these residuals
flash_add_factor = function(Y,f,K=1,init_method=c("svd","random")){
  init_method = match.arg(init_method)
  R = get_R(Y,f)
  f2 = flash_init(R,K,init_method)
  flash_combine(f,f2)
}
