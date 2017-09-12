# Contains functions related to initializing flash

#' @title  Create a flash object for a given data-set
#' @param Y an n by p data matrix
#' @return a flash object, a list that can be passed to other functions
#' @export
flash_set_data = function(Y,K=1){
  n = nrow(Y)
  p = ncol(Y)
  EL = EL2 = matrix(1,nrow=n,ncol=K)
  EF = EF2 = matrix(1,nrow=p,ncol=K)
  gl = list()
  gf = list()
  tau = matrix(1,nrow=n,ncol=p)
  list(Y = Y, EL = EL, EL2 = EL2,
       EF = EF, EF2 = EF2, gl = gl, gf = gf, tau = tau)
}

#' @title  Initialize a flash object using SVD
#' @param f a flash object
#' @param K a number of factors to initialize
#' @return a flash object, with factors initialized using first k eigenvectors and values using SVD
#' @export
flash_init_svd = function(f,K=1){
  Y.svd = svd(f$Y,nu=K,nv=K)
  f$EL = t( Y.svd$d[1:K,drop=FALSE] * t(Y.svd$u) )
  f$EF = Y.svd$v
  f$EL2 = f$EL^2
  f$EF2 = f$EF^2
  f = flash_update_precision(f)
  f$gl = list()
  f$gf = list()
  return(f)
}

#' @title  Initialize a flash object using random N(0,1) factors
#' @param f a flash object
#' @param K a number of factors to initialize
#' @return a flash object, with loadings and factors initialized randomly iid N(0,1)
#' @export
flash_init_random = function(f,K=1){
  n = get_n(f)
  p = get_p(f)
  f$EL = matrix(rnorm(n*K),ncol=K)
  f$EF = matrix(rnorm(p*K),ncol=K)
  f$EL2 = f$EL^2
  f$EF2 = f$EF^2
  f = flash_update_precision(f)
  f$gl = list()
  f$gf = list()
  return(f)
}

#' @title  Initialize a flash object with K factors
#' @param f a flash object
#' @param K a number of factors to initialize
#' @param method indicated how to initialize: can be "svd" or "random"
#' @return a flash object, with loadings and factors initialized
#' @export
flash_init = function(f,K=1,method=c("svd","random")){
  method = match.arg(method)
  if(method=="svd")
    return(flash_init_svd(f,K))
  else if (method=="random")
    return(flash_init_random(f,K))
}
