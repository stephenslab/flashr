
#' title rescale the lambda_l and lambda_f for identifiablity
#'
#' description rescale the lambda_l and lambda_f for identifiablity
#'
#' @return a list of sig2_l and sig2_f for the rescaled variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @param sig2_l variance of the kronecker product
#' @keywords internal
#'
rescale_sigmae2 = function(sig2_l,sig2_f){
  norm_l = sqrt(sum(sig2_l^2))
  norm_f = sqrt(sum(sig2_f^2))
  norm_total = norm_l * norm_f
  sig2_l = sig2_l / norm_l
  sig2_f = sig2_f / norm_f
  sig2_l = sig2_l * sqrt(norm_total)
  sig2_f = sig2_f * sqrt(norm_total)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}

#' title initial value for Bayes variance structure estimation for kronecker productor
#'
#' description initial value for Bayes variance structure estimation for kronecker productor
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @keywords internal
#'
inital_Bayes_var = function(sigmae2_v){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  # estimate the initial value of lambda_l and lambda_f
  sig2_l_pre = rowMeans(sigmae2_v, na.rm = TRUE)
  sig2_f_pre = colMeans( sigmae2_v / matrix(rep(sig2_l_pre,P), ncol = P) , na.rm = TRUE)
  sig2_pre_list = rescale_sigmae2(sig2_l_pre,sig2_f_pre)
  sig2_l_pre = sig2_pre_list$sig2_l
  sig2_f_pre = sig2_pre_list$sig2_f
  # start the iteration
  maxiter = 100
  inital_tol = 1e-3
  tau = 0
  epsilon = 1
  while(epsilon > inital_tol & tau <= maxiter){
    tau = tau + 1
    sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) , na.rm = TRUE)
    sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) , na.rm = TRUE)
    sig2_list = rescale_sigmae2(sig2_l,sig2_f)
    sig2_l = sig2_list$sig2_l
    sig2_f = sig2_list$sig2_f
    epsilon = sqrt(mean((sig2_f - sig2_f_pre)^2 )) + sqrt(mean((sig2_l - sig2_l_pre)^2))
    sig2_l_pre = sig2_l
    sig2_f_pre = sig2_f
  }
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}

#' title Bayes variance structure estimation for kronecker productor
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2 variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'
Bayes_var = function(sigmae2_v,sigmae2 = NULL){
  N = dim(sigmae2_v)[1]
  P = dim(sigmae2_v)[2]
  if( is.null(sigmae2) || !is.list(sigmae2) ){
    # we don't know the truth this is in the first iteration
    sigmae2 = inital_Bayes_var(sigmae2_v)
  }
  sig2_l_pre = sigmae2$sig2_l
  sig2_f_pre = sigmae2$sig2_f
  # this has already beed rescaled
  # here we use alpha_l = alpha_f = beta_l = beta_f = 0
  sig2_l = rowMeans( sigmae2_v / matrix(rep(sig2_f_pre,each = N),ncol = P) , na.rm = TRUE)
  sig2_f = colMeans( sigmae2_v / matrix(rep(sig2_l,P), ncol = P) , na.rm = TRUE)
  #rescaled the variance
  sig2_list = rescale_sigmae2(sig2_l,sig2_f)
  sig2_l = sig2_list$sig2_l
  sig2_f = sig2_list$sig2_f
  # sig2_out = matrix(rep(sig2_l,P),ncol = P) * matrix(rep(sig2_f,each = N),ncol = P)
  return(list(sig2_l = sig2_l,sig2_f = sig2_f))
}


#' title noisy variance structure estimation
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'

noisy_var = function(sigmae2_v,sigmae2_true){
  # a reaonable upper bound which control optimal algorithm
  upper_range = abs(mean(sigmae2_v - sigmae2_true , na.rm = TRUE)) + sd(sigmae2_v - sigmae2_true, na.rm = TRUE)
  # value of likelihood
  # I need the negative of the f_lik since optim find the minimum rather than the maximum
  f_lik <- function(x,sigmae2_v_l = sigmae2_v,sigmae2_true_l = sigmae2_true){
    (1/2)*mean( log(2*pi*(x+sigmae2_true_l)) + (1/(x+sigmae2_true_l))*sigmae2_v_l , na.rm = TRUE)
  }
  AA <- optim(mean(sigmae2_v - sigmae2_true, na.rm = TRUE), f_lik,lower = 0, upper = upper_range,  method = "Brent")
  e_sig = AA$par
  return(e_sig + sigmae2_true)
}


#' title noisy variance structure estimation with noisy variance on each column
#'
#' description prior and posterior part in objective function
#'
#' @return sig2_out estimated variance structure
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true variance structure
#' and it is a list of two vectors in this case.
#' @keywords internal
#'

noisy_var_column = function(sigmae2_v,sigmae2_true){
  # in this case, we need to think about the for each column
  # we can use noisy_var = function(sigmae2_v,sigmae2_true) with the input of each column
  P = dim(sigmae2_true)[2]
  sig2_out = sapply(seq(1,P),function(x){noisy_var(sigmae2_v[,x],sigmae2_true[,x])})
  return(sig2_out)
}


#' title module for estiamtion of the variance structure
#'
#' description estiamtion of the variance structure
#'
#' @return sigmae2 estimated variance structure
#' @param partype parameter type for the variance,
#' "constant" for constant variance,
#' "var_col" for nonconstant variance for column,
#' "known" for the kown variance,
#' "Bayes_var" for Bayes version of the nonconstant variance for row and column
#' "noisy" is noisy version s_ij + sigmae2
#' "noisy_col" noisy version for column s_ij + sigmae2_j
#' @param sigmae2_v  residual matrix
#' @param sigmae2_true  true variance structure
#' @param sigmae2_pre the previouse sigmae2 which is only used in Bayes_var
#' @keywords internal
#'
# sigma estimation function
tau_est = function(f,data){
  # to get the residual matrix 
  sigmae2_v = get_R2(data,f)
  sigmae2_v[data$missing] = NA
  N = get_n(f)
  P = get_p(f)
  if(data$noise_type == "var_col"){
    # this case, we condier the sigma_j column-wise variance
    sigmae2 = colMeans(sigmae2_v,na.rm = TRUE)
    f$tau = 1/as.matrix(rep(1,N) %*% t(sigmae2))
  } else if(data$noise_type == "noisy"){
    # in this case we consider the sigmae2tilde + sigmae2
    sigmae2 = noisy_var(sigmae2_v,data$sigmae2_tilde)
    f$tau = 1/sigmae2
  } else if(data$noise_type == "noisy_col"){
    # in this case we consider the sigmae2tilde + sigmae2_j
    sigmae2 = noisy_var_column(sigmae2_v,data$sigmae2_tilde)
    f$tau = 1/sigmae2
  }else if (data$noise_type == "Bayes_var"){
    # here sigmae2 is a list of two vectors and the variance is kronecker product of them
    f$sigmae2 = Bayes_var(sigmae2_v,f$sigmae2)
    sigmae2 = f$sigmae2$sig2_l %*% t(f$sigmae2$sig2_f)
    f$tau = 1/sigmae2
  }else if (data$noise_type == "known"){
    sigmae2 = data$sigmae2_tilde
    f$tau = 1/sigmae2
  } else {
    # this is for constant case
    sigmae2 = mean(sigmae2_v,na.rm = TRUE)
    f$tau = 1/matrix(sigmae2,ncol = P, nrow = N)
  }
  return(f)
}

