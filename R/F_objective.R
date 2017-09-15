
#currently just does it for K=1 factor
get_F = function(data,f){
  return(e_log_lik(data,f) - sum(KL_post_prior_f(data,f,1)) - sum(KL_post_prior_l(data,f,1)))
}

# to get E(log p(Y | l,f))
e_log_lik = function(data,f){
    -0.5 * sum(log((2*pi)/f$tau) + f$tau * get_R2(data,f) * !data$missing)
}



#' compute KL divergence from q (posterior approximation) to prior
#' for kth loading
KL_post_prior_f = function(data, f,k=1){

  # compute x and s
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0
  s = as.vector(sqrt(1/(t(tau) %*% f$EL2[,k])))
  Rk = get_Rk(data,f,k) #residuals excluding factor k
  x = as.vector((t(Rk*tau) %*% f$EL[,k]) * s^2)

  KL_NM(x,s,f$comp_post_f[[k]],f$gf[[k]])

}

#' compute KL divergence from q (posterior approximation) to prior
#' for kth loading
KL_post_prior_l = function(data, f,k=1){

  # compute x and s
  tau = f$tau
  if(data$anyNA){tau = tau * !data$missing} #set missing values to have precision 0
  s = as.vector(sqrt(1/(tau %*% f$EF2[,k])))
  Rk = get_Rk(data,f,k) #residuals excluding factor k
  x = as.vector(((Rk*tau) %*% f$EF[,k]) * s^2)

  KL_NM(x,s,f$comp_post_l[[k]],f$gl[[k]])
}

#' @title compute the normal means part of the objective
#' @param x vector of observed normal means data
#' @param s vector of standard errors
#' @param comp_post a list with elements (prob,mean,mean2) each of which is a matrix of
#' posterior probabilities, eg as returned from ashr::ash by flash_data (outputlevel=5)
KL_NM = function(x,s,comp_post,g){
  -0.5 * ( log(2*pi*s^2) +
          (1/s^2) * rowSums(comp_post$prob * (comp_post$mean2 - 2*x*comp_post$mean)) +
          x^2/s^2 ) -
          ashr::calc_vloglik(g,ashr::set_data(x,s))
}

KL_NM2 = function(x,s,f,g){
  -0.5 * (log(2*pi*s^2) + (1/s^2) * (f$EL2[,1] - 2*x*f$EL[,1] + x^2)) -
    ashr::calc_vloglik(g,ashr::set_data(x,s))
}

#' @title expected loglikelihood for normal means model (expectation taken over posterior on parameters)
#' @param x observations in normal means
#' @param s standard errors of x
#' @param Em posterior mean of mean for x
#' @param Em2 posterior second moment of mean for x
NM_posterior_eloglik = function(x,s,Em,Em2){
  -0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Em2 - 2*x*Em + x^2))
}

#' @title Get objective function for flash data and fit objects
#' @export
get_F2 = function(data,f){
  return(sum(unlist(f$KL_l))+sum(unlist(f$KL_f))+e_log_lik(data,f) )
}


