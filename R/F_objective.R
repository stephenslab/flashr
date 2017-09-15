
#currently just does it for K=1 factor
get_F = function(data,f){
  return(e_log_lik(data,f) - sum(KL_post_prior_f(data,f,1)) - sum(KL_post_prior_l(data,f,1)))
}

# to get E(log p(Y | l,f))
e_log_lik = function(data,f){
    -0.5 * sum(log(f$tau/(2*pi)) + f$tau * get_R2(data,f) * !data$missing)
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

  -0.5 * ( log(2*pi*s^2) +
    (1/s^2) * rowSums(f$comp_post_f[[k]]$prob *
                        (f$comp_post_f[[k]]$mean2 -
                           2*x*f$comp_post_f[[k]]$mean)) +
    x^2/s^2 ) - ashr::calc_vloglik(f$gf[[k]],ashr::set_data(x,s))

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

  -0.5 * ( log(2*pi*s^2) +
             (1/s^2) * rowSums(f$comp_post_l[[k]]$prob *
                                 (f$comp_post_l[[k]]$mean2 -
                                    2*x*f$comp_post_l[[k]]$mean)) +
             x^2/s^2 ) - ashr::calc_vloglik(f$gl[[k]],ashr::set_data(x,s))

}

