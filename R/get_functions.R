# functions for extracting useful information about the object

#' @title Return the estimated LF' matrix
#' @param f a flash fit object
#' @return  the estimated value of LF', an n by p matrix
#' @export
flash_get_lf = function(f){
  return(f$EL %*% t(f$EF))
}

#' @title  Get the residuals from a flash data and fit object,
#' excluding factor k
get_Rk = function(data,f,k){
  return( data$Y - f$EL[,-k,drop=FALSE] %*% t(f$EF[,-k,drop=FALSE]) )
}

#' @title  Get the residuals from data and a flash fit object
get_R = function(data,f){
  return( data$Y - flash_get_lf(f) )
}


#' @title  Get the expected squared residuals from a flash data and fit object
get_R2 = function(data,f){
  LF = f$EL %*% t(f$EF)
  return( data$Y^2 - 2 * data$Y * LF
          + LF^2 + f$EL2 %*% t(f$EF2) - f$EL^2 %*% t(f$EF^2) )
}

#' @title is_tiny_fl
#' @details checks whether kth factor/loading combination is tiny
is_tiny_fl=function(f,k,tol=1e-8){
    return(sum(f$EL[,k]^2)*sum(f$EF[,k]^2) < tol)
}

get_l = function(f){f$EL}
get_f = function(f){f$EF}
get_k = function(f){ncol(f$EL)}
get_n = function(f){nrow(f$EL)}
get_p = function(f){nrow(f$EF)}

#' @title flash_get_sizes
#' @details computes the factor contributions (analogous to eigenvalues)
#' for each factor in flash fit f
#' @export
flash_get_sizes = function(f){
  colSums(f$EL^2)*colSums(f$EF^2)
}

get_conv_criteria = function(f){flash_get_lf(f)}


flash_get_LB = function(data,f){
  return(get_c_likelihood(data,f)+get_prior_post(f))
}

# to get the conditional likelihood given L and F
get_c_likelihood = function(data,f){
  # to get residual square matrix
  residual_square =  data$Y^2 - 2*data$Y*(f$EL %*% t(f$EF)) + (f$EL2 %*% t(f$EF2))
  c_lik = -(1/2) * sum( - log(2*pi*f$tau) + (residual_square)*f$tau , na.rm = TRUE)
  return(c_lik)
}

# to get the posterior and prior part
get_prior_post = function(f){
  K = get_k(f)
  prior_post_val = 0
  for(k in 1:K){
    prior_post_val = prior_post_val + get_post_prior_single(f,k)
  }
  return(prior_post_val)
}

get_post_prior_single = function(f,k){
  # I use the ash_para for l and f separately for further use of different famliy of prior
  ash_para_l = list(method = "shrink", mixcompdist = "normal")
  ash_para_f = list(method = "shrink", mixcompdist = "normal")
  res_l = post_prior_func(f$matl[[k]], f$gl[[k]],ash_para_l)
  res_f = post_prior_func(f$matf[[k]], f$gf[[k]],ash_para_f)
  priopost_f = res_f$PrioPost
  priopost_l = res_l$PrioPost
  penalty_f = res_f$penalty
  penalty_l = res_l$penalty
  return(priopost_l + priopost_f + penalty_l + penalty_f)
}

post_prior_func = function(mat,fit_g,ash_para){
  if(is.null(ash_para$method)){
    # in this case the method is default "fdr"
    # make sure the first row included in the index
    nonzeroindex = unique(c(1,which(fit_g$pi != 0)))
  }else{
    # there is no zero in the shrink method
    # here we mean the shrink method without point mass!
    nonzeroindex = c(which(fit_g$pi != 0))
  }
  #prepare for the posterior part
  mat_postmean = mat$comp_postmean[nonzeroindex,]
  mat_postmean2 = mat$comp_postmean2[nonzeroindex,]
  mat_postvar = mat_postmean2 - mat_postmean^2
  mat_postprob = mat$comp_postprob[nonzeroindex,]
  # figure our the dimension
  K = dim(mat_postprob)[1]
  N = dim(mat_postprob)[2]
  if(is.vector(mat_postprob)){
    K = 1
    N = length(mat_postprob)
  }
  # prepare for the prior
  prior_pi = fit_g$pi[nonzeroindex]
  prior_var = (fit_g$sd[nonzeroindex])^2
  mat_priorvar = matrix(rep(prior_var,N),ncol = N)
  mat_priorprob = matrix(rep(prior_pi,N),ncol = N)
  # to get the value
  varodds = mat_priorvar / mat_postvar
  if(is.null(ash_para$method)){
    # only in the fdr case, which has zero
    varodds[1,] = 1 # deal with 0/0
  }
  probodds = mat_priorprob / mat_postprob
  ssrodds = mat_postmean2 / mat_priorvar
  if(is.null(ash_para$method)){
    # only in the fdr case, which has zero
    ssrodds[1,] = 1  # deal with 0/0
  }
  priorpost = mat_postprob * (log(probodds) - (1/2)*log(varodds) - (1/2)* (ssrodds -1))
  # in case of mat_postprob = 0
  priorpost[which(mat_postprob< 1e-100)] = 0
  PrioPost = sum(priorpost)

  #######
  # now we try to get the penalty term
  # I use the default of setting for lambda_k in https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/biostatistics/18/2/10.1093_biostatistics_kxw041/1/kxw041_Supp.pdf?Expires=1498001808&Signature=dtxDDFmYZ1u1jFoXfNg4mqGQnTPkQQMFZgt2-FE9r8GcgRl4DdEaP0oI6kzvUjZB5H3hLZmvAZfuFkt3GxojQJ1eX0y7I9kMJC75VTrxw1Ym~psIe8rdU5Xw4S4KZm79o1cOXBzY1iAORFkVM4z59uP1T1ltFJDu5fhlTYErU3flbiE1ivTm-BYN5P7w8dex-R5MWn98NpMVRUsoHRgbe7FhwavdnQ8QAD7~q1F4DzMyZk09TeNdEc7GNdBKbSOoS9mGGnvrqXvPz1nn~QzJ2U2CMCf04GfsgkclyUnH0aD6Vh1lSAfoEV47ZnT-9~942~f8rMs30H2Kqji7BjlctA__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q
  # this is a simple setting for fdr method we use the default lambda 10 for shrink we use nothing
  if(is.null(ash_para$method)){
    # in this case the method is default "fdr"
    #h_value = (10 - 1) * log(fit_g$pi[1])
    if(is.null(ash_para$nullweight)){
      h_value = (10 - 1) * log(fit_g$pi[1])
    }else{
      h_value = (ash_para$nullweight - 1) * log(fit_g$pi[1])
    }
  }else{
    # we assume that the only alternative method is shrink for now
    h_value = 0
  }
  return(list(PrioPost = PrioPost, penalty = h_value))
}

