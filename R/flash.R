#' @title Fit the rank1 EBFA model to data
#' @param Y an n by p data matrix
#' @param tol specify how much objective can change in a single iteration to be considered not converged
#' @return a fitted flash object
#' @export
flash_r1 = function(Y,init_method=c("svd","random"),tol=1e-2){
  init_method=match.arg(init_method)

  f = flash_set_data(Y)
  f = flash_init(f,1,init_method)
  c = get_conv_criteria(f)
  diff = 1
  while(diff > tol){
    f = flash_update_single_factor(f,1)
    f = flash_update_single_loading(f,1)
    f = flash_update_precision(f)
    cnew = get_conv_criteria(f)
    diff = mean((cnew-c)^2)
    c = cnew
  }
  return(f)
}

