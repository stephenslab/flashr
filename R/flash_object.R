construct_flash_object = function(data, fit, history, f_init) {
  flash_object = list()

  flash_object$nfactors = flash_get_nfactors(fit)
  flash_object$pve = flash_get_pve(fit)
  flash_object$fitted_values = flash_get_fitted_values(fit)
  flash_object$ldf = flash_get_ldf(fit)
  flash_object$objective = flash_get_objective(data, fit)
  if (class(f_init) == "flash") {
    flash_object$history = c(f_init$history, history)
  } else {
    flash_object$history = history
  }
  flash_object$fit = fit

  class(flash_object) = "flash"

  return(flash_object)
}

get_flash_fit = function(flash) {
  return(flash$fit)
}
