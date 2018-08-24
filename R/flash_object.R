init_flash_object_from_fit = function(fit) {
  flash_object = list(nfactors = NA,
                      pve = NA,
                      fitted_values = NA,
                      ldf = NA,
                      objective = NA,
                      history = list(),
                      fit = fit)
  class(flash_object) = "flash"
  return(flash_object)
}

update_flash_object = function(flash_object, data, fit, history) {
  flash_object$nfactors = flash_get_nfactors(fit)
  flash_object$pve = flash_get_pve(fit)
  flash_object$fitted_values = flash_get_fitted_values(fit)
  flash_object$ldf = flash_get_ldf(fit)
  flash_object$objective = flash_get_objective(data, fit)
  flash_object$history = c(flash_object$history, history)
  flash_object$fit = fit
  return(flash_object)
}

get_flash_fit = function(flash) {
  return(flash$fit)
}

set_flash_fit = function(flash, fit) {
  flash$fit = fit
  return(flash)
}

add_flash_history = function(flash, history) {
  flash$history = c(flash$history, list(history))
  return(flash)
}
