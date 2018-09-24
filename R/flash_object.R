construct_flash_object = function(data,
                                  fit,
                                  history,
                                  f_init,
                                  compute_obj = TRUE) {
  flash_object = list()

  summary_stats = compute_summary_statistics(data, fit)
  flash_object$nfactors = summary_stats$nfactors
  flash_object$pve = summary_stats$pve
  flash_object$fitted_values = summary_stats$fitted_values
  flash_object$ldf = summary_stats$ldf

  if (compute_obj) {
    flash_object$objective = flash_get_objective(data, fit)
  } else {
    flash_object$objective = NA
  }

  if (class(f_init) == "flash") {
    flash_object$fit_history = c(f_init$fit_history, history)
  } else {
    flash_object$fit_history = history
  }

  flash_object$fit = fit

  class(flash_object) = "flash"

  return(flash_object)
}

get_flash_fit = function(flash) {
  return(flash$fit)
}

get_flash_fit_history = function(flash) {
  return(flash$fit_history)
}

compute_summary_statistics = function(data, f) {
  ldf = flash_get_ldf(f)

  d = ldf$d
  nfactors = length(d)

  s = d^2
  if (is.null(f$tau)) {
    var_from_tau = 0
  } else if (is.matrix(f$tau)) {
    tau = f$tau[f$tau != 0]
    var_from_tau = sum(1/tau)
  } else { # tau is a scalar
    var_from_tau = sum(!data$missing) / f$tau
  }
  pve = s/(sum(s) + var_from_tau)

  fitted_values = flash_get_fitted_values(f)

  return(list(ldf = ldf, nfactors = nfactors, pve = pve,
              fitted_values = fitted_values))
}
