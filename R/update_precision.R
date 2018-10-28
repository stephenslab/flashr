# @title Update precision parameter
#
# @description Updates the estimated precision to increase the value of
#   the objective function.
#
# @inheritParams flash
#
# @param f A flash object.
#
# @return An updated flash object.
#
flash_update_precision = function(data,
                                  f,
                                  var_type) {
  R2 = flash_get_R2(data, f)
  f$tau = compute_precision(R2, data, var_type)
  return(f)
}


compute_precision = function(R2, data, var_type) {
  if (data$anyNA) {
    R2[data$missing] = NA
  }

  if (var_type == "by_column") {
    tau = mle_precision_by_column(R2)
  }
  else if (var_type == "by_row") {
    tau = t(mle_precision_by_column(t(R2)))
  }
  else if (var_type == "constant") {
    tau = mle_precision_constant(R2)
  }
  else if (var_type == "zero") {
    tau = 1 / data$S^2
  }

  if (is.matrix(tau) && data$anyNA) {
    tau[data$missing] = 0
  }

  return(tau)
}


# @title MLE for precision (separate parameter for each column)
#
# @param R2 An n by p matrix of squared residuals (with NAs for missing).
#
# @return An n by p matrix of precisions (separate value for each column).
#
mle_precision_by_column = function (R2) {
  sigma2 = colMeans(R2, na.rm = TRUE)  # a p vector

  # If a value of tau becomes numerically negative, set it to a
  # small positive number.
  tau = pmax(1/sigma2, .Machine$double.eps)
  return(outer(rep(1, nrow(R2)), tau))
}


# @title MLE for precision (single value)
#
# @param R2 An n by p matrix of squared residuals (with NAs for missing).
#
# @return An n by p matrix of precisions (a single value).
#
mle_precision_constant = function(R2) {
  sigma2 = mean(R2, na.rm = TRUE)  # a scalar

  tau = pmax(1/sigma2, .Machine$double.eps)
  return(tau)
}
