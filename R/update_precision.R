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
flash_update_precision = function(data, f, var_type) {
  f$tau = compute_precision(flash_get_R2(data, f), data, var_type)
  return(f)
}

compute_precision = function(R2, data, var_type) {
  # We assume that R2 correctly has NAs where data is missing.

  tau = switch(var_type,
               by_column = mle_precision_by_column(R2, data),
               by_row = mle_precision_by_row(R2, data),
               constant = mle_precision_constant(R2, data),
               zero = set_precision_zero(data))

  return(tau)
}

# If there is no missing data, store precision as a 1 x p matrix. Otherwise,
#   store as a n x p matrix with zeroes where data is missing.
#
mle_precision_by_column = function(R2, data) {
  sigma2 = colMeans(R2, na.rm = TRUE)
  # If a value of tau becomes numerically negative, set it to a
  #   small positive number:
  tau = pmax(1 / sigma2, .Machine$double.eps)

  if (data$anyNA) {
    tau = matrix(tau, nrow = nrow(R2), ncol = ncol(R2), byrow = TRUE)
    tau[data$missing] = 0
  } else {
    tau = matrix(tau, nrow = 1)
  }

  return(tau)
}

# If there is no missing data, store precision as a n x 1 matrix. Otherwise,
#   store as a n x p matrix with zeroes where data is missing.
#
mle_precision_by_row = function(R2, data) {
  sigma2 = rowMeans(R2, na.rm = TRUE)
  tau = pmax(1 / sigma2, .Machine$double.eps)

  if (data$anyNA) {
    tau = matrix(tau, nrow = nrow(R2), ncol = ncol(R2), byrow = FALSE)
    tau[data$missing] = 0
  } else {
    tau = matrix(tau, ncol = 1)
  }

  return(tau)
}

# If there is no missing data, store precision as a scalar. Otherwise,
#   store as a n x p matrix with zeroes where data is missing.
#
mle_precision_constant = function(R2, data) {
  sigma2 = mean(R2, na.rm = TRUE)
  tau = max(1 / sigma2, .Machine$double.eps)

  if (data$anyNA) {
    tau = matrix(tau, nrow = nrow(R2), ncol = ncol(R2))
    tau[data$missing] = 0
  }

  return(tau)
}

set_precision_zero = function(data) {
  tau = 1 / data$S^2

  if (data$anyNA) {
    if (!is.matrix(tau)) {
      tau = matrix(tau, nrow = nrow(data$missing), ncol = ncol(data$missing))
    }
    tau[data$missing] = 0
  }

  return(tau)
}
