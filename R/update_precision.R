#' @title Update precision parameter
#'
#' @description Updates the estimated precision to increase the value of
#'   the objective function.
#'
#' @inheritParams flash
#'
#' @param f A flash object.
#'
#' @return An updated flash object.
#'
#' @export
#'
flash_update_precision = function(data,
                                  f,
                                  var_type = c("by_column", "by_row",
                                               "constant", "zero",
                                               "kroneker")) {
  data = handle_data(data)
  var_type = match.arg(var_type)

  if (var_type == "zero" & is.null(data$S)) {
    stop(paste("Flash data object must include standard errors when",
               "var_type is zero."))
  }
  if (!is.null(data$S) & var_type != "zero") {
    stop("Standard errors are currently only used when var_type is zero.")
  }

  R2 = flash_get_R2(data, f)
  f$tau = compute_precision(R2, data$missing, var_type, data$S)
  return(f)
}


compute_precision = function(R2, missing, var_type, S) {
  R2[missing] = NA

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
    tau = 1 / S^2
  }
  else (stop("that var_type not yet implemented"))

  tau[missing] = 0
  return(tau)
}


# @title MLE or precision (separate parameter for each column).
#
# @param R2 n by p matrix of squared residuals (with NAs for missing).
#
# @return n by p matrix of precisions (separate value for each column).
#
mle_precision_by_column = function (R2) {
  sigma2 = colMeans(R2, na.rm = TRUE)  # A p vector.

  # If a value of tau becomes numerically negative, set it to a
  # small positive number.
  tau = pmax(1/sigma2, .Machine$double.eps)
  return(outer(rep(1, nrow(R2)), tau)) # An n by p matrix.
}


# @title mle for precision (separate parameter for each column).
#
# @param R2 n by p matrix of squared residuals (with NAs for missing).
#
# @return n by p matrix of precisions (separate value for each column).
#
mle_precision_constant = function(R2) {
  sigma2 = mean(R2, na.rm = TRUE)  # a scalar

  # If a value of tau becomes numerically negative, set it to a
  # small positive number.
  tau = pmax(1/sigma2, .Machine$double.eps)
  return(matrix(tau, nrow = nrow(R2), ncol = ncol(R2))) # An n by p matrix.
}
