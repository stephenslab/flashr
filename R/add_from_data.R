#' @title Add factors to a flash object based on data
#'
#' @description Computes the current residuals from \code{data} and
#'   \code{f_init} and adds \code{K} new factors based on \code{init_fn}
#'   applied to these residuals. (If \code{f_init} is \code{NULL} then
#'   the residuals are the data.)
#'
#' @inheritParams flash
#'
#' @param K The number of factors to add.
#'
#' @param ... Additional parameters to be passed to \code{flash_backfit}.
#'
#' @export
#'
flash_add_factors_from_data = function(data,
                                       K,
                                       f_init = NULL,
                                       init_fn = "udv_si",
                                       backfit = TRUE,
                                       ...) {
  f = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f)
  init_fn = handle_init_fn(init_fn)

  f = add_factors_from_data(data, K, f, init_fn)

  history = NULL
  if (backfit) {
    flash_object = flash_backfit(data, f, ...)
    f = get_flash_fit(flash_object)
    history = get_flash_fit_history(flash_object)
  }

  flash_object = construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init,
                                        compute_obj = backfit)

  return(flash_object)
}

# More efficient "private" function used by flash_add_greedy
#
add_factors_from_data = function(data,
                                 K,
                                 f,
                                 init_fn) {
  R = flash_get_R_withmissing(data, f)
  f2 = flash_init_fn(flash_set_data(R), init_fn, K)
  f = flash_combine(f, f2)

  return(f)
}
