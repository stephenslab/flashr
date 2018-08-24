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

  R = flash_get_R_withmissing(data, f)
  f2 = flash_init_fn(flash_set_data(R), init_fn, K)
  f = flash_combine(f, f2)

  history = NULL
  if (backfit) {
    flash_object = flash_backfit(data, f, ...)
    f = flash_object$fit
    history = flash_object$history
  }

  flash_object = construct_flash_object(data = data,
                                        fit = f,
                                        history = history,
                                        f_init = f_init,
                                        compute_obj = backfit)

  return(flash_object)
}


#' @title Add factor/loading pairs to a flash object
#'
#' @description Adds specified factor/loading pairs to a flash object.
#'
#' @inheritParams flash
#'
#' @param LL The loadings, an n by K matrix.
#'
#' @param FF The factors, a p by K matrix.
#'
#' @param fixl An n by K matrix of \code{TRUE}/\code{FALSE} values
#'   indicating which elements of \code{LL} should be considered fixed
#'   and not changed during updates. Useful for including a mean factor
#'   for example.
#'
#' @param fixf A p by K matrix of \code{TRUE}/\code{FALSE} values; same
#'   as \code{fixl} but for factors \code{FF}.
#'
#' @return A flash fit object, with additional factors initialized
#'   using \code{LL} and \code{FF}.
#'
#' @export
#'
flash_add_lf = function(data,
                        LL,
                        FF,
                        f_init = NULL,
                        fixl = NULL,
                        fixf = NULL) {
  f_init = handle_f(f_init, init_null_f = TRUE)
  data = handle_data(data, f_init)
  LL = handle_LL(LL, expected_nrow = flash_get_n(f_init))
  FF = handle_LL(FF, expected_nrow = flash_get_p(f_init))
  # fixl and fixf are handled by flash_init_lf

  f2 = flash_init_lf(LL, FF, fixl, fixf)
  f = flash_combine(f_init, f2)

  return(f)
}
