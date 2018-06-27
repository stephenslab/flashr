flash_default_ebnm_param = function(ebnm_fn) {
    if (identical(ebnm_fn, ebnm_ash)) {
        return(list(output      = "flash_data",
                    mixcompdist = "normal",
                    method      = "shrink"))
    } else if (identical(ebnm_fn, ebnm_pn)) {
        return(list())
    } else {
        stop(paste("No defaults available for ebnm_param for that ebnm",
                   "function---please supply them"))
    }
}

handle_ebnm_fn = function(ebnm_fn) {
  if (xor(is.null(ebnm_fn$l), is.null(ebnm_fn$f))) {
    stop("if ebnm_fn is specified for l then it must also be specified for f")
  } else if (!is.null(ebnm_fn$l)) {
    ebnm_fn_l = ebnm_fn$l
    ebnm_fn_f = ebnm_fn$f
  } else {
    ebnm_fn_l = ebnm_fn
    ebnm_fn_f = ebnm_fn
  }
  list(l = ebnm_fn_l, f = ebnm_fn_f)
}

handle_ebnm_param = function(ebnm_param, k) {
  # Check to see whether parameters are specified separately for loadings
  #   and factors:
  if (xor(is.null(ebnm_param$l), is.null(ebnm_param$f))) {
    stop("if ebnm_param is specified for l then it must also be specified for f")
  } else if (!is.null(ebnm_param$l)) {
    ebnm_param_l = ebnm_param$l
    ebnm_param_f = ebnm_param$f
  } else {
    ebnm_param_l = ebnm_param
    ebnm_param_f = ebnm_param
  }
  # Check to see whether there are different parameters for each
  #   subsequent loading/factor (k of them are required):
  if (!is.list(ebnm_param_l[[1]])) {
    ebnm_param_l = rep(list(ebnm_param_l), k)
  }
  if (!is.list(ebnm_param_f[[1]])) {
    ebnm_param_f = rep(list(ebnm_param_f), k)
  }
  list(l = ebnm_param_l, f = ebnm_param_f)
}
