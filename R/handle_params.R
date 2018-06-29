handle_kset = function(kset, f) {
  if (is.null(kset)) {
    kset = 1:flash_get_k(f)
  } else if (!is.numeric(kset) || max(kset) > flash_get_k(f)) {
    stop(paste("Invalid kset. Kset should be a vector containing the",
               "indices of the factors to be optimized."))
  }
  kset
}

handle_init_fn = function(init_fn) {
  # Either a function or the name of a function (as a string) are ok.
  if (!is.function(init_fn) && !exists(init_fn, mode="function")) {
    stop("The specified init_fn does not exist.")
  }
  init_fn
}

handle_ebnm_fn = function(ebnm_fn) {
  if (!is.list(ebnm_fn)) {
    if (length(ebnm_fn) != 1) {
      stop(paste("Either a single function or a list with fields l",
                 "and f must be specified for parameter ebnm_fn."))
    }
    ebnm_fn_l = ebnm_fn
    ebnm_fn_f = ebnm_fn
  } else if (xor(is.null(ebnm_fn$l), is.null(ebnm_fn$f))) {
    stop(paste("If ebnm_fn is specified for either loadings or factors",
               "then it must be specified for both."))
  } else if (!is.null(ebnm_fn$l)) {
    ebnm_fn_l = ebnm_fn$l
    ebnm_fn_f = ebnm_fn$f
  } else {
    stop("Invalid entry for parameter ebnm_fn.")
  }

  if (is.function(ebnm_fn_l) || is.function(ebnm_fn_f)) {
    stop(paste("Invalid entry for parameter ebnm_fn. Please supply a",
               "character string such as \"ebnm_pn\" or \"ebnm_ash\".",
               "If using a custom function then enclose the name of the",
               "function in quotes."))
  }

  if (ebnm_fn_l == "ebnm_pn" || ebnm_fn_f == "ebnm_pn") {
    if (!requireNamespace("ebnm", quietly = TRUE)) {
      message(paste("ebnm package not installed. ebnm_ash will be used",
                    "instead of ebnm_pn."))
      if (ebnm_fn_l == "ebnm_pn") {
        ebnm_fn_l = "ebnm_ash"
      }
      if (ebnm_fn_f == "ebnm_pn") {
        ebnm_fn_f = "ebnm_ash"
      }
    }
  }

  if (!exists(ebnm_fn_l, mode="function")
      || !exists(ebnm_fn_f, mode="function")) {
    stop("The specified ebnm function does not exist.")
  }

  list(l = ebnm_fn_l, f = ebnm_fn_f)
}

handle_ebnm_param = function(ebnm_param, ebnm_fn, n_expected) {
  # Check to see whether parameters are specified separately for loadings
  #   and factors:
  if (xor(is.null(ebnm_param$l), is.null(ebnm_param$f))) {
    stop(paste("if ebnm_param is specified for either loadings or",
               "factors then it must be specified for both. (Use an",
               "empty list to specify no parameters.)"))
  } else if (!is.null(ebnm_param$l)) {
    ebnm_param_l = ebnm_param$l
    ebnm_param_f = ebnm_param$f
  } else {
    ebnm_param_l = ebnm_param
    ebnm_param_f = ebnm_param
  }

  # Check to see whether there are different parameters for each
  #   subsequent loading/factor (n_expected of them are required):
  if (length(ebnm_param_l) == 0) {
    # NULL or empty list
    ebnm_param_l = rep(list(list()), n_expected)
  } else if (!is.list(ebnm_param_l[[1]])) {
    ebnm_param_l = rep(list(ebnm_param_l), n_expected)
  } else if (length(ebnm_param_l) < n_expected) {
    stop(paste("If different ebnm parameters are used for each loading",
               "then ebnm_param$l must be a list of", n_expected,
               "lists."))
  }
  if (length(ebnm_param_f) == 0) {
    ebnm_param_f = rep(list(list()), n_expected)
  } else if (!is.list(ebnm_param_f[[1]])) {
    ebnm_param_f = rep(list(ebnm_param_f), n_expected)
  } else if (length(ebnm_param_f) < n_expected) {
    stop(paste("If different ebnm parameters are used for each factor",
               "then ebnm_param$f must be a list of", n_expected,
               "lists."))
  }

  # Add defaults:
  ebnm_param_l = lapply(ebnm_param_l,
                        function(x) {add_ebnm_fn_defaults(x, ebnm_fn$l)})
  ebnm_param_f = lapply(ebnm_param_f,
                        function(x) {add_ebnm_fn_defaults(x, ebnm_fn$f)})
  list(l = ebnm_param_l, f = ebnm_param_f)
}

add_ebnm_fn_defaults = function(ebnm_param, ebnm_fn) {
  if (ebnm_fn == "ebnm_ash") {
    return(add_ebnm_ash_defaults(ebnm_param))
  } else if (ebnm_fn == "ebnm_pn") {
    return(ebnm_param)
  }
}

add_ebnm_ash_defaults = function(ebnm_param) {
  if (is.null(ebnm_param$mixcompdist)) {
    ebnm_param$mixcompdist = "normal"
  }
  if (is.null(ebnm_param$method)) {
    ebnm_param$method = "shrink"
  }
  ebnm_param
}
