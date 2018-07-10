# @title Handle kset parameter
#
# @description Checks that kset is numeric and not out of bounds for the
#   flash object. Defaults to kset = 1:flash_get_k(f).
#
# @param kset A vector of factor indices.
#
# @param f A fitted flash object.
#
# @return kset (with values possibly defaulted in).
#
handle_kset = function(kset, f) {
  if (is.null(kset)) {
    # Default:
    kset = 1:flash_get_k(f)
  } else if (!is.numeric(kset) || max(kset) > flash_get_k(f)) {
    stop(paste("Invalid kset. Kset should be a vector containing the",
               "indices of the factors to be optimized."))
  }
  kset
}


# @title Handle init_fn parameter
#
# @description Checks that init_fn is a valid function.
#
# @param init_fn An initialization function. Either the name of a
#   function or the function itself (as a character string) are
#   acceptable arguments.
#
# @return init_fn
#
handle_init_fn = function(init_fn) {
  if (!is.function(init_fn) && !exists(init_fn, mode="function")) {
    stop("The specified init_fn does not exist.")
  }
  init_fn
}


# @title Handle ebnm_fn parameter
#
# @description Checks that the argument to ebnm_fn refers to a valid
#   function. If ebnm_pn is used, checks that package ebnm is installed.
#   If not, returns "ebnm_ash" instead (with a message to the user).
#
# @param ebnm_fn The function used to solve the empirical Bayes normal
#   means problem. Either a single character string (giving the name of
#   of the function) or a list with fields l and f (specifying
#   different functions to be used for loadings and factors) are
#   acceptable arguments.
#
# @return ebnm_fn as a list with fields l and f (if a single function is
#   passed in, then it will be the value for both fields).
#
handle_ebnm_fn = function(ebnm_fn) {
  if (!is.list(ebnm_fn)) {
    if (length(ebnm_fn) != 1) {
      stop(paste("Either a single character string or a list with fields",
                 "l and f must be specified for parameter ebnm_fn."))
    }
    ebnm_fn_l = ebnm_fn
    ebnm_fn_f = ebnm_fn
  } else if (is.null(ebnm_fn$l) || is.null(ebnm_fn$f)) {
    stop(paste("If a list is specified for parameter ebnm_fn, then it",
               "must include fields l and f."))
  } else {
    ebnm_fn_l = ebnm_fn$l
    ebnm_fn_f = ebnm_fn$f
  }

  if (is.function(ebnm_fn_l) || is.function(ebnm_fn_f)) {
    stop(paste("Invalid argument for parameter ebnm_fn. Please supply a",
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


# @title Handle ebnm_param parameter
#
# @description Checks that the argument to ebnm_param makes sense and
#   adds default parameters when available.
#
# @param ebnm_param The parameters to be passed into ebnm_fn. Several
#   types of argument are possible. NULL will return default parameters
#   (or empty lists when no defaults are available). A single named list
#   will supply the parameters for both EBNM functions. A list with
#   fields l and f will separately supply parameters for the functions
#   used for loadings and factors. A list of n named lists will
#   separately specify parameters for each factor/loading combination.
#   Finally, the latter two types of argument can be combined, so that
#   the user may supply a list with fields l and f, each of which is a
#   list of n named lists.
#
# @param ebnm_fn A list with fields l and f specifying the functions
#   used to solve the EBNM problem for loadings and factors. The output
#   of handle_ebnm_fn will be the usual input here.
#
# @param n_expected The number of lists to expect if separate parameters
#   are used for each factor/loading combination. (That is, the number
#   of distinct factor/loadings to which ebnm_fn will be applied.)
#
# @return ebnm_param as a list with fields l and f, each of which is a
#   list of n_expected lists. When available, default parameters are
#   added.
#
handle_ebnm_param = function(ebnm_param, ebnm_fn, n_expected) {
  if (!is.null(ebnm_param) && !is.list(ebnm_param)) {
    stop(paste("Invalid argument for parameter ebnm_param. A list (or",
               "NULL) was expected."))
  }

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
  }
  # If the list is named, it gives parameters for all loadings:
  else if (!is.null(names(ebnm_param_l))) {
    ebnm_param_l = rep(list(ebnm_param_l), n_expected)
  }
  # And an unnamed list gives parameters separately for each loading:
  else if (length(ebnm_param_l) != n_expected) {
    stop(paste("If different ebnm parameters are used for each loading",
               "then ebnm_param$l must be a list of", n_expected,
               "lists."))
  }
  if (length(ebnm_param_f) == 0) {
    ebnm_param_f = rep(list(list()), n_expected)

  } else if (!is.null(names(ebnm_param_f))) {
    ebnm_param_f = rep(list(ebnm_param_f), n_expected)
  }
  else if (length(ebnm_param_f) != n_expected) {
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

# @title Add defaults to ebnm_param
#
# @description Adds default parameters to ebnm_param, depending on
#   which function is used for the EBNM problem. At present, defaults
#   are only added when ebnm_ash is used.
#
# @param ebnm_param The parameters to be passed into ebnm_fn, given as
#   a single list.
#
# @param ebnm_fn The name of the function used to solve the EBNM
#   problem, given as a character string.
#
# @return ebnm_param, with defaults added when available.
#
add_ebnm_fn_defaults = function(ebnm_param, ebnm_fn) {
  if (ebnm_fn == "ebnm_ash") {
    return(add_ebnm_ash_defaults(ebnm_param))
  } else if (ebnm_fn == "ebnm_pn") {
    return(ebnm_param)
  }
  ebnm_param
}

# @title Add ebnm_ash defaults to ebnm_param
#
# @description Adds default ebnm_ash parameters to ebnm_param.
#
# @param ebnm_param The parameters to be passed into ebnm_fn, given as
#   a single list.
#
# @return ebnm_param, with ash defaults added (when not specified by the
#   user).
#
add_ebnm_ash_defaults = function(ebnm_param) {
  ebnm_param$output = "flash_data"
  if (is.null(ebnm_param$mixcompdist)) {
    ebnm_param$mixcompdist = "normal"
  }
  if (is.null(ebnm_param$method)) {
    ebnm_param$method = "shrink"
  }
  ebnm_param
}
