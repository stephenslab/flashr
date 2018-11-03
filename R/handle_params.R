# Argument checks for main flash function -------------------------------


valid_var_types = c("by_column", "by_row", "constant", "zero")
valid_verbose_options = "odLFlfn"

# Handle flash object parameter
#
# @description Checks that input to f is valid and initializes a null
#   flash fit object when necessary.
#
# @param f A flash fit object.
#
# @param allow_null If TRUE, accepts NULL as a valid argument.
#
# @param init_null_f If TRUE and f is NULL, initializes a flash object
#   via a call to flash_init_null().
#
# @return f The flash object.
#
handle_f = function(f, allow_null = TRUE, init_null_f = FALSE) {
  if (!is.null(f) && is(f, "flash")) {
    return(get_flash_fit(f))
  }
  if (!allow_null && !is(f, "flash_fit")) {
    stop("f must be a flash object or a flash fit object.")
  }
  if (!is.null(f) && !is(f, "flash_fit")) {
    stop("f must be NULL, a flash object, or a flash fit object.")
  }
  if (init_null_f && is.null(f)) {
    f = flash_init_null()
  }

  return(f)
}


# Check that inputs to Y are valid
#
# @param Y A data matrix.
#
# @return Y
#
handle_Y = function(Y) {
  if(!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a matrix or a flash data object.")
  }
  if (any(is.infinite(Y))) {
    stop("Y must not contain infinite values.")
  }
  if (any(is.nan(Y))) {
    stop("Y must not contain NaNs.")
  }

  return(Y)
}


# Check that inputs to S are valid (given Y)
#
# @param S A scalar or matrix of standard errors for Y. Can be NULL.
#
# @param Y A data matrix.
#
# @return A valid flash data object.
#
handle_S = function(S, Y) {
  if (is.null(S)) {
    return(S)
  }

  if (!is.numeric(S)) {
    stop("Invalid input to S.")
  }
  if (length(S) > 1 && (!is.matrix(S) || !identical(dim(Y), dim(S)))) {
    stop("S must be a scalar or a matrix of the same dimensions as Y.")
  }
  if (requireNamespace("ebnm", quietly = TRUE) &&
      packageVersion("ebnm") < "0.1.13") {
    # Earlier versions of ebnm do not support scalar arguments for S
    S = matrix(S, nrow = nrow(Y), ncol = ncol(Y))
  }

  if (is.matrix(S)) {
    S[is.na(Y)] = Inf
  }

  if (any(is.na(S))) {
    stop("S must not contain NAs where Y has data.")
  }
  if (any(is.nan(S))) {
    stop("S must not contain NaNs.")
  }
  if (any(S < 0)) {
    stop("S must not contain negative values.")
  }

  return(S)
}


# Handle data parameter
#
# @description Checks that data object is valid (given f).
#
# @param data An n by p matrix or a flash data object.
#
# @param f A flash fit object.
#
# @param output If "flash_data", returns a flash data object. If
#   "matrix", returns a matrix.
#
# @return A matrix or flash data object.
#
handle_data = function(data, f, output = "flash_data") {
  if (!is.matrix(data) && !is(data, "flash_data")) {
    stop("Data must be a matrix or a flash data object.")
  }

  if (is.matrix(data) && output == "flash_data") {
    data = flash_set_data(data)
  }
  if (is(data, "flash_data") && output == "matrix") {
    data = get_Yorig(data)
  }

  if (!is.null(f)) {
    if (is.matrix(data)) {
      n = nrow(data)
      p = ncol(data)
    } else {
      n = nrow(data$Y)
      p = ncol(data$Y)
    }

    if ((!is.null(f$EL) && nrow(f$EL) != n) ||
        (!is.null(f$EF) && nrow(f$EF) != p) ||
        (!is.null(f$tau) && is.matrix(f$tau) &&
         !identical(dim(f$tau), c(n, p)))) {
      stop("Dimensions of data and f do not agree.")
    }
  }

  return(data)
}


# Handle backfit parameter
#
# @param backfit Can be a logical value or a vector containing the indices
#   of the factor/loading pairs to backfit.
#
# @return backfit
#
handle_backfit = function(backfit) {
  if (!(is.logical(backfit) && length(backfit) == 1)
      && !(is.numeric(backfit) && is.vector(backfit) && all(backfit > 0))) {
    stop(paste("backfit must be TRUE/FALSE or a vector containing the",
               "indices of the factor/loading pairs to backfit."))
  }

  if (identical(backfit, FALSE)) {
    return(numeric(0))
  }
  if (identical(backfit, TRUE)) {
    return(NULL) # will be subsequently set to 1:flash_get_k(fl)
  }

  return(backfit)
}


# Check that choice of var_type is valid (given data)
#
# @param var_type
#
# @param data A flash data object.
#
# @return var_type
#
handle_var_type = function(var_type, data) {
  if (!(var_type %in% valid_var_types)) {
    stop("That var_type has not yet been implemented.")
  }

  if (var_type == "zero" & is.null(data$S)) {
    stop(paste("Flash data object must include standard errors when",
               "var_type is zero."))
  }
  if (!is.null(data$S) & var_type != "zero") {
    stop("Standard errors are currently only used when var_type is zero.")
  }

  return(var_type)
}


# Handle fixed_loadings and fixed_factors
#
# @param fixed A matrix (or vector) of loadings or a list that includes
#   fields vals and is_fixed.
#
# @param expected_nrow The expected number of rows in the loadings matrix.
#
# @return A list with fields vals, is_fixed, and K.
#
handle_fixed = function(fixed, expected_nrow) {
  if (is.null(fixed)) {
    return(list(K = 0))
  }

  if (is.numeric(fixed)) {
    fixed = list(vals = fixed)
  } else if (!is.list(fixed)) {
    stop(paste("If nonnull, then fixed_ parameters must be matrices or",
               "lists."))
  }

  if (is.null(fixed$vals)) {
    stop(paste("If fixed_ parameter is a list, then it must include list",
               "element vals."))
  }

  fixed$vals = handle_LL(fixed$vals, expected_nrow)
  fixed$is_fixed = handle_fix(fixed$is_fixed,
                              fixed$vals,
                              default_val = TRUE)
  fixed$K = ncol(fixed$vals)

  return(fixed)
}


# Handle LL parameter (and FF)
#
# @description Checks that LL has the correct dimensions and converts
#   it from a vector to a matrix if necessary.
#
# @param LL A matrix of loadings. A vector can also be passed in (and is
#   treated as a single loading).
#
# @param expected_nrow The expected number of rows (can be NULL if no
#   flash object has yet been initialized).
#
# @return LL
#
handle_LL = function(LL, expected_nrow) {
  if (is.vector(LL)) {
    LL_names = names(LL)
    LL = matrix(LL, ncol = 1)
    rownames(LL) = LL_names
  }

  if (!is.numeric(LL) || !is.matrix(LL)) {
    stop("Invalid matrix of loadings/factors.")
  }
  if (any(is.infinite(LL))) {
    stop("The matrix of loadings/factors must not include infinite values.")
  }
  if (any(is.nan(LL))) {
    stop("The matrix of loadings/factors must not include NaNs.")
  }

  if (!is.null(expected_nrow) && nrow(LL) != expected_nrow) {
    stop(paste("The matrix of loadings/factors does not have the",
               "correct dimensions."))
  }

  return(LL)
}


# Handle parameter fixl (and fixf)
#
# @description Checks that fixl has the correct dimensions and converts
#   it from a vector to a matrix if necessary.
#
# @param fixl A matrix of TRUE/FALSE values indicating which entries of
#   LL are to be considered fixed.
#
# @param LL A matrix of loadings.
#
# @param default_val Indicates whether fixl should default to TRUE or
#   FALSE. If TRUE, then only NA values will not be considered fixed.
#
# @return fixl
#
handle_fix = function(fixl, LL, default_val) {
  if (default_val == FALSE && is.null(fixl)) {
    fixl = matrix(FALSE, nrow = nrow(LL), ncol = ncol(LL))
  }
  if (default_val == TRUE && is.null(fixl)) {
    fixl = !is.na(LL)
  }

  if (!identical(length(fixl), length(LL))) {
    stop("The dimensions of LL/FF and fixl/fixf do not match.")
  }

  fixl = matrix(as.logical(fixl), nrow = nrow(LL), ncol = ncol(LL))

  # If user tries to fix a NA, just unfix it:
  fixl[is.na(LL)] = FALSE

  return(fixl)
}


# Handle verbose parameter
#
# @description Converts verbose (which can be TRUE/FALSE or a character
#   string) to a vector of characters indicating which columns should be
#   printed during progress updates.
#
# @param verbose
#
# @return verbose_output
#
handle_verbose = function(verbose, stopping_rule) {
  if (identical(verbose, TRUE)) {
    verbose_output = get_default_verbose_output(stopping_rule)
  } else if (identical(verbose, FALSE)) {
    verbose_output = ""
  } else if (!is.character(verbose)) {
    stop("verbose must be TRUE/FALSE or a character string.")
  } else {
    verbose_output = verbose
  }

  verbose_output = unlist(strsplit(verbose_output, split = NULL))
  valid_options = unlist(strsplit(valid_verbose_options, split = NULL))
  if (!all(verbose_output %in% valid_options)) {
    stop("At least one verbose output option is unrecognized.")
  }

  return(verbose_output)
}


handle_output_level = function(output_level) {
  summary_stats = c("nfactors", "pve", "ldf", "objective")
  valid_fields = c("fit", summary_stats, "fitted_values", "fit_history",
                 "sampler")

  if (is.vector(output_level) && length(output_level) == 1) {
    output = switch(as.character(output_level),
                    "0" = "fit",
                    "1" = c("fit", summary_stats),
                    "2" = c("fit", summary_stats, "fitted_values"),
                    "3" = c("fit", summary_stats, "fitted_values",
                            "fit_history"),
                    "4" = valid_fields)
  } else if (is.vector(output_level)
             && all(output_level %in% valid_fields)) {
    output = output_level
  } else {
    stop("Invalid output level.")
  }

  return(output)
}


# Handle init_fn parameter
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
  if (!is.function(init_fn) && !exists(init_fn, mode = "function")) {
    stop("The specified init_fn does not exist.")
  }

  return(init_fn)
}


# Handle ebnm_fn parameter
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

  return(list(l = ebnm_fn_l, f = ebnm_fn_f))
}


# Handle ebnm_param parameter
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
handle_ebnm_param = function(ebnm_param,
                             ebnm_fn,
                             n_expected,
                             allow_lists_of_lists = TRUE) {
  if (!is.null(ebnm_param) && !is.list(ebnm_param)) {
    stop(paste("Invalid argument for parameter ebnm_param. A list (or",
               "NULL) was expected."))
  }

  # Check to see whether parameters are specified separately for loadings
  #   and factors:
  if (xor(is.null(ebnm_param[["l"]]), is.null(ebnm_param[["f"]]))) {
    stop(paste("if ebnm_param is specified for either loadings or",
               "factors then it must be specified for both. (Use an",
               "empty list to specify no parameters.)"))
  } else if (!is.null(ebnm_param[["l"]])) {
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
  else if (!allow_lists_of_lists) {
    stop(paste("Different ebnm parameters can only be used for each",
               "loading if at most one type of fit (add_fixed_loadings,",
               "add_fixed_factors, add_greedy, or backfit) is being",
               "performed."))
  }
  else if (length(ebnm_param_l) != n_expected) {
    stop(paste("If different ebnm parameters are used for each loading",
               "then ebnm_param$l must be a list of", n_expected,
               "lists."))
  }

  if (length(ebnm_param_f) == 0) {
    ebnm_param_f = rep(list(list()), n_expected)

  } else if (!is.null(names(ebnm_param_f))) {
    ebnm_param_f = rep(list(ebnm_param_f), n_expected)
  } else if (!allow_lists_of_lists) {
    stop(paste("Different ebnm parameters can only be used for each",
               "factor if at most one type of fit (add_fixed_loadings,",
               "add_fixed_factors, add_greedy, or backfit) is being",
               "performed."))
  } else if (length(ebnm_param_f) != n_expected) {
    stop(paste("If different ebnm parameters are used for each factor",
               "then ebnm_param$f must be a list of", n_expected,
               "lists."))
  }

  # Add defaults:
  ebnm_param_l = lapply(ebnm_param_l,
                        function(x) {add_ebnm_fn_defaults(x, ebnm_fn$l)})
  ebnm_param_f = lapply(ebnm_param_f,
                        function(x) {add_ebnm_fn_defaults(x, ebnm_fn$f)})
  return(list(l = ebnm_param_l, f = ebnm_param_f))
}


# Add defaults to ebnm_param
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
    return(add_ebnm_pn_defaults(ebnm_param))
  }

  return(ebnm_param)
}


# Add ebnm_ash defaults to ebnm_param
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
  if (is.null(ebnm_param$mixcompdist)) {
    ebnm_param$mixcompdist = "normal"
  }
  if (is.null(ebnm_param$method)) {
    ebnm_param$method = "shrink"
  }
  return(ebnm_param)
}


# Add ebnm_pn defaults to ebnm_param
#
# @description Adds default ebnm_pn parameters to ebnm_param.
#
# @param ebnm_param The parameters to be passed into ebnm_fn, given as
#   a single list.
#
# @return ebnm_param, with point-normal defaults added (when not
#   specified by the user).
#
#' @importFrom utils packageVersion
#'
add_ebnm_pn_defaults = function(ebnm_param) {
  if (packageVersion("ebnm") < "0.1.13"
      && (is.null(ebnm_param$warmstart)
          || ebnm_param$warmstart == TRUE)) {
    ebnm_param$warmstart = FALSE
    warning(paste("Setting warmstart = FALSE. Please update ebnm to the",
                  "latest version to use warmstarts."))
  } else if (is.null(ebnm_param$warmstart)) {
    ebnm_param$warmstart = TRUE
  }
  return(ebnm_param)
}


# Handle k parameter
#
# @description Checks that factor k exists.
#
# @param k The factor index.
#
# @param f A fitted flash object.
#
# @return k
#
handle_k = function(k, f) {
  if (flash_get_k(f) < k) {
    stop("Factor k does not exist.")
  }

  return(k)
}


# Handle kset parameter
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
  if (flash_get_k(f) == 0) {
    stop("No factors have been added to the flash object yet.")
  }

  if (is.null(kset)) {
    # Default:
    kset = 1:flash_get_k(f)
  } else if (!is.numeric(kset) || max(kset) > flash_get_k(f)) {
    stop(paste("Invalid kset. Kset should be a vector containing the",
               "indices of the factors to be optimized."))
  }

  return(kset)
}
