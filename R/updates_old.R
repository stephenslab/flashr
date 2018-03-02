# Some old code I want to keep hanging around for now, e.g., for
# testing and comparisons.

# @title Optimize a flash factor-loading combination.
# 
# @description Iteratively updates factor and loading k of f (as well
#   as residual precision) to convergence of objective (used in the
#   greedy algorithm for example).
# 
# @param data A flash data object.
# 
# @param f A flash object.
# 
# @param k The index of the factor/loading to optimize.
# 
# @param var_type Type of variance structure to assume for residuals.
# 
# @param nullcheck Flag whether to check, after running hill-climbing
#   updates, whether the achieved optimum is better than setting factor
#   to 0. If this check is performed and fails then the factor will be
#   set to 0 in the returned fit.
# 
# @param tol A tolerance for the optimization.
# 
# @param ebnm_fn Function to solve the Empirical Bayes normal means
#   problem.
# 
# @param ebnm_param Parameters to be passed to ebnm_fn when
#   optimizing.
# 
# @param verbose If TRUE, various output progress updates will be
#   printed.
# 
# @return An updated flash object.
# 
flash_optimize_single_fl_old <-
  function(data, f, k, var_type, nullcheck = TRUE, tol = 0.01,
           ebnm_fn = ebnm_ash, ebnm_param = flash_default_ebnm_param(ebnm_fn),
    verbose = FALSE) {
    f = flash_update_single_fl(data, f, k, var_type, ebnm_fn, ebnm_param)  #do an update first so that we can get a valid c
    c = flash_get_conv_criteria(data, f)
    if (verbose) {
        message("objective: ", c)
    }
    diff = 1
    while (diff > tol) {
        f = flash_update_single_fl(data, f, k, var_type, ebnm_fn, ebnm_param)
        cnew = flash_get_conv_criteria(data, f)
        diff = cnew - c
        c = cnew
        if (verbose) {
            message("objective: ", c)
        }
    }

    if (nullcheck) {
        f = perform_nullcheck(data, f, k, var_type, verbose)
    }

    return(f)
}

# @title Fit the rank1 flash model to data.
# 
# @param data An n by p matrix or a flash data object created using
#   \code{flash_set_data}.
# 
# @param f_init If supplied, a flash object to which a single new
#   factor is to be added.
# 
# @param var_type Type of variance structure to assume for residuals.
# 
# @param tol Specify how much objective can change in a single
#   iteration to be considered not converged.
# 
# @param init_fn function to be used to initialize the factor. This
#   function should take parameters (Y,K) where Y is an n by p matrix
#   of data (or a flash data object) and K is a number of factors.  It
#   should output a list with elements (u,d,v) where u is n by K matrix
#   v is a p by K matrix and d is a K vector. See \code{udv_si} for an
#   example. (If the input data includes missing values then this
#   function must be able to deal with missing values in its input
#   matrix.)
# 
# @param ebnm_fn Function to solve the Empirical Bayes Normal Means
#   problem.
# 
# @param ebnm_param Parameters to be passed to ebnm_fn when
# optimizing; defaults set by \code{flash_default_ebnm_param()}.
# 
# @param verbose If TRUE, various output progress updates will be
#   printed.
# 
# @param nullcheck Flag whether to check, after running hill-climbing
#   updates, whether the achieved optimum is better than setting factor
#   to 0. If this check is performed and fails then the factor will be
#   set to 0 in the returned fit.
# 
# @return A fitted flash object.
# 
# @examples
# 
# Y = matrix(rnorm(100),nrow=5,ncol=20)
# f = flash_r1(Y)
#
# # Run with the faster ebnm function (uses point-normal prior).
# f2 = flash_r1(ebnm_fn = ebnm_pn)
# 
flash_r1_old = function(data, f_init = NULL,
                        var_type = c("by_column", "constant"),
                        init_fn = "udv_si", tol = 0.01,
                        ebnm_fn = ebnm_ash,
                        ebnm_param = flash_default_ebnm_param(ebnm_fn),
                        verbose = FALSE, nullcheck = TRUE) {
    if (is.matrix(data)) {
        data = flash_set_data(data)
    }
    var_type = match.arg(var_type)
    f = flash_add_factors_from_data(data, f_init = f_init,
                                    init_fn = init_fn, K = 1)
    f = flash_optimize_single_fl_old(data, f, flash_get_k(f), var_type,
                                     nullcheck, tol, ebnm_fn, ebnm_param,
                                     verbose)
    return(f)
}

