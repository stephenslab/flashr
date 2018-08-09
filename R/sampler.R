#' @title Generate posterior sampling function
#'
#' @description Generates function that samples LF from a flash fit
#'   object, with either L or F fixed at its posterior mean and the
#'   columns of F or L sampled independently from their marginal
#'   posteriors.
#'
#' @param data An n by p matrix or a flash data object created using
#'   \code{flash_set_data}.
#'
#' @param f A fitted flash object.
#'
#' @param kset The indices of factors to be optimized (NULL indicates
#'   all factors).
#'
#' @param fixed Indicates whether to fix factors or loadings at their
#'   posterior mean.
#'
#' @return A function that takes a single parameter \code{nsamp}, the
#'   number of samples of LF to be produced by the sampler. Care should
#'   be used when setting \code{nsamp}, because the sampler returns a
#'   list of matrices which are each of the same size as the data matrix.
#'
#' @export
#'
flash_sampler = function(data,
                         f,
                         kset = NULL,
                         fixed = c("factors", "loadings")) {
  data = handle_data(data)
  f = handle_f(f, allow_null = FALSE)
  kset = handle_kset(kset, f)

  fixed = match.arg(fixed)
  if (fixed == "factors") {
    return(flash_lf_sampler_fixedf(data, f, kset))
  } else if (fixed == "loadings") {
    return(flash_lf_sampler_fixedl(data, f, kset))
  }
}


# @title Generates LF sampler with F fixed at its expectation
#
# @description Generates function that samples LF from a flash fit
#   object, with F fixed at its posterior mean and the columns of L
#   sampled independently from their marginal posteriors.
#
# @inheritParams flash_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of LF to be produced by the sampler.
#
flash_lf_sampler_fixedf = function(data, f, kset) {
  l_sampler = flash_l_sampler(data, f, kset)

  return(function(nsamp) {
    lsamp = l_sampler(nsamp)
    return(mapply(function(L) {L %*% t(f$EF[, kset])},
                  lsamp,
                  SIMPLIFY=FALSE))
  })
}


# @title Generates LF sampler with L fixed at its expectation
#
# @description Generates function that samples LF from a flash fit
#   object, with L fixed at its posterior mean and the columns of F
#   sampled independently from their marginal posteriors.
#
# @inheritParams flash_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of LF to be produced by the sampler.
#
flash_lf_sampler_fixedl = function(data, f, kset) {
  f_sampler = flash_f_sampler(data, f, kset)

  return(function(nsamp) {
    fsamp = f_sampler(nsamp)
    return(mapply(function(F) {f$EL[, kset] %*% t(F)},
                  fsamp,
                  SIMPLIFY=FALSE))
  })
}


# @title Generates sampler for L
#
# @description Generates function that samples L from a flash fit object.
#   The columns of L are sampled independently from their marginal
#   posteriors (conditional on F being fixed at its expectation).
#
# @inheritParams flash_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of L to be produced by the sampler. This sampler returns a
#   list of matrices.
#
#' @importFrom utils packageVersion
#'
flash_l_sampler = function(data, f, kset) {
  sampler_list = vector("list", flash_get_k(f))

  for (k in kset) {
    # Use ebnm function and parameters from flash object
    ebnm_fn = f$ebnm_fn_l[[k]]
    if (is.null(ebnm_fn)) {
      stop(paste("Factor/loading", k, "has not yet been fit."))
    }
    ebnm_param = f$ebnm_param_l[[k]]
    ebnm_param = add_l_sampler_params(ebnm_param, ebnm_fn, f, k)
    sampler_list[[k]] = flash_single_l_sampler(data, f, k, ebnm_fn,
                                               ebnm_param)
  }

  return(function(nsamp) {
    L = NULL
    for (k in kset) {
      lsamp = t(sampler_list[[k]](nsamp)) # transpose samples to columns...
      lsamp = split(lsamp, col(lsamp)) # ...and split into list
      if (is.null(L)) { # the first factor/loading
        L = lsamp
      } else { # subsequent factor/loadings
        L = mapply(cbind, L, lsamp, SIMPLIFY=FALSE)
      }
    }
    return(L)
  })
}

add_l_sampler_params = function(ebnm_param, ebnm_fn, f, k) {
  # Does not work for early versions of ebnm_pn
  if (ebnm_fn == "ebnm_ash" || (ebnm_fn == "ebnm_pn" &&
                                packageVersion("ebnm") > "0.1.11")) {
    ebnm_param$g = f$gl[[k]]
    ebnm_param$fixg = TRUE
  }

  return(ebnm_param)
}


# @title Generates sampler for a single loading
#
# @inheritParams flash_update_single_loading
#
# @param k The index of the loading to sample.
#
# @details The sampler is produced by a call to ebnm_fn with argument
#   output = "post_sampler".
#
flash_single_l_sampler = function(data, f, k, ebnm_fn, ebnm_param) {
  subset = which(!f$fixl[, k])
  if (length(subset) == 0 || all(f$EL2[subset, k] == 0)) {
    # All values are fixed or all non-fixed values are zero:
    return(sampler(rep(TRUE, length(f$EL[, k])), NULL, f$EL[, k]))
  }

  ebnm_args = calc_ebnm_l_args(data, f, k, subset)
  if (is.null(ebnm_args)) {
    stop(paste("All standard errors for either factor or loading", k,
               "are infinite. Impossible to create sampler."))
  }

  post_sampler = do.call(ebnm_fn, list(ebnm_args$x, ebnm_args$s,
                                       ebnm_param,
                                       output = "post_sampler"))

  if (class(post_sampler) != "function") {
    stop("No sampler implemented for that ebnm function.")
  }

  return(sampler(f$fixl[, k], post_sampler, f$EL[f$fixl[, k], k]))
}


# @title Generates sampler for F
#
# @description Generates function that samples F from a flash fit object.
#   The columns of F are sampled independently from their marginal
#   posteriors (conditional on L being fixed at its expectation).
#
# @inheritParams flash_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of F to be produced by the sampler. This sampler returns a
#   list of matrices.
#
flash_f_sampler = function(data, f, kset=NULL) {
  return(flash_l_sampler(flash_transpose_data(data),
                         flash_transpose(f),
                         kset))
}


# @title Generates sampler for a single factor/loading
#
# @param is_fixed A vector of Booleans that indicates which elements are
#   fixed.
#
# @param sample_fxn A function that returns samples for non-fixed
#   elements. The function should return a matrix whose rows correspond
#   to individual samples.
#
# @param fixed_vals The values of the fixed elements.
#
# @return A function that generates samples for the factor/loading. This
#   function returns a matrix whose rows correspond to individual samples.
#
sampler = function(is_fixed, sample_fxn, fixed_vals) {
  return(function(nsamp) {
    samp = matrix(0, nrow=nsamp, ncol=length(is_fixed))
    samp[, is_fixed] = fixed_vals
    if (!is.null(sample_fxn)) {
      # sample_fxn is NULL when all values are fixed
      samp[, !is_fixed] = sample_fxn(nsamp)
    }
    return(samp)
  })
}
