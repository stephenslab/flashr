#' @title Generates LF sampler
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
flash_lf_sampler = function(data,
                            f,
                            kset = NULL,
                            fixed = c("factors", "loadings")) {
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
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
# @inheritParams flash_lf_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of LF to be produced by the sampler.
#
flash_lf_sampler_fixedf = function(data, f, kset) {
  l_sampler = flash_l_sampler(data, f, kset)

  function(nsamp) {
    lsamp = l_sampler(nsamp)
    return(mapply(function(L) {L %*% t(f$EF[, kset])}, lsamp, SIMPLIFY=FALSE))
  }
}


# @title Generates LF sampler with L fixed at its expectation
#
# @description Generates function that samples LF from a flash fit
#   object, with L fixed at its posterior mean and the columns of F
#   sampled independently from their marginal posteriors.
#
# @inheritParams flash_lf_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of LF to be produced by the sampler.
#
flash_lf_sampler_fixedl = function(data, f, kset=NULL) {
  f_sampler = flash_f_sampler(data, f, kset)

  function(nsamp) {
    fsamp = f_sampler(nsamp)
    return(mapply(function(F) {f$EL[, kset] %*% t(F)}, fsamp, SIMPLIFY=FALSE))
  }
}


# @title Generates sampler for L
#
# @description Generates function that samples L from a flash fit object.
#   The columns of L are sampled independently from their marginal posteriors
#   (conditional on F being fixed at its expectation).
#
# @inheritParams flash_lf_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of L to be produced by the sampler. This sampler returns a list
#   of matrices.
#
flash_l_sampler = function(data, f, kset=NULL) {
  if (is.matrix(data)) {
    data = flash_set_data(data)
  }
  kset = handle_kset(kset, f)

  sampler_list = vector("list", flash_get_k(f))
  for (k in kset) {
    # Use ebnm function and parameters from flash object
    ebnm_fn = f$ebnm_fn[[k]]
    ebnm_param = f$ebnm_param_l[[k]]
    sampler_list[[k]] = flash_update_single_loading(data, f, k,
                                                    ebnm_fn,
                                                    ebnm_param,
                                                    return_sampler=T)
  }

  function(nsamp) {
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
  }
}


# @title Generates sampler for F
#
# @description Generates function that samples F from a flash fit object.
#   The columns of F are sampled independently from their marginal
#   posteriors (conditional on L being fixed at its expectation).
#
# @inheritParams flash_lf_sampler
#
# @return A function that takes a single parameter nsamp, the number of
#   samples of F to be produced by the sampler. This sampler returns a list
#   of matrices.
#
flash_f_sampler = function(data, f, kset=NULL) {
  if (is.matrix(data)) {data = flash_set_data(data)}
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
  function(nsamp) {
    samp = matrix(0, nrow=nsamp, ncol=length(is_fixed))
    samp[, is_fixed] = fixed_vals
    if (!is.null(sample_fxn)) { # sample_fxn is NULL when all values are fixed
      samp[, !is_fixed] = sample_fxn(nsamp)
    }
    return(samp)
  }
}
