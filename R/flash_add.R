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
flash_add_lf = function(data, LL, FF, f_init = NULL,
                        fixl = NULL, fixf = NULL) {
  data = handle_data(data)
  f_init = handle_f(f_init, init_null_f = TRUE)

  f2 = flash_init_lf(LL, FF, fixl, fixf)
  f = flash_combine(f_init, f2)

  return(f)
}


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
#' @export
#'
flash_add_factors_from_data = function(data, K, f_init = NULL,
                                       init_fn = "udv_si") {
  data = handle_data(data)
  f_init = handle_f(f_init, init_null_f = TRUE)
  init_fn = handle_init_fn(init_fn)

  R  = flash_get_R_withmissing(data, f_init)
  f2 = flash_init_fn(flash_set_data(R), init_fn, K)
  f  = flash_combine(f_init, f2)

  return(f)
}


#' @title Add a set of fixed loadings to a flash fit object.
#'
#' @inheritParams flash
#'
#' @param LL The loadings, an n by K matrix. Missing values will be
#'   initialized by the mean of the relevant column (but will generally
#'   be re-estimated when refitting the model).
#'
#' @param fixl An n by K matrix of \code{TRUE}/\code{FALSE} values
#'   indicating which elements of \code{LL} should be considered fixed
#'   and not changed during updates.  The default is to fix all
#'   non-missing values, so missing values will be updated when the
#'   flash object is updated.
#'
#' @return A flash fit object, with loadings initialized from \code{LL},
#'   and corresponding factors initialized to zero.
#'
#' @export
#'
flash_add_fixed_l = function(data, LL, f_init = NULL, fixl = NULL,
                             init_fn = "udv_si") {
  data = handle_data(data)
  f_init = handle_f(f_init, init_null_f = TRUE)

  if (is.null(fixl)) {
    fixl = !is.na(LL)
  }

  LL_init = LL
  FF_init = matrix(0, nrow=ncol(data$Y), ncol=ncol(LL))

  k_offset = ncol(f_init$EL)
  if (is.null(k_offset)) {k_offset = 0}

  # Group columns of LL into blocks, each of which has the same
  # missing data.
  blocks = find_col_blocks(is.na(LL))

  f = f_init
  for (i in 1:length(blocks)) {
    block_cols = blocks[[i]]
    missing_rows = is.na(LL[, block_cols[1]])

    # If we're only missing one element, just replace it with the
    # column mean.
    if (sum(missing_rows) == 1) {
      LL_init[missing_rows, block_cols] =
        colMeans(LL[!missing_rows, block_cols, drop=F])
  } else if (sum(missing_rows) > 1) {

      # If we're missing more, initialize via a subsetted flash object.
      subf = flash_subset_l(f, missing_rows)
      subdata = flash_subset_data(data, row_subset=missing_rows)
      subf = flash_add_factors_from_data(subdata, length(block_cols), subf,
                                         init_fn)
      LL_init[missing_rows, block_cols] = subf$EL[,k_offset + block_cols]
      FF_init[, block_cols] = subf$EF[,k_offset + block_cols]
    }

    f = flash_add_lf(data, LL_init[,block_cols, drop=F],
                     FF_init[,block_cols, drop=F],
                     f, fixl=fixl[,block_cols, drop=F])
  }

  return(f)
}


# @title Partition a matrix into blocks of identical columns.
#
# @param X the matrix to be partitioned (note that X should not have NAs).
#
# @return A list, each element of which contains the indices of a
#   single block of identical columns.
#
find_col_blocks = function(X) {
  n = nrow(X)
  K = ncol(X)

  if (K == 1) { # just one column, so just one block
    return(as.list(1))
  }

  # Check to see whether column j in X has the same data as column j+1.
  is_col_same = (colSums(X[,1:(K-1),drop=F] == X[,2:K,drop=F]) == n)

  # Group into blocks of columns; all columns in a single block have
  # the same data.
  block_ends = which(is_col_same == FALSE)
  start_idx = c(1, block_ends + 1)
  end_idx = c(block_ends, K)
  blocks = vector("list", length(start_idx))
  for (i in 1:length(start_idx)) {
    blocks[[i]] = start_idx[i]:end_idx[i]
  }

  return(blocks)
}


#' @title Add a set of fixed factors to a flash fit object.
#'
#' @inheritParams flash
#'
#' @param FF The factors, a p by K matrix. Missing values will be
#'   initialized by the mean of the relevant column (but will generally
#'   be re-estimated when refitting the model).
#'
#' @param fixf A p by K matrix of of \code{TRUE}/\code{FALSE} values
#'   indicating which elements of \code{FF} should be considered fixed
#'   and not changed during updates.  The default is to fix all
#'   non-missing values, so missing values will be updated when the
#'   flash object is updated.
#'
#' @return A flash fit object, with factors initialized from \code{FF},
#'   and corresponding loadings initialized to zero.
#'
#' @export
#'
flash_add_fixed_f = function(data, FF, f_init=NULL, fixf=NULL) {
  data = handle_data(data)
  # f_init is handled by flash_add_fixed_l

  tf = flash_add_fixed_l(flash_transpose_data(data),
                         FF,
                         flash_transpose(f_init),
                         fixf)

  return(flash_transpose(tf))
}
