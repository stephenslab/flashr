#' @title Add factors or loadings to f.
#'
#' @param data A flash data object.
#'
#' @param LL The loadings, an n by K matrix.
#'
#' @param FF The factors, a p by K matrix.
#'
#' @param f_init Description of input argument goes here.
#'
#' @param fixl An n by K matrix of TRUE/FALSE values indicating which
#'   elements of LL should be considered fixed and not changed during
#'   updates. Useful for including a mean factor for example.
#'
#' @param fixf A p by K matrix of TRUE/FALSE values; same as fixl but
#'   for factors FF.
#'
#' @return A flash fit object, with additional factors initialized
#'   using LL and FF.
#'
#' @export
#'
flash_add_lf = function(data, LL, FF, f_init=NULL, fixl=NULL, fixf=NULL) {
  if (is.null(f_init)){
    f_init = flash_init_null()
  }

  f2 = flash_init_lf(LL, FF, fixl, fixf)
  f = flash_combine(f_init, f2)

  return(f)
}

#' @title Add factors to a flash fit object based on data.
#'
#' @description Computes the current residuals from data and f_init
#'   and adds K new factors based on init_fn applied to these
#'   residuals. (If f_init is NULL then the residuals are the data).
#'
#' @param data A flash data object.
#'
#' @param K Number of factors to add.
#'
#' @param f_init An existing flash fit object to add to.
#'
#' @param init_fn The function to use to initialize new factors
#'   (typically some kind of svd-like function).
#'
#' @export
#'
flash_add_factors_from_data = function(data, K, f_init=NULL,
                                       init_fn="udv_si") {
  if (is.null(f_init)) {
    f_init = flash_init_null()
  }

  R  = flash_get_R_withmissing(data, f_init)
  f2 = flash_init_fn(flash_set_data(R), init_fn, K)
  f  = flash_combine(f_init, f2)

  return(f)
}

#' @title Add a set of fixed loadings to a flash fit object.
#'
#' @param data A flash data object.
#'
#' @param LL The loadings, an n by K matrix. Missing values will be
#'   initialized by the mean of the relevant column (but will generally
#'   be re-estimated when refitting the model).
#'
#' @param f_init A flash fit object to which loadings are to be added
#'   (if NULL then a new fit object is created).
#'
#' @param fixl An n by K matrix of TRUE/FALSE values indicating which
#'   elements of LL should be considered fixed and not changed during
#'   updates.  Default is to fix all non-missing values, so missing
#'   values will be updated when f is updated.
#'
#' @param init_fn Description of input argument goes here.
#'
#' @return A flash fit object, with loadings initialized from LL, and
#'   corresponding factors initialized to 0.
#'
#' @export
#'
flash_add_fixed_l = function(data, LL, f_init=NULL, fixl=NULL,
  init_fn="udv_si") {
  if(is.matrix(data)){data = flash_set_data(data)}
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixl)){fixl = !is.na(LL)}

  LL_init = LL
  FF_init = matrix(0, nrow=ncol(data$Y), ncol=ncol(LL))

  k_offset = ncol(f_init$EL)
  if (is.null(k_offset)) {k_offset = 0}

  # Group columns of LL into blocks, each of which has the same missing data:
  blocks = find_col_blocks(is.na(LL))

  f = f_init
  for (i in 1:length(blocks)) {
    block_cols = blocks[[i]]
    missing_rows = is.na(LL[, block_cols[1]])

    # If we're only missing one element, just replace it with the column mean:
    if (sum(missing_rows) == 1) {
      LL_init[missing_rows, block_cols] = colMeans(LL[!missing_rows, block_cols, drop=F])
    } else if (sum(missing_rows) > 1) { # If we're missing more, initialize via a subsetted flash object:
      subf = flash_subset_l(f, missing_rows)
      subdata = flash_subset_data(data, row_subset=missing_rows)
      subf = flash_add_factors_from_data(subdata, length(block_cols), subf, init_fn)
      LL_init[missing_rows, block_cols] = subf$EL[,k_offset + block_cols]
      FF_init[, block_cols] = subf$EF[,k_offset + block_cols]
    }

    f = flash_add_lf(data, LL_init[,block_cols, drop=F], FF_init[,block_cols, drop=F],
                     f, fixl=fixl[,block_cols, drop=F])
  }

  return(f)
}

# @title Partition a matrix into blocks of identical columns
#
# @param X the matrix to be partitioned (note that X should not have NAs)
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

  # Check to see whether column j in X has the same data as column j+1:
  is_col_same = (colSums(X[,1:(K-1),drop=F] == X[,2:K,drop=F]) == n)

  # Group into blocks of columns; all columns in a single block have
  # the same data:
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
#' @param data a flash data object
#'
#' @param FF The factors, a p by K matrix. Missing values will be
#'   initialized by the mean of the relevant column (but will generally
#'   be re-estimated when refitting the model).
#'
#' @param f_init A flash fit object to which factors are to be added
#'   (if NULL then a new fit object is created).
#'
#' @param fixf A p by K matrix of TRUE/FALSE values indicating which
#'   elements of FF should be considered fixed and not changed during
#'   updates. Default is to fix all non-missing values, so missing
#'   values will be updated when f is updated.
#'
#' @return A flash fit object, with factors initialized from FF, and
#'   corresponding loadings initialized to 0..
#'
#' @export
#'
flash_add_fixed_f = function(data, FF, f_init=NULL, fixf=NULL) {
  if (is.matrix(data)) {data = flash_set_data(data)}

  tf = flash_add_fixed_l(flash_transpose_data(data), FF, flash_transpose(f_init), fixf)
  return(flash_transpose(tf))
}

# @title Initialize a flash fit object by applying a function to data.
#
# @param data a flash data object.
#
# @param init_fn An initialization function, which takes as input an
#   (n by p matrix, or flash data object) and K, a number of factors,
#   and and outputs a list with elements (u,d,v).
#
# @return A flash fit object.
#
flash_init_fn = function(data, init_fn, K=1) {
  s = do.call(init_fn, list(get_Yorig(data), K))
  f = flash_init_udv(s, K)
  return(f)
}
