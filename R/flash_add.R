#' @title add factors or loadings to f
#' @details The precision parameter in f is updated after adding
#' @param data a flash data object
#' @param f a flash fit object
#' @param LL the loadings, an n by K matrix
#' @param FF the factors, a p by K matrix
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Useful for including a mean factor for example.
#' @param fixf a p by K matrix of TRUE/FALSE values; same as fixl but for factors FF.
#' @return a flash fit object, with additional factors initialized using LL and FF
#' @export
flash_add_lf = function(data,LL,FF,f_init=NULL,fixl=NULL,fixf=NULL){
  if(is.null(f_init)){f_init = flash_init_null()}
  f2 = flash_init_lf(LL,FF,fixl,fixf)
  f = flash_combine(f_init,f2)
  return(flash_update_precision(data,f))
}

#' @title add factors to a flash fit object based on data
#' @param data a flash data object
#' @param K number of factors to add
#' @param f_init an existing flash fit object to add to
#' @param init_fn the function to use to initialize new factors (typically some kind of svd-like function)
#' @details Computes the current residuals from data and f_init and adds K new factors based
#' on init_fn applied to these residuals. (If f_init is NULL then the residuals are the data)
#' @export
flash_add_factors_from_data = function(data,K,f_init=NULL,init_fn="udv_si"){
  if(is.null(f_init)){f_init = flash_init_null()}
  R = get_R_withmissing(data,f_init)
  f2 = flash_init_fn(flash_set_data(R),init_fn,K)
  f = flash_combine(f_init,f2)
  return(flash_update_precision(data,f))
}


#' @title  Add a set of fixed loadings to a flash fit object
#' @param data a flash data object
#' @param LL the loadings, an n by K matrix. Missing values will be initialized by the mean of the relevant column (but will generally be
#' re-estimated when refitting the model).
#' @param f_init a flash fit object to which loadings are to be added (if NULL then a new fit object is created)
#' @param fixl an n by K matrix of TRUE/FALSE values indicating which elements of LL should be considered fixed and not changed during updates.
#' Default is to fix all non-missing values, so missing values will be updated when f is updated.
#' @return a flash fit object, with loadings initialized from LL, and corresponding factors initialized to 0.
#' @export
flash_add_fixed_l = function(data, LL, f_init=NULL, fixl=NULL, init_fn="udv_si"){
  if(is.matrix(data)){data = flash_set_data(data)}
  if(is.null(f_init)){f_init = flash_init_null()}
  if(is.null(fixl)){fixl = !is.na(LL)}

  LL_init = LL
  FF_init = matrix(0, nrow=ncol(data$Y), ncol=ncol(LL))
  f = f_init

  f_offset = ncol(f$EL)
  if (is.null(f_offset)) {f_offset = 0}

  # Group columns of LL into blocks, each of which has the same missing data:
  blocks = find_col_blocks(is.na(LL))

  for (i in 1:nrow(blocks)) {
    block_idx = blocks[i, 1]:blocks[i, 2]
    missing = is.na(LL[, block_idx[1]])

    # If we're only missing one element, just replace it with the column mean:
    if (sum(missing) == 1) {
      LL_init[missing, block_idx] = colMeans(LL[!missing], block_idx)
    }
    # If we're missing more, initialize via a subsetted flash object:
    else if (sum(missing) > 1) {
      subf = flash_subset_l(f, missing)
      subdata = flash_subset_data(data, row_subset=missing)

      block_size = length(block_idx)
      subf = flash_add_factors_from_data(subdata, block_size, subf, init_fn)
      LL_init[missing, block_idx] = subf$EL[,f_offset + block_idx]
      FF_init[, block_idx] = subf$EF[,f_offset + block_idx]
    }

    f_new = flash_init_lf(LL_init[,block_idx, drop=F], FF_init[,block_idx, drop=F], fixl=fixl)
    f = flash_combine(f, f_new)
  }

  return(f)
}


#' @title partition matrix into blocks of identical columns
#' @param X a matrix
#' @return a matrix Y with 2 columns. Each row of Y corresponds to a block of columns in X.
#' The first column of Y gives the index of the first column in the block; the second
#' column gives the index of the last column. That is, block i is composed of columns
#' Y[i, 1]:Y[i, 2].
find_col_blocks = function(X) {
  n = nrow(X)
  K = ncol(X)

  if (K == 1) { # just one column, so just one block
    return(matrix(1, nrow=1, ncol=2))
  }

  # Check to see whether column j in X has the same data as column j+1:
  is_col_same = (colSums(X[,1:(K-1),drop=F] == X[,2:K,drop=F]) == n)

  # Group into blocks of columns, each of which have the same data:
  block_ends = which(is_col_same == FALSE)
  start_idx = c(1, block_ends + 1)
  end_idx = c(block_ends, K)
  cbind(start_idx, end_idx)
}


#' @title  Add a set of fixed factors to a flash fit object
#' @param data a flash data object
#' @param FF the factors, a p by K matrix. Missing values will be initialized by the mean of the relevant column (but will generally be
#' re-estimated when refitting the model).
#' @param f_init a flash fit object to which factors are to be added (if NULL then a new fit object is created)
#' @param fixf a p by K matrix of TRUE/FALSE values indicating which elements of FF should be considered fixed and not changed during updates.
#' Default is to fix all non-missing values, so missing values will be updated when f is updated.
#' @return a flash fit object, with factors initialized from FF, and corresponding loadings initialized to 0.
#' @export
flash_add_fixed_f = function(data,FF,f_init=NULL,fixf=NULL){
  if (is.matrix(data)) {data = flash_set_data(data)}

  tf = flash_add_fixed_l(flash_transpose_data(data), FF, flash_transpose(f_init), fixf)
  return(flash_transpose(tf))
}


#' @title  Initialize a flash fit object by applying a function to data
#' @param data a flash data object
#' @param init_fn an initialization function, which takes as input an (n by p matrix, or flash data object)
#' and K, a number of factors, and and outputs a list with elements (u,d,v)
#' @return a flash fit object
flash_init_fn = function(data,init_fn,K=1){
  s = do.call(init_fn,list(get_Yorig(data),K))
  f = flash_init_udv(s,K)
  f = flash_update_precision(data,f)
  return(f)
}
