flash_add_l_sampler = function(data, f, kset=NULL, ebnm_fn=ebnm_ash) {
  if (is.null(kset)) {kset = 1:get_k(f)}
  f$l_sampler = vector("list", get_k(f))

  # Get a sampler by updating the loadings (this can change the fit if not converged!)
  for (k in kset) {
    f = flash_update_single_loading(data, f, k, ebnm_fn, return_sampler=T)
  }
  return(f)
}

flash_add_f_sampler = function(data, f, kset=NULL, ebnm_fn=ebnm_ash) {
  if (is.null(kset)) {kset = 1:get_k(f)}
  f$f_sampler = vector("list", get_k(f))

  # Get a sampler by updating the loadings (this can change the fit if not converged!)
  for (k in kset) {
    f = flash_update_single_factor(data, f, k, ebnm_fn, return_sampler=T)
  }
  return(f)
}

# Careful not to set nsamp too large!
flash_fixl_samplef = function(data, f, kset=NULL, nsamp) {
  if (is.null(kset)) {kset = 1:get_k(f)}
  if (is.null(f$f_sampler)) {
    stop("Sampler for F needs to be added first!")
  }
  LF = NULL

  # Sampler each factor/loading independently
  for (k in kset) {
    l = f$EL[, k]
    if (is.null(f$f_sampler[[k]])) {
      stop(paste0("Sampler does not exist for factor ", k))
    }
    fsamp = f$f_sampler[[k]](nsamp)
    lfsamp = lapply(split(fsamp, 1:nrow(fsamp)), function(x) {outer(l, x)})
    if (is.null(LF)) { # begin a running total with the first factor/loading
      LF = lfsamp
    } else { # add subsequent factor/loadings to this running total
      LF = mapply(`+`, LF, lfsamp, SIMPLIFY=FALSE)
    }
  }

  return(LF)
}
# warn when we don't get sampler (b/c all values infinite) but proceed anyway
# need to deal with fixed vals

eval_fxn_using_sampler = function(data, f, kset=NULL, fxn, nsamp, batch_size=10) {
  vals = NULL
  n_batch = ceiling(nsamp / batch_size)
  for (i in 1:n_batch) {
    vals = c(vals, sapply(flash_fixl_samplef(data, f, kset, batch_size), fxn))
  }
  names(vals) = NULL
  return(vals)
}
