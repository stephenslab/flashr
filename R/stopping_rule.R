# Determines whether convergence has occurred when greedily adding a factor/
#   loading or when backfitting.
#
is_converged = function(stopping_rule, tol, obj_diff, max_chg_l, max_chg_f) {
  if (stopping_rule == "objective") {
    if (obj_diff < 0) {
      verbose_obj_decrease_warning()
    }
    return(obj_diff < tol)
  } else if (stopping_rule == "loadings") {
    return(max_chg_l < tol)
  } else if (stopping_rule == "factors") {
    return(max_chg_f < tol)
  } else { # stopping_rule == "all_params"
    return(max(max_chg_l, max_chg_f) < tol)
  }
}

# Normalizes EL and EF before changes in parameter values are calculated.
#
normalize_lf = function(EL, EF) {
  # normalize both factors and loadings to have l_1 norm equal to 1:
  if (is.matrix(EL)) {
    lnorms = apply(abs(EL), 2, max)
    EL = as.vector(sweep(EL, 2, lnorms, `/`))
    fnorms = apply(abs(EF), 2, max)
    EF = as.vector(sweep(EF, 2, fnorms, `/`))
  } else {
    EL = EL / max(abs(EL))
    EF = EF / max(abs(EF))
  }

  # # normalize loadings to have l_2 norm equal to 1:
  # if (is.matrix(EL)) {
  #   norms = sqrt(colSums(EL^2))
  #   EL = as.vector(sweep(EL, 2, norms, `/`))
  #   EF = as.vector(sweep(EF, 2, norms, `*`))
  # } else {
  #   norm = sqrt(sum(EL^2))
  #   EL = EL / norm
  #   EF = EF * norm
  # }

  return(list(EL = EL, EF = EF))
}

calc_max_chg = function(new_vals, old_vals) {
  # Absolute difference:
  chgs = abs(new_vals - old_vals)
  max_chg = max(chgs[!is.nan(chgs)], 0)

  # # Relative difference:
  # pct_chgs = abs(new_vals/old_vals - 1)
  # max_chg = max(pct_chgs[!is.nan(pct_chgs)], 0)

  return(max_chg)
}
