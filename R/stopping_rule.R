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

# Normalizes EL and EF before changes in parameter values are calculated. In
#   order for `tol` to have a similar meaning regardless of the size of the
#   data set, we normalize each loading (i.e., each column of EL) so that the
#   mean absolute value of the elements of the loading is equal to one (and
#   similarly for the factors).
#
normalize_lf = function(EL, EF) {
  if (is.matrix(EL)) {
    lnorms = apply(abs(EL), 2, mean)
    EL = as.vector(sweep(EL, 2, lnorms, `/`))
    fnorms = apply(abs(EF), 2, mean)
    EF = as.vector(sweep(EF, 2, fnorms, `/`))
  } else {
    EL = EL / mean(abs(EL))
    EF = EF / mean(abs(EF))
  }

  return(list(EL = EL, EF = EF))
}

# Calculates the maximum change in parameter values. Since EL and EF are
#   normalized, we use absolute changes rather than relative ones.
#
calc_max_chg = function(new_vals, old_vals) {
  # Absolute difference:
  chgs = abs(new_vals - old_vals)
  max_chg = max(chgs[!is.nan(chgs)], 0)

  # # Relative difference:
  # pct_chgs = abs(new_vals/old_vals - 1)
  # max_chg = max(pct_chgs[!is.nan(pct_chgs)], 0)

  return(max_chg)
}
