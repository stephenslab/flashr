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
  if (is.matrix(EL)) {
    lnorms = sqrt(apply(EL^2, 2, sum))
    fnorms = sqrt(apply(EF^2, 2, sum))
    EL = as.vector(sweep(EL, 2, lnorms, `/`))
    EF = as.vector(sweep(EF, 2, fnorms, `/`))
  } else {
    EL = EL / sqrt(sum(EL^2))
    EF = EF / sqrt(sum(EF^2))
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
