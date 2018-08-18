is_converged = function(stopping_rule, tol, obj_diff, max_chg) {
  if (stopping_rule == "objective") {
    if (obj_diff < 0) {
      verbose_obj_decrease_warning()
    }
    return(obj_diff < tol)
  } else {
    return(max_chg < tol)
  }
}

calc_max_chg = function(EL, EF, old_EL, old_EF, track_param_chg) {
  old_vals = new_vals = numeric(0)
  if (track_param_chg != "factors") {
    old_vals = old_EL
    new_vals = EL
  }
  if (track_param_chg != "loadings") {
    old_vals = c(old_vals, old_EF)
    new_vals = c(new_vals, EF)
  }

  pct_chgs = abs(new_vals/old_vals - 1)
  # Ignore entries where both old and new values are zero.
  max_chg = max(pct_chgs[!is.nan(pct_chgs)], 0)

  return(max_chg)
}
