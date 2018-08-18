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

calc_max_chg = function(new_vals, old_vals) {
  pct_chgs = abs(new_vals/old_vals - 1)
  # Ignore entries where both old and new values are zero.
  max_chg = max(pct_chgs[!is.nan(pct_chgs)], 0)

  return(max_chg)
}
