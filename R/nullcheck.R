# @title Zeros out factors when that improves the objective.
#
# @description Sometimes zeroing out a factor can improve the
#   objective. This function iterates over factors with indices in
#   kset and checks whether zeroing it out will improve the objective;
#   if so then that factor is set to 0 (and precision is updated).
#   Returns the final flash fit object obtained when this iterative
#   process stops (ie a complete pass is performed with no factor being
#   zeroed).
#
# @param data A flash data object.
#
# @param f A flash object.
#
# @param kset The indices of the factor/loading to check.
#
# @param var_type Type of variance structure to assume for residuals.
#
# @param verbose If TRUE, various output progress updates will be
#   printed.
#
# @return A flash object.
#
perform_nullcheck = function(data, f, kset, var_type, verbose) {
  if (verbose) {
    message("Performing nullcheck...")
  }

  f_changed = TRUE  # We are going to iterate until f does not change.
  while (f_changed) {

    f_changed = FALSE
    for (k in kset) {

      f0 = flash_zero_out_factor(f, k)
      f0 = flash_update_precision(data, f0, var_type)
      F0 = flash_get_objective(data, f0)
      F1 = flash_get_objective(data, f)

      if (F0 > F1) {
        if (verbose) {
          message("  Deleting factor ", k,
                  " increases objective by ",
                  formatC(F0 - F1, format="e", digits=2),
                  ". Factor zeroed out.")
        }
        f = f0
        f_changed = TRUE
      } else if (F1 > F0) {
        if (verbose) {
          message ("  Deleting factor ", k,
                   " decreases objective by ",
                   formatC(F1 - F0, format="e", digits=2),
                   ". Factor retained.")
        }
      }

    }
  }
  if (verbose) {
    message("  Nullcheck complete. Objective: ",
            round(flash_get_objective(data, f), digits=3))
  }
  return(f)
}
