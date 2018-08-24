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
    verbose_nullcheck_announce()
  }

  zeroed_out = integer(0)

  f_changed = TRUE  # We are going to iterate until f does not change.
  while (f_changed) {

    f_changed = FALSE
    for (k in kset) {
      f0 = flash_zero_out_factor(f, k)
      if (!identical(f, f0)) {
        f0 = flash_update_precision(data, f0, var_type)
        obj0 = flash_get_objective(data, f0)
        obj1 = flash_get_objective(data, f)

        if (obj0 >= obj1) {
          if (verbose) {
            verbose_nullcheck_delete_fl(k, obj0 - obj1)
          }
          f = f0
          zeroed_out = c(zeroed_out, k)
          f_changed = TRUE
        } else if (obj1 > obj0) {
          if (verbose) {
            verbose_nullcheck_keep_fl(k, obj1 - obj0)
          }
        }
      }
    }
  }

  if (verbose) {
    verbose_nullcheck_complete(flash_get_objective(data, f))
  }

  return(list(f=f, zeroed_out=zeroed_out))
}
