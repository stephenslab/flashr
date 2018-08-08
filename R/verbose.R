verbose_greedy_next_fl = function(k) {
  message("Fitting factor/loading ", k, "...")
}

verbose_backfit_announce = function() {
  message("Backfitting flash object...")
}

verbose_obj_table_header = function() {
  message("  Iteration          Objective")
}

verbose_diff_table_header = function() {
  message("  Iteration         Difference")
}

verbose_obj_table_entry = function(iteration, obj) {
  message(sprintf("%11d", iteration), sprintf("%19.3f", obj))
}

verbose_diff_table_entry = function(iteration, diff) {
  message(sprintf("%11d", iteration), sprintf("%19.3f", diff))
}

verbose_obj_decrease_warning = function() {
  warning(paste("An iteration decreased the objective.",
                "This happens occasionally, perhaps due to",
                "numeric reasons. You could ignore this",
                "warning, but you might like to check out",
                "https://github.com/stephenslab/flashr/issues/26",
                "for more details."))
}

verbose_nullcheck_announce = function() {
  message("Performing nullcheck...")
}

verbose_nullcheck_delete_fl = function(k, diff) {
  message("  Deleting factor ", k,
          " increases objective by ",
          formatC(diff, format="e", digits=2),
          ". Factor zeroed out.")
}

verbose_nullcheck_keep_fl = function(k, diff) {
  message ("  Deleting factor ", k,
           " decreases objective by ",
           formatC(diff, format="e", digits=2),
           ". Factor retained.")
}

verbose_nullcheck_complete = function(obj) {
  message("  Nullcheck complete. Objective: ",
          round(obj, digits=3))
}
