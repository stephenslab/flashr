# Messages displayed when verbose = TRUE.

verbose_greedy_next_fl = function(k, stopping_rule, tol) {
  message("Fitting factor/loading ", k, " (",
          stopping_criterion_string(stopping_rule, tol), "):")
}

verbose_backfit_announce = function(n, stopping_rule, tol) {
  message("Backfitting ", n, " factor/loading(s) (",
          stopping_criterion_string(stopping_rule, tol), "):")
}

stopping_criterion_string = function(stopping_rule, tol) {
  rule_string = ifelse(stopping_rule == "objective",
                       "difference in objective",
                       "maximum parameter change")
  tol_string = ifelse(stopping_rule == "objective",
                      formatC(tol, format = "e", digits = 2),
                      paste0(100 * tol, "%"))
  return(paste("stop when", rule_string, "is <", tol_string))
}

verbose_obj_table_header = function(stopping_rule,
                                    track_obj,
                                    track_param_chg) {
  header_string = "  Iteration"
  if (track_obj) {
    header_string = paste0(header_string,
                           sprintf("%17s", "Objective"))
  }
  if (stopping_rule == "objective") {
    header_string = paste0(header_string,
                           sprintf("%17s", "Obj Diff"))
  }
  if (track_param_chg != "none") {
    header_string = paste0(header_string,
                           sprintf("%17s", "Max Param Chg"))
  }
  message(header_string)
}

verbose_obj_table_entry = function(iteration,
                                   obj,
                                   obj_diff,
                                   max_chg,
                                   stopping_rule) {
  entry_string = sprintf("%11d", iteration)
  if (!is.null(obj)) {
    entry_string = paste0(entry_string,
                          sprintf("%17.2f", obj))
  }
  if (stopping_rule == "objective") {
    diff_string = formatC(obj_diff, format="e", digits=2)
    entry_string = paste0(entry_string,
                          sprintf("%17s", diff_string))
  }
  if (!is.null(max_chg)) {
    entry_string = paste0(entry_string,
                          sprintf("%16.2f", 100 * max_chg), "%")
  }

  message(entry_string)
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
          round(obj, digits=2))
}
