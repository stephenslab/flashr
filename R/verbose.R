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
                       "difference in obj.",
                       "max. parameter change")
  tol_string = ifelse(stopping_rule == "objective",
                      formatC(tol, format = "e", digits = 2),
                      paste0(100 * tol, "%"))
  return(paste("stop when", rule_string, "is <", tol_string))
}

verbose_obj_table_header = function(stopping_rule,
                                    track_obj,
                                    track_param_chg) {
  header_string = "  Iteration"

  # if (ebnm_fn_l == "ebnm_pn") {
  #   header_string = paste0(header_string,
  #                          sprintf("%10s", "pi0 (l)"))
  # }
  # if (ebnm_fn_f == "ebnm_pn") {
  #   header_string = paste0(header_string,
  #                          sprintf("%10s", "pi0 (f)"))
  # }

  if (track_param_chg != "none") {
    header_string = paste0(header_string,
                           sprintf("%16s", "Max Param Chg"))
  }

  if (track_obj) {
    header_string = paste0(header_string,
                           sprintf("%16s", "Objective"))
  }
  if (stopping_rule == "objective") {
    header_string = paste0(header_string,
                           sprintf("%11s", "Obj Diff"))
  }

  message(header_string)
}

verbose_obj_table_entry = function(iteration,
                                   obj,
                                   obj_diff,
                                   max_chg,
                                   stopping_rule) {
  entry_string = sprintf("%11d", iteration)

  # l_sparsity = verbose_sparsity(gl)
  # f_sparsity = verbose_sparsity(gf)
  # if (!is.null(l_sparsity)) {
  #   entry_string = paste0(entry_string, sprintf("%10.3f", l_sparsity))
  # }
  # if (!is.null(f_sparsity)) {
  #   entry_string = paste0(entry_string, sprintf("%10.3f", f_sparsity))
  # }

  if (!is.null(max_chg)) {
    entry_string = paste0(entry_string,
                          sprintf("%15.2f", 100 * max_chg), "%")
  }

  if (!is.null(obj)) {
    entry_string = paste0(entry_string,
                          sprintf("%16.2f", obj))
  }

  if (stopping_rule == "objective") {
    diff_string = formatC(obj_diff, format="e", digits=2)
    entry_string = paste0(entry_string,
                          sprintf("%11s", diff_string))
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

# # At present, only returns nonnull for ebnm_pn.
# verbose_sparsity = function(g) {
#   if (is.null(g[[1]]$pi0)) {
#     return(NULL)
#   } else {
#     return(mean(sapply(g, function(k) {k$pi0})))
#   }
# }
