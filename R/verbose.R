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
  if (stopping_rule == "objective") {
    rule_string = "difference in obj. is"
    tol_string = formatC(tol, format = "e", digits = 2)
  } else {
    tol_string = paste0(100 * tol, "%")
    if (stopping_rule == "loadings") {
      rule_string = "max loading change is"
    } else if (stopping_rule == "factors") {
      rule_string = "max factor change is"
    } else { # stopping_rule == "any_param"
      rule_string = "max parameter change is"
    }
  }
  return(paste("stop when", rule_string, "<", tol_string))
}

verbose_obj_table_header = function(verbose_output) {
  header_string = "  Iteration"
  if ("l" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%9s", "pi0 (l)"))
  }
  if ("f" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%9s", "pi0 (f)"))
  }
  if ("L" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%13s", "Max Chg (l)"))
  }
  if ("F" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%13s", "Max Chg (f)"))
  }
  if ("o" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%15s", "Objective"))
  }
  if ("d" %in% verbose_output) {
    header_string = paste0(header_string,
                           sprintf("%11s", "Obj Diff"))
  }
  message(header_string)
}

verbose_obj_table_entry = function(verbose_output, iter, obj, obj_diff,
                                   max_chg_l, max_chg_f, gl, gf) {
  entry_string = sprintf("%11d", iter)
  if ("l" %in% verbose_output) {
    l_sparsity = verbose_sparsity(gl)
    entry_string = paste0(entry_string, sprintf("%9s", l_sparsity))
  }
  if ("f" %in% verbose_output) {
    f_sparsity = verbose_sparsity(gf)
    entry_string = paste0(entry_string, sprintf("%9s", f_sparsity))
  }
  if ("L" %in% verbose_output) {
    entry_string = paste0(entry_string,
                          sprintf("%12.2f", 100 * max_chg_l), "%")
  }
  if ("F" %in% verbose_output) {
    entry_string = paste0(entry_string,
                          sprintf("%12.2f", 100 * max_chg_f), "%")
  }
  if ("o" %in% verbose_output) {
    entry_string = paste0(entry_string,
                          sprintf("%15.2f", obj))
  }
  if ("d" %in% verbose_output) {
    diff_string = formatC(obj_diff, format="e", digits=2)
    entry_string = paste0(entry_string,
                          sprintf("%11s", diff_string))
  }
  message(entry_string)
}

# At present, only returns nonnull for ebnm_pn.
verbose_sparsity = function(g) {
  if (is.null(g[[1]]$pi0)) {
    return("NA")
  } else {
    s = mean(sapply(g, function(k) {k$pi0}))
    return(formatC(s, format = "f", digits = 3))
  }
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
