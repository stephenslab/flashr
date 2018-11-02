get_default_verbose_output = function(stopping_rule) {
  return(switch(stopping_rule,
                factors = "odFn",
                loadings = "odLn",
                "odn"))
}

get_control_defaults <- function(method) {
  defaults = list(stopping_rule = "objective",
                  tol = 1e-2,
                  r1opt_maxiter = 500,
                  backfit_maxiter = 200,
                  nullcheck = TRUE,
                  seed = 1)

  return(c(defaults, get_method_defaults(method)))
}

get_method_defaults = function(method) {
  nn_param <- list(mixcompdist = "+uniform",
                   optmethod = "mixSQP")
  ash_param <- list(mixcompdist = "normal",
                    optmethod = "mixSQP")

  return(switch(method,
                fastest = list(init_fn = "udv_si",
                               ebnm_fn = "ebnm_pn",
                               ebnm_param = NULL),
                normalmix = list(init_fn = "udv_si",
                                 ebnm_fn = "ebnm_ash",
                                 ebnm_param = ash_param),
                nonnegative = list(init_fn = "udv_nn",
                                   ebnm_fn = "ebnm_ash",
                                   ebnm_param = list(l = nn_param,
                                                     f = nn_param)),
                nnloadings = list(init_fn = "udv_nnloadings",
                                  ebnm_fn = "ebnm_ash",
                                  ebnm_param = list(l = nn_param,
                                                    f = ash_param)),
                nnfactors = list(init_fn = "udv_nnfactors",
                                 ebnm_fn = "ebnm_ash",
                                 ebnm_param = list(l = ash_param,
                                                   f = nn_param))))
}
