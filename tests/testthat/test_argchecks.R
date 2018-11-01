context("argument checking")

set.seed(1)

Y = matrix(5, nrow = 5, ncol = 10) + rnorm(50)
fl = flash(Y, greedy_Kmax = 1, var_type = "constant", nullcheck = FALSE)
fl2 = flash(Y, f_init = fl, fixed_loadings = rep(1, 5),
            var_type = "constant", nullcheck = FALSE)

test_that("bad inputs to Y and S are rejected", {
  Y = rep(1, 4) # vector
  expect_error(flash_set_data(Y, S = NULL))

  Y = matrix(Y, nrow = 2, ncol = 2)
  Y[1, 1] = Inf
  expect_error(flash_set_data(Y, S = NULL))

  Y[1, 1] = NaN
  expect_error(flash_set_data(Y, S = NULL))

  Y[1, 1] = NA # OK
  expect_is(flash_set_data(Y, S = NULL), "flash_data")

  S = 1 # OK to pass in scalar
  expect_is(flash_set_data(Y, S), "flash_data")

  S = rep(1, 2) # vector
  expect_error(flash_set_data(Y, S))

  S = matrix(1, nrow = 2, ncol = 3) # incorrect dimensions
  expect_error(flash_set_data(Y, S))

  S = matrix(1, nrow = 2, ncol = 2)
  expect_is(flash_set_data(Y, S), "flash_data")

  S[2, 2] = -1
  expect_error(flash_set_data(Y, S))

  S[2, 2] = NA
  expect_error(flash_set_data(Y, S))

  S[2, 2] = Inf # OK
  expect_is(flash_set_data(Y, S), "flash_data")
})

test_that("Y is checked against f_init", {
  bad_Y = matrix(rnorm(25), nrow = 5, ncol = 5)
  expect_error(flash(bad_Y, f_init = fl))
})

test_that("current restrictions on var_type are respected", {
  expect_error(flash(Y, greedy_Kmax = 1, var_type = "zero"))
  expect_error(flash(Y, greedy_Kmax = 1, var_type = "kroneker"))
  expect_error(flash(Y, S = 1, greedy_Kmax = 1, var_type = "constant"))
})

test_that("init_fn handling works as expected", {
  expect_identical(handle_init_fn("udv_si"), "udv_si")
  expect_identical(handle_init_fn(udv_si), udv_si)
  expect_error(handle_init_fn("sdkfjsldkfjlsj"))
})

test_that("ebnm_fn handling works as expected", {
  expect_identical(handle_ebnm_fn("ebnm_pn"),
                   list(l = "ebnm_pn", f = "ebnm_pn"))
  expect_identical(handle_ebnm_fn(list(l = "ebnm_pn", f = "ebnm_ash")),
                   list(l = "ebnm_pn", f = "ebnm_ash"))
  expect_error(handle_ebnm_fn(ebnm_ash))
  expect_error(handle_ebnm_fn(list(l = "ebnm_ash")))
  expect_error(handle_ebnm_fn(c("ebnm_ash", "ebnm_pn")))
  expect_error(handle_ebnm_fn("sfdskjfhskdjhfk"))
})

test_that("kset handling works as expected", {
  fit = fl2$fit

  expect_identical(handle_kset(NULL, fit), 1:2)
  expect_identical(handle_kset(1:2, fit), 1:2)

  expect_error(handle_kset(1:3, fit))
  expect_error(handle_kset(list(1, 2), fit))

  nullf = flash_init_null()
  expect_error(handle_kset(NULL, nullf))
})

test_that("ebnm_param handling works as expected", {
  # ebnm_param:
  ash_defaults = list(mixcompdist = "normal", method = "shrink")
  pn_defaults = list(warmstart = TRUE)

  # defaults:
  expect_identical(handle_ebnm_param(NULL,
                                     handle_ebnm_fn("ebnm_ash"),
                                     1),
                   list(l = list(ash_defaults), f = list(ash_defaults)))
  expect_identical(handle_ebnm_param(NULL,
                                     handle_ebnm_fn("ebnm_pn"),
                                     2),
                   list(l = list(pn_defaults, pn_defaults),
                        f = list(pn_defaults, pn_defaults)))
  new_ash = list(mixcompdist = "uniform", method = "shrink")

  # different ebnm_param for l and f:
  expect_identical(handle_ebnm_param(list(l = list(mixcompdist = "uniform"),
                                          f = list()),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l = list(new_ash, new_ash),
                        f = list(ash_defaults, ash_defaults)))
  expect_identical(handle_ebnm_param(list(l = list(),
                                          f = list(mixcompdist = "uniform")),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l = list(ash_defaults, ash_defaults),
                        f = list(new_ash, new_ash)))

  # different ebnm_param for each loading/factor:
  expect_identical(handle_ebnm_param(list(list(mixcompdist = "uniform"),
                                          list()),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l = list(new_ash, ash_defaults),
                        f = list(new_ash, ash_defaults)))

  # if specify either l or f then must specify both:
  expect_error(handle_ebnm_param(list(l = list(mixcompdist = "uniform")),
                                 handle_ebnm_fn("ebnm_ash"),
                                 1))

  # if specify for different factors then must specify for all:
  expect_error(handle_ebnm_param(list(list(mixcompdist = "uniform"),
                                      list()),
                                 handle_ebnm_fn("ebnm_ash"),
                                 3))

  expect_error(handle_ebnm_param(c(mixcompdist = "uniform", method = "fdr")))
})

test_that("ebnm_param only accepts lists of lists when appropriate", {
  ebnm_param = list(list(warmstart = TRUE), list(warmstart = FALSE))
  expect_error(flash(Y, f_init = fl2, greedy_Kmax = 1, backfit = TRUE,
                     control = list(ebnm_param = ebnm_param)))
  expect_error(flash(Y, f_init = fl2, fixed_factors = rep(1, 10),
                     backfit = TRUE,
                     control = list(ebnm_param = ebnm_param)))
  expect_error(flash(Y, f_init = fl2, fixed_factors = rep(1, 10),
                     fixed_loadings = 1:5,
                     control = list(ebnm_param = ebnm_param)))

  fl3 = flash(Y, f_init = fl2, greedy_Kmax = 2,
              control = list(ebnm_param = ebnm_param))
  expect_is(fl3, "flash")

  fl3 = flash(Y, f_init = fl2, backfit = TRUE,
              control = list(ebnm_param = ebnm_param))
  expect_is(fl3, "flash")

  fl3 = flash(Y, f_init = fl2, fixed_factors = matrix(1:20, ncol = 2),
              control = list(ebnm_param = ebnm_param))
  expect_is(fl3, "flash")
})
