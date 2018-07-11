test_that("argument checking works", {
  set.seed(1)
  l1 = 5*c(rnorm(2), rep(0, 3))
  l2 = 5*c(rep(0, 3), rnorm(2))
  f1 = c(rnorm(10), rep(0, 10))
  f2 = c(rep(0, 10), rnorm(10))
  LF = outer(l1,f1) + outer(l2,f2)
  Y = LF + rnorm(5*20)
  f = flash(Y, 2, nullcheck=F)

  # kset:
  expect_identical(handle_kset(NULL, f), 1:2)
  expect_identical(handle_kset(1:2, f), 1:2)
  expect_error(handle_kset(1:3, f))
  expect_error(handle_kset(list(1, 2), f))
  nullf = flash_init_null()
  expect_error(handle_kset(NULL, nullf))

  # init_fn:
  expect_identical(handle_init_fn("udv_si"), "udv_si")
  expect_identical(handle_init_fn(udv_si), udv_si)
  expect_error(handle_init_fn("sdkfjsldkfjlsj"))

  # ebnm_fn:
  expect_identical(handle_ebnm_fn("ebnm_pn"),
                   list(l="ebnm_pn", f="ebnm_pn"))
  expect_identical(handle_ebnm_fn(list(l="ebnm_pn", f="ebnm_ash")),
                   list(l="ebnm_pn", f="ebnm_ash"))
  expect_error(handle_ebnm_fn(ebnm_ash))
  expect_error(handle_ebnm_fn(list(l="ebnm_ash")))
  expect_error(handle_ebnm_fn(c("ebnm_ash", "ebnm_pn")))
  expect_error(handle_ebnm_fn("sfdskjfhskdjhfk"))

  # ebnm_param:
  ash_defaults=list(mixcompdist="normal",
                    method="shrink")
  # defaults:
  expect_identical(handle_ebnm_param(NULL,
                                     handle_ebnm_fn("ebnm_ash"),
                                     1),
                   list(l=list(ash_defaults), f=list(ash_defaults)))
  expect_identical(handle_ebnm_param(NULL,
                                     handle_ebnm_fn("ebnm_pn"),
                                     2),
                   list(l=list(list(), list()), f=list(list(), list())))
  new_ash=list(mixcompdist="uniform", method="shrink")
  # different ebnm_param for l and f:
  expect_identical(handle_ebnm_param(list(l=list(mixcompdist="uniform"),
                                          f=list()),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l=list(new_ash, new_ash),
                        f=list(ash_defaults, ash_defaults)))
  expect_identical(handle_ebnm_param(list(l=list(),
                                          f=list(mixcompdist="uniform")),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l=list(ash_defaults, ash_defaults),
                        f=list(new_ash, new_ash)))
  # different ebnm_param for each loading/factor:
  expect_identical(handle_ebnm_param(list(list(mixcompdist="uniform"),
                                          list()),
                                     handle_ebnm_fn("ebnm_ash"),
                                     2),
                   list(l=list(new_ash, ash_defaults),
                        f=list(new_ash, ash_defaults)))
  # if specify either l or f then must specify both:
  expect_error(handle_ebnm_param(list(l=list(mixcompdist="uniform")),
                                 handle_ebnm_fn("ebnm_ash"),
                                 1))
  # if specify for different factors then must specify for all:
  expect_error(handle_ebnm_param(list(list(mixcompdist="uniform"),
                                      list()),
                                 handle_ebnm_fn("ebnm_ash"),
                                 3))

  expect_error(handle_ebnm_param(c(mixcompdist="uniform", method="fdr")))
})

