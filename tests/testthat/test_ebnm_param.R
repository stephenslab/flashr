context("test_ebnm_param")

test_that("passing in and saving ebnm parameters works", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  f = flash_add_greedy(Y, Kmax=1, nullcheck=FALSE, ebnm_fn="ebnm_ash")$fit
  flash_defaults = list(mixcompdist="normal", method="shrink")
  expect_equal(f$ebnm_param_l[[1]], flash_defaults)
  expect_equal(f$ebnm_param_f[[1]], flash_defaults)
  ebnm_param_l = list(mixcompdist="uniform", method="shrink", g=f$gl[[1]],
                      fixg=TRUE)
  f = flash_add_greedy(Y, Kmax=1, nullcheck=FALSE,
                       ebnm_fn=list(l="ebnm_ash", f="ebnm_pn"),
                       ebnm_param=list(l=ebnm_param_l, f=list()))$fit
  expect_equal(f$ebnm_fn_l[[1]], "ebnm_ash")
  expect_equal(f$ebnm_fn_f[[1]], "ebnm_pn")
  expect_equal(f$ebnm_param_l[[1]], ebnm_param_l)
  expect_equal(f$ebnm_param_f[[1]], list(warmstart = TRUE))
})
