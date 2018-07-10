test_that("passing in and saving ebnm parameters works", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  f = flash_add_greedy(Y, Kmax=1, nullcheck=FALSE, ebnm_fn="ebnm_ash")
  flash_defaults = list(output="flash_data", mixcompdist="normal", method="shrink")
  expect_equal(f$ebnm_param_l[[1]], flash_defaults)
  expect_equal(f$ebnm_param_f[[1]], flash_defaults)
  ebnm_param_l = list(mixcompdist="uniform", method="shrink", g=f$gl[[1]],
                      fixg=TRUE)
  f = flash_add_greedy(Y, Kmax=1, nullcheck=FALSE,
                       ebnm_fn=list(l="ebnm_ash", f="ebnm_pn"),
                       ebnm_param=list(l=ebnm_param_l, f=list()))
  expect_equal(f$ebnm_param_l[[1]], c(ebnm_param_l, output="flash_data"))
  expect_equal(f$ebnm_param_f[[1]], list())
})
