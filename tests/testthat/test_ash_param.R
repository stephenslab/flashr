test_that("passing in and saving ash parameters works", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  f = flash_r1(Y)
  expect_equal(f$ash_param_l[[1]],flash_default_ash_param())
  expect_equal(f$ash_param_f[[1]],flash_default_ash_param())
  f = flash_r1(Y,ash_param=list(method="fdr"))
  expect_equal(f$ash_param_l[[1]],
               modifyList(flash_default_ash_param(),list(method="fdr")))
}
)
