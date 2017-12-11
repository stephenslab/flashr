test_that("passing in and saving ebnm parameters works", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  f = flash_r1(Y,nullcheck=FALSE)
  expect_equal(f$ebnm_param_l[[1]],flash_default_ebnm_param())
  expect_equal(f$ebnm_param_f[[1]],flash_default_ebnm_param())
  f = flash_r1(Y,ebnm_param=list(method="fdr"),nullcheck=FALSE)
  expect_equal(f$ebnm_param_l[[1]],
               modifyList(flash_default_ebnm_param(),list(method="fdr")))
}
)
