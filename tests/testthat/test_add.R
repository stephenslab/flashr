test_that("various additions work", {
  set.seed(1)
  l = rep(1,5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  data = flash_set_data(Y)
  f1 = flash_init_fn(data,"udv_svd",1)
  f2 = flash_init_fn(data,"udv_svd",2)
  f3 = flash_add_factors_from_data(data,1,f1,"udv_svd")
  expect_equal(flash_get_lf(f3),flash_get_lf(f2))

  f2 = flash_init_null()
  f3 = flash_combine(f2,f1)
  expect_equal(f3,f1)

  f1 = flash_add_fixed_l(data, cbind(rep(1,5)))
  f1 = flash_backfit(data,f1)
  expect_equal(f1$EL,cbind(rep(1,5)))

  l = rnorm(20)
  f = rep(1,20)
  LF = outer(l,f)
  Y = LF + rnorm(20*20)
  data = flash_set_data(Y)
  f1 = flash_add_fixed_f(data, cbind(rep(1,20)))
  f1 = flash_backfit(data,f1,verbose=TRUE)
  expect_equal(f1$EF,cbind(rep(1,20)))

}
)
