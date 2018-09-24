test_that("flash update improves mean squared error in simple situation", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  data = flash_set_data(Y)
  f = flash_init_fn(data,udv_random,1)
  e1 = mean((LF - flash_get_fitted_values(f))^2) # compute mean squared error
  f = flash_update_precision(data,f,"by_column")
  Rk = flash_get_Rk(data, f, 1)
  f = flash_update_single_factor(data,f,1,"ebnm_pn",NULL,Rk)
  f = flash_update_single_loading(data,f,1,"ebnm_pn",NULL,Rk)
  e2 = mean((LF - flash_get_fitted_values(f))^2) # compute mean squared error
  expect_true(e1 > e2)
})




# just used to run flash.. at the moment don't have a good test!
# test_that("flash r1 approximates solution", {
#   set.seed(10)
#   ll = rnorm(20)
#   ff = rnorm(100)
#   LF = outer(ll,ff)
#   Y = LF + rnorm(20*100)
#   data = flash_set_data(data)
#   f = flash_r1(data,"random")
#   f = flash_r1(data,"svd")
# }
# )
