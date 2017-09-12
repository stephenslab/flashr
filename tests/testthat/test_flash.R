test_that("flash update improves mean squared error in simple situation", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  f = flash_set_data(Y)
  f = flash_init_random(f,1)
  e1 = mean((LF- flash_get_lf(f))^2) # compute mean squared error
  f = flash_update_single_factor(f,1)
  f = flash_update_single_loading(f,1)
  e2 = mean((LF- flash_get_lf(f))^2) # compute mean squared error
  expect_true(e1>e2)
}
)
