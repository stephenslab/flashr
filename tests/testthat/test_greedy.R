test_that("greedy works a simple example as expected", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  data = flash_set_data(Y)  # note that some of these expectations might fail for some seeds
  f = flash_add_greedy(data,2)  # they are just expected based on the true model
  expect_equal(get_k(f),2)  # and do pass for the seed i selected
  f2 = flash_add_greedy(data,3,f_init=f)
  expect_equal(flash_get_objective(data,f),flash_get_objective(data,f2))
}
)
