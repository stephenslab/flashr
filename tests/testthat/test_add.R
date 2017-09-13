test_that("adding factor based on svd does same as using 2 svd factors", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  data = set_flash_data(Y)
  f1 = flash_init(data,1,"svd")
  f2 = flash_init(data,2,"svd")
  f3 = flash_add_factor(data,f1,1)
  expect_equal(flash_get_lf(f3),flash_get_lf(f2))
}
)
