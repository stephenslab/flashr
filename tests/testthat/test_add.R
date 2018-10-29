# TODO: move these

test_that("various additions work", {
  f2 = flash_init_fn(data,"udv_svd",2)
  f3 = flash_add_factors_from_data(data,1,f1,"udv_svd",backfit=FALSE)
  expect_s3_class(f3, "flash")
  expect_equal(flash_get_fitted_values(f3),flash_get_fitted_values(f2))

  f2 = flash_init_null()
  f3 = flash_combine(f2,f1)
  expect_equal(f3,f1)
})

