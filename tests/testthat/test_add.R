context("test_add")

test_that("various additions work", {
  set.seed(1)
  l = rep(1,5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  l_names = c("a", "b", "c", "d", "e")
  rownames(Y) = l_names

  data = flash_set_data(Y)
  expect_identical(rownames(data$Yorig), l_names)
  expect_identical(rownames(data$missing), l_names)
  expect_identical(rownames(data$Y), l_names)

  f1 = flash_init_fn(data,"udv_svd",1)
  expect_s3_class(f1, "flash_fit")
  expect_identical(rownames(f1$EL), l_names)
  expect_identical(rownames(f1$EL2), l_names)
  expect_identical(rownames(f1$fixl), l_names)

  f2 = flash_init_fn(data,"udv_svd",2)
  f3 = flash_add_factors_from_data(data,1,f1,"udv_svd",backfit=FALSE)
  # add_factors_from_data is "public", so it produces a flash object:
  expect_s3_class(f3, "flash")
  expect_equal(flash_get_fitted_values(f3),flash_get_fitted_values(f2))

  f2 = flash_init_null()
  f3 = flash_combine(f2,f1)
  expect_equal(f3,f1)
  expect_s3_class(f2, "flash_fit")
  expect_s3_class(f3, "flash_fit")

  fixed_l = matrix(1, nrow=5, ncol=1)
  rownames(fixed_l) = l_names
  f1 = flash_add_fixed_loadings(data, fixed_l)
  expect_identical(rownames(f1$fit$EL), l_names)
  f1 = flash_backfit(data,f1)
  expect_s3_class(f1, "flash")
  expect_identical(rownames(f1$fit$EL), l_names)
  expect_equal(f1$fit$EL, fixed_l)

  f_names = as.character(20:1)

  l = rnorm(20)
  f = rep(1,20)
  LF = outer(l,f)
  Y = LF + rnorm(20*20)
  colnames(Y) = f_names

  data = flash_set_data(Y)
  expect_identical(colnames(data$Yorig), f_names)
  expect_identical(colnames(data$missing), f_names)
  expect_identical(colnames(data$Y), f_names)

  fixed_f = matrix(1, nrow=20, ncol=1)
  rownames(fixed_f) = f_names
  f1 = flash_add_fixed_factors(data, fixed_f)
  expect_equal(f1$fit$EF, fixed_f)
  expect_s3_class(f1, "flash")
  expect_identical(rownames(f1$fit$EF), f_names)
  expect_identical(rownames(f1$fit$EF2), f_names)
  expect_identical(rownames(f1$fit$fixf), f_names)
})

