test_that("test R2 computation matches R^2 when EL2=EL^2 and EF2=EF^2 ", {
  set.seed(1)
  Y = matrix(nrow=5,ncol=20,rnorm(100))
  f = flash_init_svd(Y,1)
  R = Y- flash_get_lf(f)
  R2 = R^2
  expect_equal(R2,get_R2(Y,f))
}
)
