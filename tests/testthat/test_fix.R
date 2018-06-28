test_that("fixing a loading it does not change", {
  set.seed(1)
  LL = matrix(1,nrow=10,ncol=1)
  FF = matrix(rnorm(5),nrow=5,ncol=1)
  fixl = matrix(TRUE,nrow=10,ncol=1)
  f = flash_init_lf(LL,FF,fixl = fixl)
  Y = matrix(rnorm(50),nrow=10)
  data = flash_set_data(Y)
  f = flash_backfit(data,f,1,"constant")
  expect_equal(as.vector(f$EL),rep(1,10))
  fixl[1:4,]=FALSE
  f = flash_init_lf(LL,FF,fixl = fixl)
  f = flash_backfit(data,f,1,"constant")
  expect_equal(as.vector(f$EL)[5:10],rep(1,6))

  fixf = matrix(TRUE, nrow=5,ncol=1)
  f = flash_init_lf(LL,FF=matrix(rep(1,5),nrow=5),fixf = fixf)
  f = flash_backfit(data,f,1,"constant")
  expect_equal(as.vector(f$EF),rep(1,5))
  fixf[3:5,] = FALSE
  f = flash_init_lf(LL,FF=matrix(rep(1,5),nrow=5),fixf = fixf)
  f = flash_backfit(data,f,1,"constant")
  expect_equal(as.vector(f$EF)[1:2],rep(1,2))
})
