test_that("flash_fill works", {
  l = rep(1,5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)
  missing = rbinom(5*20, 1, .2)
  Y[missing == 1] = NA

  fl = flash(Y)
  Yhat = flash_fill(Y, fl)
  # non-missing entries are the same
  expect_identical(Y[missing == 0], Yhat[missing == 0])
  # no missing entries in Yhat
  expect_equal(sum(is.null(Yhat)), 0)
})
