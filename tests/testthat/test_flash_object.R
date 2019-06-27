context("test_flash_object")

test_that("flash object interface works as expected", {
  set.seed(666)
  n = 10
  p = 50
  k = 5

  LL = matrix(rnorm(n*k), nrow=n, ncol=k)
  FF = matrix(0, nrow=p, ncol=k)
  for (k in 1:5) {
    FF[(1 + 5*(k - 1)):(5*(k + 1)), k] = 1
  }
  Y = LL %*% t(FF) + rnorm(n * p)

  fo = flash_add_greedy(Y, Kmax = 2)
  expect_s3_class(fo, "flash")
  expect_equal(length(fo$fit_history), 2)

  fo = flash_backfit(Y, fo)
  expect_s3_class(fo, "flash")
  expect_equal(length(fo$fit_history), 3)

  fo = flash_add_fixed_loadings(Y, rep(1, 10), fo, backfit = FALSE)
  expect_identical(fo$objective, NA)

  fo = flash_zero_out_factor(Y, fo, 2)
  expect_identical(fo$objective, NA)

  fo = flash_add_fixed_factors(Y, c(rep(1, 10), rep(0, 40)), fo,
                               nullcheck = FALSE)
  expect_equal(fo$nfactors, 2) # two are zero factors
  expect_equal(length(fo$fit_history), 4)

  fo = flash_zero_out_factor(Y, fo, 1)
  expect_s3_class(fo, "flash")
  expect_equal(fo$nfactors, 1)
  expect_false(is.na(fo$objective))

  fo = flash_add_factors_from_data(Y, 2, fo)
  expect_s3_class(fo, "flash")
  expect_equal(length(fo$fit_history), 6) # nullcheck causes repeated backfit

  fo = flash_add_greedy(Y, 1, fo)
  expect_s3_class(fo, "flash")
  expect_equal(length(fo$fit_history), 7)
  expect_equal(fo$fit_history[[7]]$zeroed_out, 7)

  Ymiss_idx = sample(1:(n * p), floor(n * p / 4), replace = FALSE)
  Y[Ymiss_idx] = NA
  data = flash_set_data(Y)
  fo2 = flash_add_factors_from_data(data, 3)
  Yfill = flash_fill(data, fo2)
  expect_equal(dim(Yfill), c(n, p))
})
