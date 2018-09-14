test_that("various var_types work as expected", {
  set.seed(1)
  l = rep(5,5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + 2*rnorm(5*20)

  # var_type = zero
  data <- flash_set_data(Y, S=2)
  fl <- flash_add_greedy(data, 1, var_type="zero")$fit
  expect_identical(fl$tau, matrix(0.25, nrow=nrow(Y), ncol=ncol(Y)))

  # var_type = by_column, by_row
  Y = LF + c(rnorm(5 * 10), 2 * rnorm(5 * 10))
  flcol <- flash_add_greedy(Y, 1, var_type="by_column")$fit
  expect_setequal(flcol$tau[2:nrow(Y), 1], flcol$tau[1, 1])
  flrow <- flash_add_greedy(t(Y), 1, var_type="by_row")$fit
  expect_equal(flcol$tau, t(flrow$tau), tolerance=0.001)
})
