context("var_type")

set.seed(1)

test_that("zero and constant var_types work as expected", {
  # rank-zero data with sd = 2
  Y = matrix(2 * rnorm(5 * 20), nrow = 5, ncol = 20)
  fl = flash(Y, S = 2, var_type = "zero", greedy_Kmax = 1)
  expect_equal(fl$fit$tau, 0.25)

  flconst = flash(Y, var_type = "constant", greedy_Kmax = 1)
  expect_equal(flconst$fit$tau, 1 / sd(Y)^2, tol = 0.01)
})

test_that("by_column and by_row var_types work as expected", {
  # rank-zero data: first column has sd = 1; second has sd = 2
  Y = matrix(rnorm(100 * 2, sd = rep(c(1, 2), each = 100)),
             nrow = 100, ncol = 2)

  flcol = flash(Y, var_type = "by_column", greedy_Kmax = 1)
  expect_setequal(flcol$fit$tau[, 1], flcol$fit$tau[1, 1])
  expect_equal(flcol$fit$tau[1, 1], 1 / sd(Y[, 1])^2, tol = 0.01)

  flrow = flash(t(Y), var_type="by_row", greedy_Kmax = 1)
  expect_equal(t(flrow$fit$tau), flcol$fit$tau, tol = 0.001)
})

# TODO: move to param check tests
# check for error when trying to pass bad arguments to S
# expect_error(flash_set_data(Y, S=1:5))
# expect_error(flash_set_data(Y, matrix(1, nrow=2, ncol=2)))
