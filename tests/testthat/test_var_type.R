context("var_type")

set.seed(1)

n = 100
p = 5

sd = 2
Y_r0 = matrix(rnorm(n * p, sd = sd), nrow = n, ncol = p)

test_that("zero and constant var_types work with no missing data", {
  fl_zero = flash(flash_set_data(Y_r0, S = sd), Kmax = 1, var_type = "zero")
  expect_equal(fl_zero$fit$tau, 1 / sd^2)

  fl_const = flash(Y_r0, Kmax = 1, var_type = "constant")
  expect_equal(fl_const$fit$tau, 1 / sd(Y_r0)^2, tol = 0.01)
})

missing = matrix(FALSE, nrow = n, ncol = p)
missing[1, 1] = TRUE
Y_r0[missing] = NA

test_that("zero and constant var_types work with missing data", {
  fl_zero = flash(flash_set_data(Y_r0, S = sd), Kmax = 1, var_type = "zero")
  expect_setequal(fl_zero$fit$tau[missing], 0)
  expect_setequal(fl_zero$fit$tau[!missing], 1 / sd^2)

  fl_const = flash(Y_r0, Kmax = 1, var_type = "constant")
  expect_setequal(fl_const$fit$tau[missing], 0)
  expect_setequal(fl_const$fit$tau[!missing], fl_const$fit$tau[!missing][1])
  expect_equal(fl_const$fit$tau[!missing][1], 1 / sd(Y_r0, na.rm = TRUE)^2,
               tol = 0.01)
})

sd = matrix(rep(1:p, n), nrow = n, ncol = p, byrow = TRUE)
Y_r0 = matrix(rnorm(n * p, sd = sd), nrow = n, ncol = p)

test_that("by_column and by_row var_types work with no missing data", {
  fl_col = flash(Y_r0, Kmax = 1, var_type = "by_column")
  expect_equivalent(fl_col$fit$tau, 1 / apply(Y_r0, 2, sd)^2, tol = 0.01)

  fl_row = flash(t(Y_r0), Kmax = 1, var_type = "by_row")
  expect_equal(fl_row$fit$tau, t(fl_col$fit$tau), tol = 0.01)
})

Y_r0[missing] = NA

test_that("by_column and by_row var_types work with missing data", {
  fl_col = flash(Y_r0, Kmax = 1, var_type = "by_column")
  expected_tau = matrix(1 / apply(Y_r0, 2, function(x) sd(x, na.rm = TRUE))^2,
                        nrow = n, ncol = p, byrow = TRUE)
  expected_tau[missing] = 0
  expect_equal(fl_col$fit$tau, expected_tau, tol = 0.01)

  fl_row = flash(t(Y_r0), Kmax = 1, var_type = "by_row")
  expect_equal(fl_row$fit$tau, t(fl_col$fit$tau), tol = 0.01)
})

# TODO: move this to set_data tests
# check for error when trying to pass bad arguments to S
# expect_error(flash_set_data(Y, S=1:5))
# expect_error(flash_set_data(Y, matrix(1, nrow=2, ncol=2)))
