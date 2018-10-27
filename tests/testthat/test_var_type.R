context("var_type")

set.seed(1)

n = 50
p = 5
sd = 2
l = rep(5, n)
f = rnorm(p)
LF = outer(l, f)
E = rnorm(n * p)
Y = LF + sd * E

test_that("zero var_type works as expected", {
  data = flash_set_data(Y, S = sd)
  fl = flash_add_greedy(data, 1, var_type = "zero")
  expect_equal(fl$fit$tau, 1 / sd^2)
})

test_that("constant var_type works as expected", {
  fl = flash_add_greedy(Y, 1, var_type = "constant")
  expect_equivalent(fl$fit$tau, 1 / sd^2, tol = 0.05)
})

test_that("by_column var_type works as expected", {
  Y = LF + rep(1:p, each = n) * E
  fl = flash_add_greedy(Y, 1, var_type = "by_column")
  expect_equivalent(fl$fit$tau, 1 / (1:p)^2, tol = 0.1)
})

test_that("by_row var_type works as expected", {
  Y = LF + rep(1:n) * E
  fl = flash_add_greedy(Y, 1, var_type = "by_row")
  expect_equivalent(fl$fit$tau, 1 / (1:n)^2, tol = 0.1)
})

# TODO: move this to set_data tests
# check for error when trying to pass bad arguments to S
# expect_error(flash_set_data(Y, S=1:5))
# expect_error(flash_set_data(Y, matrix(1, nrow=2, ncol=2)))
