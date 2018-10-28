context("greedy")

set.seed(1)
l = rnorm(5)
f = rnorm(20)
LF = outer(l, f)
Y = LF + rnorm(5 * 20)

# Only one factor should be added here:
fl = flash(Y, greedy_Kmax = 2)
# An additional call to flash_add_greedy should not add any more factors:
fl2 = flash(Y, greedy_Kmax = 1, f_init = fl)

test_that("greedy works a simple rank-1 example as expected", {
  expect_equal(fl$nfactors, 1)
  expect_equal(fl2$nfactors, 1)
  expect_equal(fl$objective, fl2$objective)
})
