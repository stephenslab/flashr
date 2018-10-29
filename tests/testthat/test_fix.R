context("fixing factors and loadings")

set.seed(1)

Y = matrix(rnorm(50), nrow = 10, ncol = 5)

LL = rep(1, 10)
fixl = c(rep(TRUE, 5), rep(FALSE, 5))
fl = flash(Y, fixed_loadings = list(vals = LL, is_fixed = fixl),
           greedy_Kmax = 0, nullcheck = FALSE)

test_that("a fixed loading does not change during fitting", {
  expect_setequal(fl$fit$EL[1:5, 1], 1)
})

FF = 1:5
fixf = c(FALSE, rep(TRUE, 4))
fl = flash(Y, f_init = fl, fixed_factors = list(vals = FF, is_fixed = fixf),
           greedy_Kmax = 0, nullcheck = FALSE)

test_that("a fixed factor does not change during fitting", {
  expect_equal(fl$fit$EF[2:5, 2], 2:5)
})

LF = outer(rep(5, 5), c(rep(1, 5), rep(0, 15)))
Y = LF + rnorm(5 * 20)
FF = c(rep(NA, 5), rep(0, 15))
fl = flash(Y, fixed_factors = list(vals = FF), greedy_Kmax = 0)

test_that("adding sparse factors works as expected", {
  expect_setequal(fl$ldf$f[6:20, 1], 0)
  expect_false(any(fl$ldf$f[1:5, 1] == 0))
})
