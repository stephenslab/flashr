context("fixed factors and loadings")

set.seed(1)

# rank-one matrix with sparse loading and some missing data:
LL = c(rep(1, 3), rep(0, 7))
FF = rep(1, 10)
Y = outer(LL, FF) + rnorm(100)
Y[2, 2] = NA
Y[7, 8] = NA

test_that("a loading with all elements fixed is successfully added", {
  fl = flash(Y, fixed_loadings = LL, var_type = "constant")
  expect_equivalent(fl$fit$EL, LL)
})

test_that("a loading with not all elements fixed is successfully added", {
  is_fixed = c(rep(TRUE, 5), rep(FALSE, 5))
  fl = flash(Y,
             fixed_loadings = list(vals = LL, is_fixed = is_fixed),
             var_type = "constant")
  expect_true(all(fl$fit$EL[is_fixed, 1] == LL[is_fixed]))
  expect_false(any(fl$fit$EL[!is_fixed, 1] == LL[!is_fixed]))
})

test_that("a loading with missing elements is successfully added", {
  missing = c(rep(FALSE, 2), rep(TRUE, 3), rep(FALSE, 5))
  LL_missing = LL
  LL_missing[missing] = NA
  fl = flash(Y, fixed_loadings = LL_missing, var_type = "constant")
  expect_true(all(fl$fit$EL[!missing, 1] == LL[!missing]))
})

test_that("a sparse loading is successfully added", {
  LL_sparse = c(rep(NA, 3), rep(0, 7))
  fl = flash(Y, fixed_loadings = LL_sparse, var_type = "constant")
  expect_true(all(fl$fit$EL[!is.na(LL_sparse), 1] == 0))
  expect_false(any(fl$fit$EL[is.na(LL_sparse), 1] == 0))
})

test_that("a very sparse loading is successfully added", {
  LL_sparse = c(NA, rep(0, 9))
  fl = flash(Y, fixed_loadings = LL_sparse, var_type = "constant")
  expect_true(all(fl$fit$EL[!is.na(LL_sparse), 1] == 0))
  expect_false(any(fl$fit$EL[is.na(LL_sparse), 1] == 0))
})

test_that("a fixed loading is successfully added to an existing object", {
  # TODO
})

# TODO : test factors

LL = rep(1, 10)
fixl = c(rep(TRUE, 5), rep(FALSE, 5))
fl = flash(Y, fixed_loadings = list(vals = LL, is_fixed = fixl),
           greedy_Kmax = 0)

test_that("a fixed loading does not change during fitting", {
  expect_setequal(fl$fit$EL[1:5, 1], 1)
})

FF = 1:5
fixf = c(FALSE, rep(TRUE, 4))
fl = flash(Y, f_init = fl, fixed_factors = list(vals = FF, is_fixed = fixf),
           greedy_Kmax = 0, var_type = "constant")

test_that("a fixed factor does not change during fitting", {
  expect_equal(fl$fit$EF[2:5, 2], 2:5)
})

LF = outer(rep(5, 5), c(rep(1, 5), rep(0, 15)))
Y = LF + rnorm(5 * 20)
FF = c(rep(NA, 2), rep(0, 18))
fl = flash(Y, fixed_factors = FF, greedy_Kmax = 0, var_type = "constant")

test_that("adding sparse factors works as expected", {
  expect_setequal(fl$ldf$f[3:20, 1], 0)
  expect_false(any(fl$ldf$f[1:2, 1] == 0))
})

LL = c(rep(NA, 2), rep(0, 3))
fl2 = flash(Y, f_init = fl, fixed_loadings = LL, greedy_Kmax = 0,
            var_type = "constant")

test_that("adding sparse loadings works as expected", {
  expect_setequal(fl2$ldf$l[3:5, 2], 0)
  expect_false(any(fl$ldf$l[1:2, 1] == 0))
})

Y[, 1] = Y[, 1] + 50
FF = c(NA, rep(0, 19))
fl = flash(Y, fixed_factors = FF, greedy_Kmax = 0, var_type = "constant")

test_that("factors with only one missing element can be added", {

})
