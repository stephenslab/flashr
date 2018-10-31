context("fixed factors and loadings")

set.seed(1)

# rank-one matrix with sparse loading and some missing data:
LL = c(rep(1, 3), rep(0, 7))
FF = rep(1, 10)
Y = outer(LL, FF) + rnorm(100)
Y[2, 2] = NA
Y[7, 8] = NA

test_that("fixed loadings and factors are successfully added", {
  fl = flash(Y, fixed_loadings = LL, var_type = "constant")
  expect_equivalent(fl$fit$EL, LL)

  fl2 = flash(t(Y), fixed_factors = LL, var_type = "constant")
  expect_equivalent(fl2$fit$EF, LL)
  expect_equal(fl$fit$EF, fl2$fit$EL)
})

test_that("fixed loadings and factors respect is_fixed parameter", {
  is_fixed = c(rep(TRUE, 5), rep(FALSE, 5))

  fl = flash(Y,
             fixed_loadings = list(vals = LL, is_fixed = is_fixed),
             var_type = "constant")
  expect_true(all(fl$fit$EL[is_fixed, 1] == LL[is_fixed]))
  expect_false(any(fl$fit$EL[!is_fixed, 1] == LL[!is_fixed]))

  fl2 = flash(t(Y),
              fixed_factors = list(vals = LL, is_fixed = is_fixed),
              var_type = "constant")
  expect_equal(fl2$fit$EF, fl$fit$EL)
  expect_equal(fl2$fit$EL, fl$fit$EF)
})

test_that("loadings and factors with missing elements are successfully added", {
  missing = c(rep(FALSE, 2), rep(TRUE, 3), rep(FALSE, 5))
  LL_missing = LL
  LL_missing[missing] = NA

  fl = flash(Y, fixed_loadings = LL_missing, var_type = "constant")
  expect_true(all(fl$fit$EL[!missing, 1] == LL[!missing]))

  fl2 = flash(t(Y), fixed_factors = LL_missing, var_type = "constant")
  expect_equal(fl2$fit$EF, fl$fit$EL)
  expect_equal(fl2$fit$EL, fl$fit$EF)
})

test_that("sparse loadings and factors are successfully added", {
  LL_sparse = c(rep(NA, 3), rep(0, 7))

  fl = flash(Y, fixed_loadings = LL_sparse, var_type = "constant")
  expect_true(all(fl$fit$EL[!is.na(LL_sparse), 1] == 0))
  expect_false(any(fl$fit$EL[is.na(LL_sparse), 1] == 0))

  fl2 = flash(t(Y), fixed_factors = LL_sparse, var_type = "constant")
  expect_equal(fl2$fit$EF, fl$fit$EL)
  expect_equal(fl2$fit$EL, fl$fit$EF)
})

test_that("very sparse loadings and factors are successfully added", {
  LL_sparse = c(NA, rep(0, 9))

  fl = flash(Y, fixed_loadings = LL_sparse, var_type = "constant")
  expect_true(all(fl$fit$EL[!is.na(LL_sparse), 1] == 0))
  expect_false(any(fl$fit$EL[is.na(LL_sparse), 1] == 0))

  fl2 = flash(t(Y), fixed_factors = LL_sparse, var_type = "constant")
  expect_equal(fl2$fit$EF, fl$fit$EL)
  expect_equal(fl2$fit$EL, fl$fit$EF)
})

test_that("a fixed loading is successfully added to an existing object", {
  # TODO
})



