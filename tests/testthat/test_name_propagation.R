context("Propagation of column and row names")

set.seed(1)

LF = matrix(5, nrow = 5, ncol = 20)
Y = LF + rnorm(5 * 20)

l_names = c("a", "b", "c", "d", "e")
rownames(Y) = l_names
f_names = as.character(1:20)
colnames(Y) = f_names

fl = flash(Y, greedy_Kmax = 1)

test_that("row and column names are correctly propagated from data", {
  expect_identical(rownames(fl$ldf$l), l_names)
  expect_identical(rownames(fl$ldf$f), f_names)
  expect_identical(rownames(fl$fitted_values), l_names)
  expect_identical(colnames(fl$fitted_values), f_names)
})

names(Y) = NULL
l = 1:5
names(l) = l_names
f = matrix(1:20, nrow = 20, ncol = 1)
rownames(f) = f_names

fl = flash(Y, fixed_loadings = l, fixed_factors = f)

test_that("names are correctly propagated from fixed factors and loadings", {
  expect_identical(rownames(fl$ldf$l), l_names)
  expect_identical(rownames(fl$ldf$f), f_names)
  expect_identical(rownames(fl$fitted_values), l_names)
  expect_identical(colnames(fl$fitted_values), f_names)
})

fl = flash(Y, f_init = fl, greedy_Kmax = 1)

test_that("adding additional factors does not remove names", {
  expect_identical(rownames(fl$ldf$l), l_names)
  expect_identical(rownames(fl$ldf$f), f_names)
  expect_identical(rownames(fl$fitted_values), l_names)
  expect_identical(colnames(fl$fitted_values), f_names)
})
