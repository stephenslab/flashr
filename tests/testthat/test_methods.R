context("flash methods")

set.seed(1)

test_that("nonnegative method works as expected", {
  Y = matrix(5 + abs(rnorm(20)), nrow = 5, ncol = 10)

  fl_nn = flash(Y, greedy_Kmax = 1, method = "nonnegative",
                var_type = "constant")
  expect_true(all(fl_nn$ldf$l >= 0))
  expect_true(all(fl_nn$ldf$f >= 0))
})

test_that("nnloadings and nnfactors methods work as expected", {
  Y = outer(c(rep(-3, 2), rep(3, 3)), c(rep(-1, 5), rep(1, 5))) + rnorm(50)

  fl_nnl = flash(Y, greedy_Kmax = 1, method = "nnloadings",
                 var_type = "constant")
  expect_true(all(fl_nnl$ldf$l >= 0))
  expect_false(all(fl_nnl$ldf$f >= 0))

  fl_nnf = flash(Y, greedy_Kmax = 1, method = "nnfactors",
                 var_type = "constant")
  expect_true(all(fl_nnf$ldf$f >= 0))
  expect_false(all(fl_nnf$ldf$l >= 0))
})
