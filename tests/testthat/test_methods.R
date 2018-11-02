context("flash methods")

set.seed(3)

Y = outer(c(rep(-1, 2), rep(3, 3)), c(rep(-1, 5), rep(3, 5))) + rnorm(50)

test_that("fastest and normalmix methods work as expected", {
  fl_pn = flash(Y, greedy_Kmax = 1, method = "fastest",
                var_type = "constant")
  expect_equal(fl_pn$fit$ebnm_fn_l[[1]], "ebnm_pn")
  expect_equal(fl_pn$fit$ebnm_fn_f[[1]], "ebnm_pn")

  fl_ash = flash(Y, greedy_Kmax = 1, method = "normalmix",
                 var_type = "constant")
  expect_equal(fl_ash$fit$ebnm_fn_l[[1]], "ebnm_ash")
  expect_equal(fl_ash$fit$ebnm_fn_f[[1]], "ebnm_ash")
})

test_that("nnloadings and nnfactors methods work as expected", {
  fl_nnl = flash(Y, greedy_Kmax = 1, method = "nnloadings",
                 var_type = "constant")
  expect_true(all(fl_nnl$ldf$l >= 0))
  expect_false(all(fl_nnl$ldf$f >= 0))

  fl_nnf = flash(Y, greedy_Kmax = 1, method = "nnfactors",
                 var_type = "constant")
  expect_true(all(fl_nnf$ldf$f >= 0))
  expect_false(all(fl_nnf$ldf$l >= 0))
})

test_that("nonnegative method works as expected", {
  Y = matrix(5 + abs(rnorm(20)), nrow = 5, ncol = 10)

  fl_nn = flash(Y, greedy_Kmax = 1, method = "nonnegative",
                var_type = "constant")
  expect_true(all(fl_nn$ldf$l >= 0))
  expect_true(all(fl_nn$ldf$f >= 0))
})
