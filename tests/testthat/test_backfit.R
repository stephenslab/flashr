context("backfit")

set.seed(1)

LL = cbind(c(rep(3, 3), rep(0, 2)), c(rep(0, 2), rep(5, 3)))
FF = matrix(rnorm(20), ncol = 2)
Y = LL %*% t(FF) + rnorm(50)

fl = flash(Y, greedy_Kmax = 2, nullcheck = FALSE, var_type = "constant")

test_that("backfit = TRUE backfits the entire flash object", {
  fl_b = flash(Y, f_init = fl, backfit = TRUE, var_type = "constant")
  expect_false(any(fl$ldf$l == fl_b$ldf$l))
})

test_that("setting backfit to a vector works as expected", {
  fl_b = flash(Y, f_init = fl, backfit = 1:2, var_type = "constant")
  expect_false(any(fl$ldf$l == fl_b$ldf$l))

  fl_b = flash(Y, f_init = fl, backfit = 1, var_type = "constant")
  expect_false(any(fl$ldf$l[, 1] == fl_b$ldf$l[, 1]))
  expect_identical(fl$ldf$l[, 2], fl_b$ldf$l[, 2])
})
