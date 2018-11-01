context("flash object")

set.seed(666)
n = 10
p = 50
k = 5

LL = matrix(rnorm(n * k), nrow=n, ncol=k)
FF = matrix(0, nrow=p, ncol=k)
for (k in 1:5) {
  FF[(1 + 5 * (k - 1)):(5 * (k + 1)), k] = 1
}
Y = LL %*% t(FF) + rnorm(n * p)

fo = flash(Y, greedy_Kmax = 2, var_type = "constant")

test_that("flash object classes are set correctly", {
  expect_s3_class(fo, "flash")
  expect_s3_class(fo$fit, "flash_fit")
})

test_that("flash object fields are set as expected", {
  expect_equal(length(fo$fit_history), 2)

  # Backfit the two factors:
  fo2 = flash(Y, f_init = fo, backfit = TRUE, var_type = "constant")
  expect_equal(length(fo2$fit_history), 3)

  # Add a single fixed loading:
  fo3 = flash(Y, f_init = fo2, fixed_loadings = rep(1, n),
              var_type = "constant")
  expect_equal(fo3$nfactors, 3)
  expect_equal(length(fo3$fit_history), 4)

  # Zero out the first factor:
  fo4 = flashr:::flash_zero_out_factor(flash_set_data(Y), f_init = fo3,
                                       k = 1)
  expect_equal(fo$nfactors, 2)
  expect_false(is.na(fo$objective))

  # Add two factors without optimizing and backfit together:
  fo5 = flash(Y, f_init = fo4, greedy_Kmax = 2, backfit = TRUE,
              var_type = "constant",
              control = list(r1opt_maxiter = 0))
  expect_equal(length(fo5$fit_history), 8) # nullcheck causes second backfit

  # Attempt to add another factor (it will be zeroed out):
  fo6 = flash(Y, f_init = fo5, greedy_Kmax = 1, var_type = "constant")
  expect_equal(length(fo6$fit_history), 9)
  expect_equal(fo6$fit_history[[9]]$zeroed_out, fo6$fit_history[[9]]$kset)
})

test_that("flash object objective is set to NA when appropriate", {
  fo2 = flash(Y, fixed_loadings = rep(1, 10), f_init = fo,
              control = list(r1opt_maxiter = 0))
  expect_true(is.na(fo2$objective))

  fo3 = flashr:::flash_zero_out_factor(flash_set_data(Y), f_init = fo2, k = 3)
  expect_true(is.na(fo3$objective))
})
