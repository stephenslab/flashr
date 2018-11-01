context("flash sanity checks")

set.seed(666)
l = rnorm(5)
f = rnorm(20)
LF = outer(l, f)
Y = LF + rnorm(5 * 20)

# Add a single factor/loading pair using udv_random. Don't optimize.
fl = flash(Y,
           method = "fastest",
           greedy_Kmax = 1,
           nullcheck = FALSE,
           verbose = FALSE,
           control = list(init_fn = "udv_random", r1opt_maxiter = 0))
mse1 = mean((LF - fl$fitted_values)^2)

# Do a single backfitting iteration.
fl2 = flash(Y,
            f_init = fl,
            method = "fastest",
            backfit = TRUE,
            nullcheck = FALSE,
            verbose = FALSE,
            control = list(backfit_maxiter = 1))
mse2 = mean((LF - fl2$fitted_values)^2)

test_that("flash update improves mean squared error in simple situation", {
  expect_true(mse1 > mse2)
})
