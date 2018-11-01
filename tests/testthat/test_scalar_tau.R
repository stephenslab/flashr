context("storing tau as a scalar")

set.seed(1)
l = rnorm(5)
f = rnorm(20)
LF = outer(l,f)
Y = LF + rnorm(5*20)

test_that("the objective is the same whether S is a matrix or a scalar", {
  fl_scalar = flash(Y, S = 1, var_type = "zero", greedy_Kmax = 1)
  fl_matrix = flash(Y, S = matrix(1, nrow=5, ncol=20), var_type = "zero",
                    greedy_Kmax = 1)
  expect_equal(fl_scalar$objective, fl_matrix$objective, tol = 1e-4)
})

fl1 = flash(Y, var_type = "constant", greedy_Kmax = 1)

test_that("tau is indeed stored as a scalar when var_type = constant", {
  expect_length(fl1$fit$tau, 1)
})

test_that(paste("calc_ebnm_args gives same results whether we pass in a",
                "scalar or matrix with var_type = constant"), {
  data_scalar = flash_set_data(Y, S = 1)
  largs1 = calc_ebnm_l_args(data_scalar, fl1$fit, 1, 1:5, FALSE,
                            flash_get_Rk(data_scalar, fl1$fit, 1))
  expect_length(largs1$s, 1)
  fl2 = fl1
  fl2$fit$tau = matrix(fl1$fit$tau, nrow=5, ncol=20)
  data_matrix = flash_set_data(Y, S = matrix(1, nrow=5, ncol=20))
  largs2 = calc_ebnm_l_args(data_matrix, fl2$fit, 1, 1:5, FALSE,
                            flash_get_Rk(data_matrix, fl2$fit, 1))
  expect_equal(largs1$x, largs2$x)
  expect_equal(largs1$s, largs2$s[1])

  fargs1 = calc_ebnm_f_args(data_scalar, fl1$fit, 1, 1:20, FALSE,
                            flash_get_Rk(data_scalar, fl1$fit, 1))
  expect_length(fargs1$s, 1)
  fargs2 = calc_ebnm_f_args(data_matrix, fl2$fit, 1, 1:20, FALSE,
                            flash_get_Rk(data_matrix, fl2$fit, 1))
  expect_equal(fargs1$x, fargs2$x)
  expect_equal(fargs1$s, fargs2$s[1])

  expect_equal(flash_get_pve(data_scalar, fl1),
               flash_get_pve(data.matrix, fl2))
})

test_that("storing tau as a scalar works with missing data", {
  Y[sample(1:100, 10, replace=FALSE)] = NA
  fl_scalar = flash(Y, S = 1, greedy_Kmax = 1, var_type = "zero")
  fl_matrix = flash(Y, S = matrix(1, nrow=5, ncol=20), greedy_Kmax = 1,
                    var_type="zero")
  expect_equal(fl_scalar$objective, fl_matrix$objective, tol = 1e-4)

  data_scalar = flash_set_data(Y, S = 1)
  largs1 = calc_ebnm_l_args(data_scalar, fl1$fit, 1, 1:5, FALSE,
                            flash_get_Rk(data_scalar, fl1$fit, 1))
  expect_length(largs1$s, 5) # SEs differ where there is missing data
  fl2 = fl1
  fl2$fit$tau = matrix(fl1$fit$tau, nrow=5, ncol=20)
  data_matrix = flash_set_data(Y, S = matrix(1, nrow=5, ncol=20))
  largs2 = calc_ebnm_l_args(data_matrix, fl2$fit, 1, 1:5, FALSE,
                            flash_get_Rk(data_matrix, fl2$fit, 1))
  expect_equal(largs1$x, largs2$x)
  expect_equal(largs1$s, largs2$s)

  fargs1 = calc_ebnm_f_args(data_scalar, fl1$fit, 1, 1:20, FALSE,
                            flash_get_Rk(data_scalar, fl1$fit, 1))
  expect_length(fargs1$s, 20)
  fargs2 = calc_ebnm_f_args(data_matrix, fl2$fit, 1, 1:20, FALSE,
                            flash_get_Rk(data_matrix, fl2$fit, 1))
  expect_equal(fargs1$x, fargs2$x)
  expect_equal(fargs1$s, fargs2$s)
})
