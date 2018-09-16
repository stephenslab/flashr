test_that("storing tau as a scalar works", {
  set.seed(1)
  l = rnorm(5)
  f = rnorm(20)
  LF = outer(l,f)
  Y = LF + rnorm(5*20)

  # test that objective is same whether S is matrix or scalar:
  data_scalar = flash_set_data(Y, S = 1)
  fl_scalar = flash(data_scalar, var_type="zero")
  data_matrix = flash_set_data(Y, S = matrix(1, nrow=5, ncol=20))
  fl_matrix = flash(data_matrix, var_type="zero")

  fl1 = flash(Y, var_type="constant")
  expect_length(fl1$fit$tau, 1)

  # Test that calc_ebnm_args gives same results whether we pass in a
  #   scalar or matrix with var_type = "constant":
  largs1 = calc_ebnm_l_args(data, fl1$fit, 1, 1:5,
                            flash_get_Rk(data, fl1$fit, 1))
  expect_length(largs1$s, 1)
  fl2 = fl1
  fl2$fit$tau = matrix(fl1$fit$tau, nrow=5, ncol=20)
  largs2 = calc_ebnm_l_args(data, fl2$fit, 1, 1:5,
                            flash_get_Rk(data, fl2$fit, 1))
  expect_equal(largs1$x, largs2$x)
  expect_equal(largs1$s, largs2$s[1])

  fargs1 = calc_ebnm_f_args(data, fl1$fit, 1, 1:20,
                            flash_get_Rk(data, fl1$fit, 1))
  expect_length(fargs1$s, 1)
  fargs2 = calc_ebnm_f_args(data, fl2$fit, 1, 1:20,
                            flash_get_Rk(data, fl2$fit, 1))
  expect_equal(fargs1$x, fargs2$x)
  expect_equal(fargs1$s, fargs2$s[1])

  expect_equal(flash_get_pve(data, fl1),
               flash_get_pve(data, fl2))

  # Repeat with missing data:
  Y[sample(1:100, 10, replace=FALSE)] = NA
  data = flash_set_data(Y)

  largs1 = calc_ebnm_l_args(data, fl1$fit, 1, 1:5,
                            flash_get_Rk(data, fl1$fit, 1))
  expect_length(largs1$s, 5) # SEs differ where there is missing data
  fl2 = fl1
  fl2$fit$tau = matrix(fl1$fit$tau, nrow=5, ncol=20)
  largs2 = calc_ebnm_l_args(data, fl2$fit, 1, 1:5,
                            flash_get_Rk(data, fl2$fit, 1))
  expect_equal(largs1$x, largs2$x)
  expect_equal(largs1$s, largs2$s)

  fargs1 = calc_ebnm_f_args(data, fl1$fit, 1, 1:20,
                            flash_get_Rk(data, fl1$fit, 1))
  expect_length(fargs1$s, 20)
  fargs2 = calc_ebnm_f_args(data, fl2$fit, 1, 1:20,
                            flash_get_Rk(data, fl2$fit, 1))
  expect_equal(fargs1$x, fargs2$x)
  expect_equal(fargs1$s, fargs2$s)
})
