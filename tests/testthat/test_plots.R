test_that("plotting functions work", {
  set.seed(1)
  l = cbind(c(1, 1, 1, 0, 0), c(0, 0, 1, 1, 1))
  f = cbind(rep(5, 20), 1:20)
  LF = l %*% t(f)
  Y = LF + rnorm(5*20)
  fl = flash(Y)

  pveplot = flash_plot_pve(fl)
  expect_s3_class(pveplot, "ggplot")

  factorplots = flash_plot_factors(data, fl)
  expect_named(factorplots, c("plot_f", "plot_l"))
  expect_s3_class(factorplots$plot_f[[1]], "ggplot")
})
