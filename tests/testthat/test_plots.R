test_that("plotting functions work", {
  set.seed(1)
  l = cbind(c(1, 1, 1, 0, 0), c(0, 0, 1, 1, 1))
  f = cbind(rep(5, 20), 1:20)
  LF = l %*% t(f)
  Y = LF + rnorm(5*20)
  fl = flash(Y)

  pveplot = plot_pve(fl)
  expect_s3_class(pveplot, "ggplot")

  factorplot = plot_kset(fl, 1, factors=TRUE)
  expect_s3_class(factorplot, "ggplot")

  # test options:
  factorplot = plot_kset(fl, 1, factors=FALSE,
                         bar_colors=rep("blue", 5),
                         legend_size = 10,
                         plot_grid_ncol = 1)
  expect_s3_class(factorplot, "ggplot")
})
