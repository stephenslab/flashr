test_that("summary produces expected output", {
  l = rep(1, 5)
  f = rnorm(20)
  LF = outer(l, f)
  Y = LF + rnorm(5*20)
  fl = flash(Y)

  smry = summary(fl)
  expect_s3_class(smry, "summary.flash")
  smry_elems = c("nfactors", "pve", "fitted.values", "ldf", "objective")
  expect_named(smry, smry_elems)
  expect_equal(smry$objective, NA)

  smry2 = summary(fl, Y)
  expect_type(smry2$objective, "double")

  #expect_output(print(fl), "Summary of flash object")
})
