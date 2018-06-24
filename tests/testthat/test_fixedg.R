test_that("fixing gl and gf works as expected", {
  set.seed(2)
  LF = 3 * outer(rnorm(5),rnorm(20))
  Y = LF + rnorm(5*20)

  gl = ashr::normalmix(c(0.5, 0.5), c(0, 0), c(1, 2))
  data = flash_set_data(Y)
  fl = flash_add_greedy(data, 1, ebnm_fn=ebnm_ash, gl=gl, fixgl=TRUE)
  expect_identical(flash_get_gl(fl)[[1]], gl)

  # This test takes a few seconds...
  # gf = list()
  # gf[[1]] = ashr::normalmix(1, 0, 1)
  # gf[[2]] = ashr::normalmix(1, 0, 10)
  # fl = flash_add_greedy(data, 2, ebnm_fn=ebnm_ash, gf=gf, nullcheck=F)
  # expect_identical(flash_get_gf(fl), gf)

  gf = list()
  gf[[1]] = ashr::normalmix(1, 0, 0.25)
  fl = flash_add_greedy(data, 1, ebnm_fn=ebnm_ash)
  fl2 = flash_backfit(data, fl, 1, ebnm_fn=ebnm_ash, gf=gf, fixgf=TRUE)
  expect_identical(flash_get_gf(fl2), gf)
})
