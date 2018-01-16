test_that("adding sparse factors/loadings works", {
  set.seed(1)
  n = 20
  p = 50

  n_ll = 1
  ll1 = matrix(0, ncol=n_ll, nrow=n)
  ll1[1:5, 1] <- c(1, 2, 3, 4, 5)
  ff1 <- matrix(rnorm(n_ll*p), ncol=n_ll, nrow=p)
  Y <- ll1 %*% t(ff1) + rnorm(n*p)

  n_ff = 2
  ff2 = matrix(0, ncol=n_ff, nrow=p)
  ff2[1:5, 1] <- c(1, 2, 3, 4, 5)
  ff2[6:10, 2] <- rep(10, 5)
  ll2 <- matrix(rnorm(n_ff*n), ncol=n_ff, nrow=n)
  Y <- ll2 %*% t(ff2) + rnorm(n*p)

  # Test sparse loadings:
  fl <- flash_add_sparse_l(Y, !(ll1 == 0))
  # Check that sparsity pattern is maintained and fixl is set correctly:
  expect_equal(sum(fl$EL[ll1 == 0] != 0), 0)
  expect_equal(sum(fl$EL[ll1 != 0] == 0), 0)
  expect_equal(sum(fl$fixl[ll1 == 0] == FALSE), 0)
  expect_equal(sum(fl$fixl[ll1 != 0] == TRUE), 0)
  # Test fix_nulls parameter:
  fl <- flash_add_sparse_l(Y, !(ll1 == 0), fix_nulls=F)
  expect_equal(sum(fl$fixl == TRUE), 0)

  # Test sparse factors (add to existing matrix):
  fl <- flash_add_sparse_f(Y, !(ff2 == 0), f_init=fl)
  # Check that sparsity pattern is maintained and fixf is set correctly:
  expect_equal(sum(fl$EF[,(n_ll+1):ncol(fl$EF)][ff2 == 0] != 0), 0)
  expect_equal(sum(fl$EF[,(n_ll+1):ncol(fl$EF)][ff2 != 0] == 0), 0)
  expect_equal(sum(fl$fixf[,(n_ll+1):ncol(fl$fixf)][ff2 == 0] == FALSE), 0)
  expect_equal(sum(fl$fixf[,(n_ll+1):ncol(fl$fixf)][ff2 != 0] == TRUE), 0)
  # Test fix_nulls parameter (initializing with null flash object):
  fl <- flash_add_sparse_f(Y, !(ff2 == 0), fix_nulls=F)
  expect_equal(sum(fl$fixf == TRUE), 0)
})
