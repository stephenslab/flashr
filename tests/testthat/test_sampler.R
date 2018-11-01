context("sampler")

test_that(paste("sampling functions produce objects of correct",
                "dimensions and give reasonable results"), {
  set.seed(1)

  for (ebnm_fn in c("ebnm_ash", "ebnm_pn")) {

    # No fixed elements.
    l = seq(0.5, 10, by=0.5)
    f = rep(1, 30)
    Y = outer(l, f)
    l = rep(2, 20)
    f = seq(1.1, 4, by=0.1)
    Y = Y + outer(l, f)
    Y = Y + rnorm(20*30)
    fl = flash(Y, greedy_Kmax = 2, var_type = "constant",
               control = list(ebnm_fn = ebnm_fn))

    data = flash_set_data(Y)
    fit = fl$fit
    # Check flash_l_sampler.

    lsampler = flash_l_sampler(data, fl$fit, 1:flash_get_k(fit))
    lsamp = lsampler(50)
    expect_length(lsamp, 50)
    expect_equal(dim(lsamp[[1]]), c(20, 2))
    var_l = fit$EL2 - fit$EL^2
    sampled_means_l = Reduce(`+`, lsamp) / 50
    expect_equal(sampled_means_l, fit$EL, tolerance = sqrt(var_l))

    # Test LF sampler.
    lfsampler = flash_sampler(data, fit, fixed="factors")
    lfsamp = lfsampler(10)
    expect_length(lfsamp, 10)
    expect_equal(dim(lfsamp[[1]]), c(20, 30))

    lf_means = Reduce(`+`, lfsamp) / 10

    # Need to suppress warning about scale parameter being a matrix.
    LF = flash_get_fitted_values(fit)
    suppressWarnings(expect_equal(lf_means, LF, tolerance=0.1, scale=LF))

    # Fix some elements.
    fl = flash(Y, fixed_factors = list(vals = rep(1, 30),
                                        is_fixed = c(rep(TRUE, 10), rep(FALSE, 20))),
                var_type = "constant")
    fit = fl$fit

    # Check flash_f_sampler.
    fsampler = flash_f_sampler(data, fit, 1:flash_get_k(fit))
    fsamp = fsampler(50)
    expect_length(fsamp, 50)
    expect_length(fsamp[[1]], 30)

    # Make sure the fixed elements stay fixed.
    expect_equal(fsamp[[1]][1:10], rep(1, 10))
    expect_false(fsamp[[1]][11] == 1)
    var_f = fit$EF2 - fit$EF^2
    sampled_means_f = Reduce(`+`, fsamp) / 50
    expect_equal(sampled_means_f[1:10], fit$EF[1:10])
    expect_equal(sampled_means_f[11:30], fit$EF[11:30],
                 tolerance = sqrt(var_f[11:30]))

    flsampler = flash_sampler(data, fit, fixed="l")
    flsamp = flsampler(10)
    expect_length(flsamp, 10)
    expect_equal(dim(flsamp[[1]]), c(20, 30))
  }

  # Test fixed = "none" (and sampled variances):
  Y = matrix(5, nrow=5, ncol=20) + rnorm(100)
  fl = flash(Y, greedy_Kmax = 1,
             control = list(ebnm_fn = "ebnm_pn",
                                 ebnm_param = list(g = list(pi0 = 0, a = 1))))
  # Check posterior distribution of LF'_{1, 1}
  pmean_l = fl$fit$EL[1, 1]
  pvar_l = fl$fit$EL2[1, 1] - fl$fit$EL[1, 1]^2
  pmean_f = fl$fit$EF[1, 1]
  pvar_f = fl$fit$EF2[1, 1] - fl$fit$EF[1, 1]^2

  flsampler_l = flash_sampler(Y, fl, fixed="l")
  flsamp_l = flsampler_l(1000)
  samp_l = sapply(flsamp_l, function(x) x[1,1])
  var_l = var(samp_l)
  expect_equal(var_l, pmean_l^2 * pvar_f, tolerance = .01)

  flsampler_n = flash_sampler(Y, fl, fixed="none")
  flsamp_n = flsampler_n(1000)
  samp_n = sapply(flsamp_n, function(x) x[1,1])
  var_n = var(samp_n)
  expect_equal(mean(samp_n), pmean_l * pmean_f, tolerance = .01)
  expect_equal(var_n,
               pmean_l^2 * pvar_f + pmean_f^2 * pvar_l + pvar_l * pvar_f,
               tolerance = .01)

})
