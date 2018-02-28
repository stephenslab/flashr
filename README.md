# flashr

Methods for matrix factorization based on
[Empirical Bayes Matrix Factorization](https://arxiv.org/abs/1802.06931).
The name of the package, "flashr," comes from "Factors and Loadings by
Adaptive SHrinkage in R".
	
*Note:* This code is in development. The interface is fairly stable
but not guaranteed to stay the same.

## License

Copyright (c) 2017-2018, Matthew Stephens and Wei Wang.

All source code and software in this repository are made available
under the terms of the [BSD 3-Clause
License](https://opensource.org/licenses/BSD-3-Clause). See
the [LICENSE](LICENSE) file for the full text of the license.

## Citing this work

If you find that this R package is useful for your work, please cite
our paper:

> W. Wang and M. Stephens, 2018. *Empirical Bayes matrix factorization.* 
[arXiv:1802.06931](https://arxiv.org/abs/1802.06931).

## Quick start

Follow these steps to quickly get started using `flashr`.

1. 

For a minimal example, install the `devtools` package and then:

```R
devtools::install_github("stephenslab/ebnm") # installs ebnm package
devtools::install_github("stephenslab/flashr") # installs flashr package
library("flashr")

## Set up some simulated data

set.seed(1) # for reproducibility
ftrue = matrix(rnorm(200),ncol=2)
ltrue = matrix(rnorm(40),ncol=2)
ltrue[1:10,1] = 0 # set up some sparsity
ltrue[11:20,2] = 0
Y = ltrue %*% t(ftrue)+rnorm(2000) # set up a simulated matrix

## Run flash

f = flash(Y)
ldf = flash_get_ldf(f)$d # show the weights, analogous to singular values showing importance of each factor
plot(ltrue[,1],ldf$l[,1]) # plot true l against estimated l (note estimate is normalized);
plot(ftrue,ldf$f) # plot true f against estimated f (note estimate is normalized)
plot(ltrue %*% t(ftrue), flash_get_lf(f)) #plot true lf' against estimated lf'; the scale of the estimate matches the data

# example to use the more flexible ebnm function in ashr; show how to pass parameters to
f2 = flash_r1(Y,ebnm_fn = ebnm_ash)
# example to show how to pass parameters to ashr
f3= flash_r1(Y,ebnm_fn = ebnm_ash, ebnm_param = list(mixcompdist = "normal",method="fdr"))

```

For more see [flash_intro.Rmd](vignettes/flash_intro.Rmd).

## How to build static HTML documentation

These are the R commands to build the website (make sure you are
connected to Internet while running these commands):

```R
library(pkgdown)
build_site(mathjax = FALSE)
```

## Credits

This software was developed by
[Matthew Stephens](http://stephenslab.uchicago.edu)
[Wei Wang](https://github.com/NKweiwang) and
[Jason Willwerscheid](https://github.com/willwerscheid) at the University
of Chicago.
