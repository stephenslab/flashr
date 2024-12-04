# flashr: Empirical Bayes Factor Analysis in R

[![Travis Build Status](https://travis-ci.org/stephenslab/flashr.svg?branch=master)](https://travis-ci.org/stephenslab/flashr) 
[![Appveyor Build status](https://ci.appveyor.com/api/projects/status/ge410qhikk6j8iks?svg=true)](https://ci.appveyor.com/project/pcarbo/flashr)
[![codecov](https://codecov.io/gh/stephenslab/flashr/branch/master/graph/badge.svg)](https://codecov.io/gh/stephenslab/flashr)

Methods for matrix factorization based on
[Empirical Bayes Matrix Factorization](https://arxiv.org/abs/1802.06931).
The name of the package, "flashr," comes from "Factors and Loadings by
Adaptive SHrinkage in R".

## Important note

We are no longer supporting or developing flashr: 
we have switched to supporting a replacement R package,
[flashier](https://github.com/willwerscheid/flashier), which contains
most of the features of flashr, plus many important improvements,
including a new, user-friendly interface, more modeling options,
faster algorithms, and more detailed statistical outputs. **Please use
flashier rather than flashr**!  

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

1. Install the flashr using
   [devtools](https://github.com/r-lib/devtools). Please note that it
   can take several minutes to build the vignettes:

   ```R
   install.packages("devtools")
   library(devtools)
   install_github("stephenslab/flashr",build_vignettes = TRUE)
   ```

   This command should automatically retrieve and install the `ashr`
   and `ebnm` packages from GitHub (and possibly other packages). If
   it does not, install `ashr` and `ebnm` separately using devtools:

   ```R
   install_github("stephens999/ashr")
   install_github("stephenslab/ebnm")
   ```

   *Note:* If you are interested in attempting to reproduce the
   results in the Wang and Stephens (2018) manuscript, the flashr
   release that most closely matches the package used in the paper is
   [version
   0.4-10](https://github.com/stephenslab/flashr/releases/tag/v0.4-10).
   This release can be installed by running the following in R:

   ```R
   install_github("stephenslab/flashr@v0.4-10")
   ```
   
2. Optionally, install MOSEK and the Rmosek package, for faster model
   fitting. See the
   [ashr GitHub repository](https://github.com/stephens999/ashr) for
   details.

3. Run a few toy examples illustrating the `flash` function:

   ```R
   example("flash")
   ```

4. Explore the [introductory flashr
   vignette](https://stephenslab.github.io/flashr/articles/flash_intro.html):

   ```R
   vignette("flash_intro")
   ```

5. Explore the
   [vignette illustrating some more advanced features of
   flashr](https://stephenslab.github.io/flashr/articles/flash_advanced.html):

   ```R
   vignette("flash_advanced")
   ```

6. See the [online documentation](https://stephenslab.github.io/flashr) 
   to learn more about the `flashr` package.

## Developer notes

+ Run `pkgdown::build_site(mathjax = FALSE)` in R to build the website
using [pkgdown](https://github.com/r-lib/pkgdown). Make sure you are
connected to the Internet while running these commands.

## Credits

This software was developed by
[Matthew Stephens](http://stephenslab.uchicago.edu),
[Wei Wang](https://github.com/NKweiwang),
[Jason Willwerscheid](https://github.com/willwerscheid) and
[Peter Carbonetto](http://pcarbo.github.io) at the University
of Chicago.
