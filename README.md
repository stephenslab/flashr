# flashr

This code is in development.  The interface is fairly stable but not guaranteed not to change.

For a minimal example, install the `devtools` package and then:
```
devtools::install_github("stephenslab/ebnm") # installs ebnm package
devtools::install_github("stephenslab/flashr") # installs flashr package
library("flashr")

# Set up some simulated data
set.seed(1) # for reproducibility
ftrue = matrix(rnorm(200),ncol=2)
ltrue = matrix(rnorm(40),ncol=2)
ltrue[1:10,1] = 0 # set up some sparsity
ltrue[11:20,2] = 0
Y = ltrue %*% t(ftrue)+rnorm(2000) # set up a simulated matrix

# Run flash

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



For more see 
(https://github.com/stephenslab/flashr/blob/master/vignettes/flash_intro.Rmd)
