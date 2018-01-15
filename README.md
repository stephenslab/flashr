# flashr

This code is in development, and the interface may change
under you! 

For a minimal example, install the `devtools` package and then:
```
devtools::install_github("stephenslab/flashr")
library("flashr")
set.seed(1)
l = rnorm(5)
f = rnorm(20)
LF = outer(l,f)
Y = LF + rnorm(5*20)

data = flash_set_data(Y)
f1 = flash_r1(data)
fg = flash_add_greedy(data,10)
fb = flash_backfit(data,fg)
flash_get_l(fb)
flash_get_l(fg)
```

If you want to use the faster `ebnm_pn` function (instead of the default adaptive shrinkage approach) you will have to install 
[*ebnm*](https://github.com/stephenslab/ebnm/)

For more see 
(https://github.com/stephenslab/flashr/blob/master/vignettes/flash_intro.Rmd)
