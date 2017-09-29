# flashr2
repo for refactoring flashr code

This code is currently very much in development, and the interface may change
under you! 

For a minimal example
```
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

For more see 
(https://github.com/stephenslab/flashr2/blob/master/vignettes/flash_intro.Rmd)
