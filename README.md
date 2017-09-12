# flashr2
repo for refactoring flashr code

A minimal example - code is currently very much untested!
```
set.seed(1)
l = rnorm(5)
f = rnorm(20)
LF = outer(l,f)
Y = LF + rnorm(5*20)
f1 = flash_r1(Y,"svd")
fg = flash_greedy(Y,10)
fb = flash_backfit(Y,fg)
```
