# flashr2
repo for refactoring flashr code

A minimal example - code is currently very much in development and untested! Do not use for applications!

```
set.seed(1)
l = rnorm(5)
f = rnorm(20)
LF = outer(l,f)
Y = LF + rnorm(5*20)

data = set_flash_data(Y)
f1 = flash_r1(data,"svd")
fg = flash_greedy(data,10)
fb = flash_backfit(data,fg)
```
