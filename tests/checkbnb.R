q <- sample(0:10,1)
k <- runif(1,0.1,10)
a <- runif(1,0.1,10)
b <- runif(1,0.1,10)

library('brr')
p1 <- pbeta_nbinom(q,k,a,b)

library('bayescount')
p2 <- bayescount:::pbnb1(q,k,a,b)
p3 <- bayescount:::pbnb2(q,k,a,b)


abs(p1-p2)
abs(p1-p3)
