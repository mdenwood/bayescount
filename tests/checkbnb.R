library('bayescount')

if(FALSE){

pre <- rnbinom(100, 1, mu=100)
post <- rnbinom(100, 1, mu=20)
pre <- rnbinom(100, 1, mu=10)
post <- rnbinom(100, 1, mu=1)
bayescount:::fecrt_analyses(pre, post)


post <- round(pre * 0.1)
cov(pre,post)
bayescount:::fecrt_analyses(pre, post)


## Outstanding checks:

# C++ BNB, WAAVP, Dobson and Assymptotic methods for fixed dataset

# fecrt_power vs R script generating the same data (first copy power_comparison code to power)

# fecrt_power_comparison


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

}
