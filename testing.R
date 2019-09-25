df <- fecrt_sim_unpaired(as.integer(100), as.double(c(0.1, 0.2)), as.integer(10), as.integer(10), as.double(20), as.double(1.0), as.double(0.8), matrix(as.double(c(0.9,0.95,0.8,0.95)), ncol=2, byrow=TRUE), as.double(c(1,1)), as.integer(1), as.integer(10000), as.integer(1), as.double(0.025))
df
summary(df)
table(df$Classification, useNA='always')


# Test c++ function:
estimate_k(as.double(10.0), as.double(5.0), as.double(10.0), as.double(5.0), as.double(2.0), as.logical(TRUE))

d1 <- rnbinom(10, 1, mu=10)
d2 <- rnbinom(10, 1, mu=5)
estimate_k(as.double(mean(d1)), as.double(var(d1)), as.double(mean(d2)), as.double(var(d2)), as.double(cov(d1,d2)), as.logical(TRUE))

estimate_k(as.double(mean(d1)), as.double(var(d1)), as.double(mean(d2)), as.double(var(d2)), as.double(-cov(d1,d2)), as.logical(TRUE))


# Recover true non-correlated variance:
v1 <- 12
v2 <- 10
m1 <- 200
m2 <- 100
vc <- 8
vcov <- matrix(c(v1,vc,vc,v2), ncol=2)
mdat <- MASS::mvrnorm(10^6, c(m1,m2), vcov)
vcv <- var(mdat)
mdat1 <- mdat[,1]
mdat2 <- mdat[,2]

cor(mdat[,1], mdat[,2])
cor(mdat[,1]*10, mdat[,2])
vcv

pdat1 <- rpois(nrow(mdat), mdat[,1])
pdat2 <- rpois(nrow(mdat), mdat[,2])
var(cbind(pdat1,pdat2))

cor(pdat1, pdat2)
cov(pdat1, pdat2)
cov(mdat1, mdat2)

(correctedcorr <- cov(pdat1,pdat2) / (sqrt(var(pdat1) - mean(pdat1)) * sqrt(var(pdat2) - mean(pdat2))))
cor(mdat1, mdat2)

cor(mdat1, mdat2)
cor(mdat1/mean(mdat1), mdat2/mean(mdat2))

cov(mdat1/mean(mdat1), mdat2/mean(mdat2)) / (sqrt(var(pdat1) - mean(pdat1))/mean(mdat1) * sqrt(var(pdat2) - mean(pdat2))/mean(mdat2))


cvar1 <- (var(pdat1) - mean(pdat1))
cvar2 <- (var(pdat2) - mean(pdat2))
(uncvar1 <- cvar1 - (correctedcorr * (sqrt(cvar1)*sqrt(cvar2))))
(v1-vc)
(uncvar2 <- cvar2 - (correctedcorr * (sqrt(cvar1)*sqrt(cvar2))))
(v2-vc)


dat_ctl <- pdat1
dat_tx <- pdat2
# This wasn't generated as gamma but this assumes the dist is gamma:
estimate_k(as.double(mean(dat_ctl)), as.double(var(dat_ctl)), as.double(mean(dat_tx)), as.double(var(dat_tx)), as.double(cov(dat_ctl/mean(dat_ctl), dat_tx/mean(dat_tx))))
# Still works comparing to the coefficient of variation though:
1 / (sd(mdat1)/mean(mdat1))^2
1 / (sd(mdat2)/mean(mdat2))^2


# TODO: work out kc/k1/k2 from real data and check sims here:
powersim_unpaired(rvs=1-seq(0.6,0.8,by=0.001), N=91, mu=74, k1=1.27, k2=1.57, ta=0.65, ti=0.7, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(0.6,0.8,by=0.001), N=91, mu=74, kc=1, k1=0.44, k2=0.54, ta=0.65, ti=0.7, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_unpaired(rvs=1-seq(0.4,0.6,by=0.001), N=91, mu=162, k1=2.56, k2=0.97, ta=0.45, ti=0.5, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(0.4,0.6,by=0.001), N=91, mu=162, kc=1.2, k1=0.82, k2=0.30, ta=0.45, ti=0.5, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_unpaired(rvs=1-seq(0.8,1,by=0.001), N=91, mu=1255, k1=0.18, k2=0.18, ta=0.9, ti=0.95, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(0.8,1,by=0.001), N=91, mu=1255, kc=0.5, k1=0.18, k2=0.18, ta=0.9, ti=0.95, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.8,1,by=0.001), N=50, mu=50, k1=0.18, k2=0.18, ta=0.9, ti=0.95, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))


powersim_unpaired(rvs=1-seq(0.6,0.8,by=0.001), N=91, mu=74, k1=1.27, k2=1.57, ta=0.65, ti=0.7, iters=10^4, approx=1, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.8,1,by=0.001), N=10, mu=20, k1=0.18*3, k2=0.18*3, ta=0.9, ti=0.95, iters=10^3, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.94,0.96,by=0.001), N=20, mu=20, k1=0.18*3, k2=0.18*3, ta=0.95, ti=0.95, iters=10^4, approx=0, priors=c(0.0001,0.0001), pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_unpaired(rvs=1-seq(0.94,0.96,by=0.001), N=20, mu=20, k1=0.18*3, k2=0.18*3, ta=0.95, ti=0.95, iters=10^4, approx=0, priors=c(0,10), pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.94,0.96,by=0.001), N=20, mu=20, k1=0.18*3, k2=0.18*3, ta=0.95, ti=0.95, iters=10^4, approx=0, priors=c(0,10), pm=c('BNB','Levecke','MLE','WAAVP'))


powersim_unpaired(rvs=1-seq(0.94,0.96,by=0.001), N=20, mu=20, k1=0.18*3, k2=0.18*3, ta=0.95, ti=0.95, iters=10^4, approx=2, pm=c('BNB','Levecke','MLE','WAAVP'))


powersim_unpaired(rvs=1-seq(0.48,0.52,by=0.001), N=30, mu=30, k1=0.18*3, k2=0.18*3, ta=0.5, ti=0.5, iters=10^4, approx=0, priors=c(0,0), pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.48,0.52,by=0.001), N=30, mu=30, k1=0.18*3, k2=0.18*3, ta=0.5, ti=0.5, iters=10^4, approx=2, pm=c('BNB','Levecke','MLE','WAAVP'))


powersim_unpaired(rvs=1-seq(-0.02,0.02,by=0.001), N=30, mu=30, k1=0.18*3, k2=0.18*3, ta=0, ti=0, iters=10^4, approx=0, priors=c(0.01,0.01), pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.48,0.52,by=0.001), N=30, mu=30, k1=0.18*3, k2=0.18*3, ta=0.5, ti=0.5, iters=10^4, approx=0, priors=c(0.01,0.01), pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.98,1,by=0.001), N=30, mu=30, k1=0.18*3, k2=0.18*3, ta=0.99, ti=0.99, iters=10^4, approx=0, priors=c(0.01,0.01), pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(0.98,1,by=0.001), N=10, mu=10, k1=0.18*3, k2=0.18*3, ta=0.99, ti=0.99, iters=10^4, approx=0, priors=c(0.01,0.01), pm=c('BNB','Levecke','MLE','WAAVP'))


findtheta(rnbinom(20,0.1,mu=20))
findtheta(c(0,0,0))
findtheta(c(1,10000,10^6))

d1 <- rnbinom(20,0.1,mu=20)
d2 <- rnbinom(20,1,mu=2)
estimate_k_ml(as.integer(d1),as.double(mean(d1)), as.double(var(d1)), as.integer(d2), as.double(mean(d2)), as.double(var(d2)), as.double(cov(d1,d2)), as.logical(TRUE))
estimate_k_ml(as.integer(d1),as.double(mean(d1)), as.double(var(d1)), as.integer(d2), as.double(mean(d2)), as.double(var(d2)), as.double(cov(d1,d2)), as.logical(FALSE))

powersim_unpaired(rvs=1-seq(0.48,0.52,by=0.001), N=10000, mu=10000, k1=1, k2=1, ta=0.5, ti=0.5, iters=10^3, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(0.48,0.52,by=0.001), N=1000, mu=1000, kc=1.5, k1=1, k2=1, ta=0.5, ti=0.5, iters=10^3, approx=1, pm=c('BNB','Levecke','MLE','WAAVP','Dobson'))


## THIS SHOWS A PROBLEM WITH ESTIMATING KC/K1/K2 FOR PAIRED:
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=30, mu=30, kc=1.5, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=30, mu=30, kc=11.6, k1=9.2, k2=5.6, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=30, mu=30, kc=0.5, k1=0.4, k2=0.4, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=30, mu=30, kc=5, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(0.48,0.52,by=0.001), N=30, mu=30, kc=1.5, k1=1, k2=1, ta=0.5, ti=0.5, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=100, mu=100, kc=5, k1=2, k2=2, ta=0, ti=0, iters=10^4, approx=1, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_unpaired(rvs=1-seq(-0.02,0.02,by=0.001), N=100, mu=100, k1=2, k2=2, ta=0, ti=0, iters=10^4, approx=1, pm=c('BNB','Levecke','MLE','WAAVP'))


powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=100, mu=100, kc=2, k1=1, k2=1, ta=0, ti=0, iters=10^3, approx=1, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=20, mu=10, kc=2, k1=0.1, k2=0.1, ta=0, ti=0, iters=10^3, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'), useml=TRUE)

powersim_unpaired(rvs=1-seq(-0.02,0.02,by=0.001), N=20, mu=10, k1=0.1, k2=0.1, ta=0, ti=0, iters=10^3, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'), useml=TRUE)



# Test correlations:
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=50, mu=50, kc=1.2, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=50, mu=50, kc=1.5, k1=1.2, k2=0.8, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=91, mu=50, kc=2, k1=1.2, k2=0.8, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=50, mu=50, kc=2, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_paired(rvs=1-seq(-0.02,0.02,by=0.001), N=50, mu=50, kc=200, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))
powersim_unpaired(rvs=1-seq(-0.02,0.02,by=0.001), N=50, mu=50, k1=1, k2=1, ta=0, ti=0, iters=10^4, approx=0, pm=c('BNB','Levecke','MLE','WAAVP'))

powersim_unpaired(rvs=1-seq(-0.02,0.02,by=0.001), N=10000, mu=1000, k1=1, k2=1, ta=0, ti=0, iters=10^3, approx=1, pm=c('BNB','Levecke','MLE','WAAVP'))


## Testing pbnb:

for(i in 0:10){
	noninc <- ifelse(i==0, 0, pbnb(i-1, 1, 20, 20, TRUE, TRUE))
	inc <- pbnb(i, 1, 20, 20, TRUE, TRUE)
	print(inc - noninc)
}
for(i in 0:10){
	noninc <- pbnb(i, 1, 20, 20, TRUE, FALSE)
	inc <- pbnb(i, 1, 20, 20, TRUE, TRUE)
	print(inc - noninc)
}

for(i in 0:10){
	inc <- 1 - ifelse(i==0, 0, pbnb(i-1, 1, 20, 20, TRUE, TRUE))
	noninc <- 1 - pbnb(i, 1, 20, 20, TRUE, TRUE)
	print(inc - noninc)
}
for(i in 0:10){
	noninc <- pbnb(i, 1, 20, 20, FALSE, FALSE)
	inc <- pbnb(i, 1, 20, 20, FALSE, TRUE)
	print(inc - noninc)
}

pbnb_upper(0, 1, 20, 20)
pbnb_lower(1, 1, 20, 20)

pbnb_upper(1, 1, 20, 20) - pbnb_lower(1, 1, 20, 20)
pbnb_upper(2, 1, 20, 20) - pbnb_lower(2, 1, 20, 20)

#powersim_unpaired(N=100)

powersim_unpaired(N=1000, approx=2)

powersim_unpaired(rvs=seq(0.25,0.75,by=0.001), ta=0.45, ti=0.5, N=500)

powersim_paired(kc=1.2, mu=50, k1=1, k2=1, N=1000)

powersim_paired(kc=0.8, mu=5, k1=0.5, k2=0.5, N=5, approx=0)
powersim_unpaired(mu=5, k1=0.5, k2=0.5, N=5, approx=0)

powersim_paired(kc=0.8, mu=5, k1=0.5, k2=0.5, N=5, approx=2)

# This is a nice illustration of continuous vs discrete:
powersim_unpaired(k1=1, k2=1, N=10, mu=10, approx=0)

powersim_unpaired(k1=1, k2=1, N=10, mu=10, approx=2)
powersim_paired(N=10, mu=10, approx=0)
powersim_paired(N=10, mu=10, approx=2)

powersim_unpaired(k1=1, k2=1, N=10, mu=50, approx=0)

powersim_unpaired(k1=1, k2=1, N=10, mu=50, approx=2)
powersim_paired(N=10, mu=50, approx=0)
powersim_paired(N=10, mu=50, approx=2)

powersim_unpaired(k1=1, k2=1, N=50, mu=50, approx=0)
powersim_unpaired(k1=1, k2=1, N=50, mu=50, approx=2)
powersim_paired(N=50, mu=50, approx=0)
powersim_paired(N=50, mu=50, approx=2)

powersim_unpaired(k1=1, k2=1, N=500, mu=500, approx=1)
powersim_paired(N=500, mu=500, approx=1)


powersim_unpaired(k1=1, k2=1, N=20, mu=50, approx=1, ta=0.95, ti=0.99)

# Performance is more similar with ta around 0.5:
powersim_unpaired(rvs=seq(0.25,0.75,by=0.01), N=50, mu=50, approx=1, ta=0.5, ti=0.55)
powersim_paired(rvs=seq(0.25,0.75,by=0.01), N=20, mu=50, approx=1, ta=0.5, ti=0.55)


powersim_unpaired(k1=2, k2=2, N=50, mu=20, approx=2, pm='BNB')
powersim_unpaired(k1=0.26, k2=0.26, N=50, mu=20, approx=0, pm='BNB')


## How to generate correlated gamma variates?

mu <- 200
red <- 0.1
N <- 10^4
k_ctl <- 0.8
k_tx <- 0.6
c <- 0.8
b1 <- c - k_ctl
b2 <- c - k_tx
gams <- rgamma(N, c, scale=1)
mval_ctl <- rbeta(N, k_ctl, b1) * gams * mu/k_ctl
mval_tx <- rbeta(N, k_tx, b2) * gams * (mu*red)/k_tx
mean(mval_ctl)^2/var(mval_ctl)
mean(mval_tx)^2/var(mval_tx)
dat_ctl <- rpois(N, mval_ctl)
dat_tx <- rpois(N, mval_tx)

mean(dat_ctl)^2 / (var(dat_ctl)-mean(dat_ctl))
mean(dat_tx)^2 / (var(dat_tx)-mean(dat_tx))
cor(dat_ctl, dat_tx)

estimate_k(as.double(mean(dat_ctl)), as.double(var(dat_ctl)), as.double(mean(dat_tx)), as.double(var(dat_tx)), as.double(cov(dat_ctl/mean(dat_ctl), dat_tx/mean(dat_tx))))

var_eff_ctl = var(dat_ctl) - mean(dat_ctl);
var_eff_tx = var(dat_tx) - mean(dat_tx);
scaled_cov = as.double(cov(dat_ctl/mean(dat_ctl), dat_tx/mean(dat_tx)))
sdprod = sqrt(var_eff_ctl) * sqrt(var_eff_tx)
c_cor <- (scaled_cov * mean(dat_ctl) * mean(dat_tx)) / sdprod
c_cor <- scaled_cov / sdprod

var_eff_ctl
var_eff_ctl * (1-cor(dat_ctl, dat_tx))
var_eff_tx
var_eff_tx * (1-cor(dat_ctl, dat_tx))


var(dat_tx)-mean(dat_tx)
estimate_k(as.double(mean(dat_ctl)), as.double(var(dat_ctl)), as.double(mean(dat_tx)), as.double(var(dat_tx)), as.double(cov(dat_ctl, dat_tx)))
k_ctl
k_tx


# Approximating Beta_NB with gamma:
r=k <- 125
a <- 167
b <- 436

it <- 10^6
# Swap successes and failures:
sim <- rnbinom(it, k, prob=rbeta(it, a, b))
mean(sim)
var(sim)

p = (a/(a+b))
bb <- a
aa <- b
(totalvar = ((r*(aa+r-1)*bb*(aa+bb-1))/((aa-2)*(aa-1)^2)))
var(sim)
(totalmean = (r*bb)/(aa-1))
mean(sim)

(scale = totalvar/totalmean)
(shape = totalmean/scale)

gsim <- rgamma(it, shape, scale=scale)

plot.ecdf(sim)
plot.ecdf(gsim, col='red', add=TRUE)

sum <- 54
pbnb(sum, k, aa, bb, TRUE)
pgamma(sum+0.5, shape, scale=scale, lower=TRUE)
pbnb(sum, k, aa, bb, FALSE)
pgamma(sum-0.5, shape, scale=scale, lower=FALSE)


# Using a Poisson approximation to the hypergeometric:
library('hypergeo')
hypergeo(r,a,a+b+r,exp(sum))



var(rbeta(10^4,a,b)*rgamma(10^4,c,1))
var(rgamma(10^4,a,1))


df <- fecrt_sim_unpaired(as.integer(1000), 0.049, as.integer(20), as.integer(20), as.double(20), as.double(1.0), as.double(0.8), matrix(as.double(c(0.9,0.95)), ncol=2, byrow=TRUE), as.double(c(1,1)), as.integer(1), as.integer(10000), as.integer(1), as.double(0.025))

summary(df$ObsReduction)

# If the reduction is identical in all animals, should we have a reduced CI as for the Levecke and WAAVP methods ...?  This does kind of make sense.
# 

pre <- rnbinom(10, 1, mu=100)
post <- round(pre*0.1)
methodcomp(as.integer(sum(pre)), as.integer(length(pre)), as.double(1.0), as.double(mean(pre)), as.double(var(pre)), as.integer(sum(post)), as.integer(length(post)), as.double(1.0), as.double(mean(post)), as.double(var(post)), as.double(cov(pre,post)), as.double(1.0), as.double(0.9), as.double(0.95), as.double(c(1.0,1.0)), as.integer(1), as.integer(10000), as.integer(1), as.double(0.025))





otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=TRUE, paired_data=TRUE, paired_analysis=FALSE, preN=20, premean=50, animalk=1, efficacyk=2, prek=2, postk=2)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=TRUE, paired_data=TRUE, paired_analysis=TRUE, preN=20, premean=50, animalk=1, efficacyk=2, prek=2, postk=2)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)


otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=FALSE, paired_data=TRUE, paired_analysis=FALSE, preN=20, premean=50, animalk=0.001, efficacyk=0.5, prek=1000, postk=1000)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=FALSE, paired_data=TRUE, paired_analysis=TRUE, preN=20, premean=50, animalk=0.001, efficacyk=0.5, prek=1000, postk=1000)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)



otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=FALSE, paired_data=FALSE, paired_analysis=FALSE, preN=50, premean=50, animalk=1, efficacyk=1, prek=1, postk=0.1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=TRUE, paired_data=FALSE, paired_analysis=FALSE, preN=50, premean=50, animalk=1, efficacyk=1, prek=1, postk=0.1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=1, efficacyk=1, prek=1, postk=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

# This works well when non-paired knows the true k
# So is it possible to guess the animal/efficacyk from the covariance and then add this to observed post k?


otu <- bayescount:::fecrt_power_wrap(reduction=0.9, use_truek=TRUE, k_change=0.001, paired_data=TRUE, paired_analysis=FALSE, preN=50, premean=50, animalk=0.5, efficacyk=0.15, prek=1, postk=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=0.15, prek=1, postk=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)


otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=FALSE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)


otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, rep_pre=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, rep_pre=2)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, rep_pre=0.5)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, edt_pre=1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, edt_pre=10)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired_data=TRUE, paired_analysis=TRUE, preN=50, premean=50, animalk=0.5, efficacyk=1, prek=1, postk=1, edt_pre=0.1)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)

stop('rep and edt pre dont work properly')



pre <- c(139,0,0,0,0,0,0,0,0,0)
post <- c(6,0,0,0,0,0,0,0,0,0)
res <- bayescount:::fecrt_analyses(pre, post, true_k=1, k_change=1)
res

139, 10, 1.000000, 13.900000, 411.211111, 6, 10, 1.000000, 0.600000, 0.488889, 0.000000, 1.000000, 0.900000, 0.950000, 1.000000, 1.000000, 1, 0, 10000, 0.016286, 0.601419

otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired=FALSE, iterations=1)
pre = c(7, 12, 8, 31, 37, 68, 74, 40, 1, 53);post = c(2, 0, 9, 1, 0, 0, 4, 1, 10, 2)
mean(pre)
mean(pre)^2 / (var(pre)-mean(pre))
mean(post)
mean(post)^2 / (var(post)-mean(post))

library('bayescount')
otu <- bayescount:::fecrt_power_wrap(reduction=0.9, paired=FALSE)
otp <- bayescount:::fecrt_power_wrap(reduction=0.9, paired=TRUE)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
sum(otp$p_1 < 0.025) / nrow(otp)
sum(otp$p_2 < 0.025) / nrow(otp)

getstats <- function(H0_1=0.9, H0_2=0.95, ...){
	otu <- bayescount:::fecrt_power_wrap(reduction=H0_1, H0_1=H0_1, H0_2=H0_2, paired=FALSE, ...)
	otp <- bayescount:::fecrt_power_wrap(reduction=H0_1, H0_1=H0_1, H0_2=H0_2, paired=TRUE, ...)
	power_1 <- c(sum(otu$p_1 < 0.025) / nrow(otu), sum(otp$p_1 < 0.025) / nrow(otp))
	err_2 <- c(sum(otu$p_2 < 0.025) / nrow(otu), sum(otp$p_2 < 0.025) / nrow(otp))
	
	otu <- bayescount:::fecrt_power_wrap(reduction=H0_2, H0_1=H0_1, H0_2=H0_2, paired=FALSE, ...)
	otp <- bayescount:::fecrt_power_wrap(reduction=H0_2, H0_1=H0_1, H0_2=H0_2, paired=TRUE, ...)
	err_1 <- c(sum(otu$p_1 < 0.025) / nrow(otu), sum(otp$p_1 < 0.025) / nrow(otp))
	power_2 <- c(sum(otu$p_2 < 0.025) / nrow(otu), sum(otp$p_2 < 0.025) / nrow(otp))

	ret <- expand.grid(Rate=NA, Method=c('Paired','Unpaired'), Hypothesis=1:2, Type=c('Power','Error'))[,4:1]
	ret$Rate <- c(power_1, power_2, err_1, err_2)
	return(ret)
}
getstats()


it <- 1000
pvals1=pvals2=obsred <- numeric(it)
premu <- 20
pb <- txtProgressBar(style=3)
for(i in 1:it){
	pre <- rnbinom(20, 1, mu=premu)
	post <- rnbinom(20, 1, mu=premu*0.11)
	#post <- rnbinom(20, 1, mu=pre*0.1)
	
	res <- bayescount:::fecrt_analyses(pre, post)
	pvals1[i] <- res$p1[res$Method=='BNB unpaired']
	pvals2[i] <- res$p2[res$Method=='BNB unpaired']
#	pvals1[i] <- res$p1[res$Method=='BNB paired']
#	pvals2[i] <- res$p2[res$Method=='BNB paired']
	
	obsred[i] <- 1- mean(post)/mean(pre)
	setTxtProgressBar(pb, i/it)
}
close(pb)
sum(pvals1 < 0.025) / it
sum(pvals2 < 0.025) / it
hist(obsred, breaks='fd')
hist(pvals1, breaks='fd')


library('bayescount')

preN <- 20
premean <- 20
prek <- 1
postk <- 1

otu <- bayescount:::fecrt_power_wrap(reduction=0.99, paired=FALSE, preN=preN, premean=premean, prek=prek, postk=postk)
otp <- bayescount:::fecrt_power_wrap(reduction=0.99, paired=TRUE, preN=preN, premean=premean, prek=prek, postk=postk)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
sum(otp$p_1 < 0.025) / nrow(otp)
sum(otp$p_2 < 0.025) / nrow(otp)


otu <- bayescount:::fecrt_power_wrap(reduction=0.95, paired=FALSE, preN=preN, premean=premean, prek=prek, postk=postk)
otp <- bayescount:::fecrt_power_wrap(reduction=0.95, paired=TRUE, preN=preN, premean=premean, prek=prek, postk=postk)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
sum(otp$p_1 < 0.025) / nrow(otp)
sum(otp$p_2 < 0.025) / nrow(otp)


otu <- bayescount:::fecrt_power_wrap(reduction=0.89, paired=FALSE, preN=preN, premean=premean, prek=prek, postk=postk)
otp <- bayescount:::fecrt_power_wrap(reduction=0.89, paired=TRUE, preN=preN, premean=premean, prek=prek, postk=postk)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
sum(otp$p_1 < 0.025) / nrow(otp)
sum(otp$p_2 < 0.025) / nrow(otp)

hist(otu$ObsEfficacy, breaks='fd', xlim=c(0.5,1))
hist(otu$p_1, breaks='fd')


otu <- bayescount:::fecrt_power_wrap(reduction=0.85, paired=FALSE, preN=preN, premean=premean, prek=prek, postk=postk)
otp <- bayescount:::fecrt_power_wrap(reduction=0.85, paired=TRUE, preN=preN, premean=premean, prek=prek, postk=postk)
sum(otu$p_1 < 0.025) / nrow(otu)
sum(otu$p_2 < 0.025) / nrow(otu)
sum(otp$p_1 < 0.025) / nrow(otp)
sum(otp$p_2 < 0.025) / nrow(otp)


it <- 10000
pvals1=pvals2=obsred <- numeric(it)
premu <- 20
pb <- txtProgressBar(style=3)
for(i in 1:it){
	pre <- rnbinom(20, 1, mu=premu)
	post <- rnbinom(20, 1, mu=premu*0.11)
	#post <- rnbinom(20, 1, mu=pre*0.1)
	
	obsred[i] <- 1- mean(post)/mean(pre)
	setTxtProgressBar(pb, i/it)
}
close(pb)

hist(obsred, breaks='fd', xlim=c(0.5,1))




if(FALSE){
	pow <- function(x,i) x^i
var1 <- var(pre)
var2 <- var(post)
cov12 <- cov(pre,post)
mu1 <- mean(pre)
mu2 <- mean(post)
varred = pow(mu2/mu1, 2) * (var1/pow(mu1,2) + var2/pow(mu2,2))
varred
varred = pow(mu2/mu1, 2) * (var1/pow(mu1,2) + var2/pow(mu2,2) - 2.0 * cov12 * sqrt(var1) * sqrt(var2) / (mu1 * mu2))
varred
Nd = length(pre)
r = mu2/mu1;
shape = (pow(r, 2) * Nd) / varred;
scale = varred / (r * Nd);
1-qgamma(c(0.975,0.025),shape,scale=scale)
}
bayescount:::fecrt_analyses(pre, post)
cor(pre,post)
bayescount:::fecrt_analyses(pre, post, replicates_ratio=1.5, edt_ratio=0.7)
bayescount:::fecrt_analyses(pre, post, edt_ratio=0.5)
bayescount:::fecrt_analyses(pre, c(0,0))

bayescount:::waavp_ci(pre, post, tail=0.025)
bayescount:::dobson_ci(pre, post)
bayescount:::conjbeta_ci(pre, post)
bayescount:::asymptotic_ci(pre, post)
bayescount:::fecrt_bnb(pre, post)

pre <- rnbinom(10, 1, mu=10)
post <- rnbinom(10, 1, mu=0.01)

bayescount:::waavp_ci(pre, post)
bayescount:::dobson_ci(pre, post)
bayescount:::conjbeta_ci(pre, post)
bayescount:::asymptotic_ci(pre, post)
bayescount:::fecrt_bnb(pre, post)

bayescount:::fecrt_power_comparison_wrap(reduction=0.99)
bayescount:::fecrt_power_comparison_wrap(reduction=0.95)
bayescount:::fecrt_power_comparison_wrap(reduction=0.9)
bayescount:::fecrt_power_comparison_wrap(reduction=0.8)

bayescount:::fecrt_power_wrap(reduction=0.99, H0_1=c(0.8,0.9), H0_2=c(0.99, 0.95))
bayescount:::fecrt_power_wrap(reduction=0.95, H0_1=c(0.8,0.9), H0_2=c(0.99, 0.95))
bayescount:::fecrt_power_wrap(reduction=0.9, H0_1=c(0.8,0.9), H0_2=c(0.99, 0.95))
bayescount:::fecrt_power_wrap(reduction=0.8, H0_1=c(0.8,0.9), H0_2=c(0.99, 0.95))

# bnb_pres(preN, presum, preK, postN, postsum, postK, H0_1, H0_2, prob_priors, delta, beta_iters){
bayescount:::bnb_ps(preN=100, presum=556724, preK=0.011865, postN=100, postsum=5405, postK=0.028920, H0_1=0.9, H0_2=0.95, prob_priors=c(1,1), delta=1, beta_iters=10^4)
mu = 0.9999
s2 = 0.000100
if(s2 > (mu * (1-mu)))
	s2 <- (mu * (1-mu))*(1-10^-6)
ab = mu * (1.0 - mu) / s2 -1
mu * ab
(1-mu)*ab




# The asymptotic p-values seem to be quite close although asymptotic ones are always smaller:
predata <- rnbinom(1000, 0.1, mu=100)
postdata <- rnbinom(1000, 0.1, mu=6)
H0_1 <- 0.9
H0_2 <- 0.95
ps <- bayescount:::fecrt_bnb(predata, postdata, H0_1=H0_1, H0_2=H0_2)
format(ps[c(1,3)], scientific=FALSE)
format(ps[c(2,4)], scientific=FALSE)
# Something is broken with BNB when preN!=postN

red_se <- sqrt(((mean(postdata)^2) / (mean(predata)^4) * var(predata) + (mean(predata)^-2) * var(postdata)) / min(length(predata), length(postdata)))
meanchange = (1.0 - H0_1);
pnorm(meanchange,mean(postdata)/mean(predata),red_se,FALSE);
meanchange = (1.0 - H0_2);
pnorm(meanchange,mean(postdata)/mean(predata),red_se,TRUE);



# Small k is fine whatever:
bayescount:::pbnb1(10^4, 0.1, 16865.085493, 19481.403311)
bayescount:::pbnb1(10^6, 0.1, 16865.085493, 19481.403311)

# Small post sum is fine whatever:
bayescount:::pbnb1(100, 1, 16865.085493, 19481.403311)
bayescount:::pbnb1(100, 100, 16865.085493, 19481.403311)
bayescount:::pbnb1(100, 1000, 16865.085493, 19481.403311)

# Big k and big post sum is a problem:
bayescount:::pbnb1(1000, 1000, 16865.085493, 19481.403311)
# So use the asymptotic unpaired method when effK is big and obs sum is big
# Can still get p-value from asymptotic method

k <- 1000
alpha <- 168
beta <- 194
sum <- 100
bayescount:::pbnb1(sum, k, alpha, beta)
mean <- k*beta / (alpha-1)
var <- k*(alpha+k-1)*beta*(alpha+beta-1) / ((alpha-2)*(alpha-1)^2)
pnorm(sum, mean, sqrt(var))

predata <- rnbinom(10000, 10, mu=1000)
postdata <- rnbinom(10000, 10, mu=mean(predata)*0.05)
postdata <- rnbinom(10000, 10, mu=0.05)
sum(postdata)
bayescount:::fecrt_bnb(predata, postdata, H0_1=0.9, H0_2=0.95)

testing <- bayescount:::fecrt_power_wrap(reduction=0.95, preN=10^4, iterations=10, premean=1000)
testing

bayescount:::pbnb2(1036, 1000, 16865.085493, 19481.403311)
bayescount:::pbnb2(1036, 1000, 16865.085493, 19481.403311)

r=471
a=-7475.245189
n=-753.642449
bn=9352.048569
rp=r+1
((r-a)*(r-n))/(rp*(bn+rp))


A <- 10
dN <- 34
dn <- 12
dr <- 0.01
A * log((A*dN)/(dn*dr))
A * ((log(A)+log(dN)) - (log(dn)+log(dr)))

A * (log(A) + log(dn) - log(dN) - log(dr))

library('tidyverse')

ans <- bayescount:::fecrt_power_wrap(reduction=0.95, H0_1=c(0.8,0.9), H0_2=c(0.99, 0.95), iterations=10^4)
ans %>% group_by(H0_1, H0_2, TrueEfficacy, Classification) %>% tally

ans <- bayescount:::fecrt_power_wrap(reduction=0.95, H0_1=c(0.8), H0_2=c(0.99), iterations=10^4)
ans %>% group_by(H0_1, H0_2, Classification) %>% tally

ans <- bayescount:::fecrt_power_wrap(reduction=0.95, H0_1=c(0.9), H0_2=c(0.95), iterations=10^4)
ans %>% group_by(H0_1, H0_2, Classification) %>% tally

length(rv$p_1)

rv=ans

### OLDER


tt <- bayescount:::asymptotic_ci_wrap(1:10, 1:10)


fpr <- bayescount:::fecrt_power_wrap(paired=FALSE, iterations=10000, trues=0.90, rep_pre=1, rep_post=10, edt_pre=1, edt_post=0.1, preN=40, postN=40)

fpr <- bayescount:::fecrt_power_wrap(iterations=100000, red=0.9, pair_type=0)
fpr$classifications[1] / 100000
fpr$classifications[4] / 100000
fpr <- bayescount:::fecrt_power_wrap(iterations=10000, red=0.9, pair_type=1, delta=0)
fpr$classifications[1] / 10000
fpr$classifications[4] / 10000


# TOFIX:
	# Failure values:  p1: nan, p2: nan  (presum: 235, preN: 20, prek: 0.183922, postsum: 16, postN: 20, postk: 6.080000, lt: 0.900000, ut: 0.950000)
rv <- bayescount:::fecrt_pee_direct(235, 20, 0.183922, 16, 20, 6.080000, 0.9, 0.95, delta=0)
rv$p_1; rv$p_2
rv <- bayescount:::fecrt_pee_direct(235, 20, 0.183922, 16, 20, 6.080000, 0.9, 0.95, delta=1)
rv$p_1; rv$p_2
rv <- bayescount:::fecrt_pee_direct(235, 20, 0.183922, 16, 20, 6.080000, 0.9, 0.95, delta=2)
rv$p_1; rv$p_2

	# Failure values:  p1: 0.626118, p2: nan  (presum: 148, preN: 20, prek: 0.147283, postsum: 14, postN: 20, postk: 0.721705, lt: 0.900000, ut: 0.950000)
	rv <- bayescount:::fecrt_pee_direct(148, 20, 0.147283, 14, 20, 0.721705, 0.9, 0.95)
	rv$p_1; rv$p_2
	
	rv <- bayescount:::fecrt_pee_direct(148, 20, 0.147283, 100, 20, 0.721705, 0.9, 0.95)
	rv$p_1; rv$p_2
	rv <- bayescount:::fecrt_pee_direct(148, 20, 0.147283, 0, 20, 0.721705, 0.9, 0.95)
	rv$p_1; rv$p_2

	pre <- rnbinom(20,1,mu=20)
	post <- rnbinom(20,1,mu=1.5)
	
	rv <- bayescount:::fecrt_pee_wrap(pre, post, delta=0)	
	rv$p_1; rv$p_2
	rv <- bayescount:::fecrt_pee_wrap(pre, post, delta=2)	
	rv$p_1; rv$p_2

	post <- rnbinom(20,1,mu=0.1)
	
	rv <- bayescount:::fecrt_pee_wrap(pre, post, delta=0)	
	rv$p_1; rv$p_2
	rv <- bayescount:::fecrt_pee_wrap(pre, post, delta=2)	
	rv$p_1; rv$p_2
	
		
fpr$classifications[6] / 10000
fpr$classifications[9] / 10000
summary(fpr$obsred)

##### Christian's approach

r_mle <- function(pre, post){
  stopifnot(length(post)==length(pre))
  return(sum(post) / sum(pre))
}
mu1_mle <- function(pre){
  return(sum(pre) / length(pre))
}

sig2 <- function(pre, post){
  stopifnot(length(pre)==length(post))
  N <- length(pre)
  r <- r_mle(pre, post)
  mu1 <- mu1_mle(pre)
  p1 <- 1/(r*mu1)
  p2 <- 1/(r*mu1) * 1/(N-1) * sum((pre-mean(pre))^2)
  p3 <- - (mu1 / (r*mu1)^2)
  p4 <- - (mu1 / (r*mu1)^2) * 1/(N-1) * sum((post-mean(post))^2)
  return(p1*p2 + p3*p4)
}

# Wrong but kind of works:
sig2 <- function(pre, post){
  stopifnot(length(pre)==length(post))
  N <- length(pre)
  r <- r_mle(pre, post)
  mu1 <- mu1_mle(pre)
  p1 <- 1/(r*mu1)
  p2 <- 1/(N-1) * sum((pre-mean(pre))^2) * 1/(r*mu1)   # Can also be / 1/(r*mu1) so OK when 1/(N-1)(r*mu1) == (r*mu1)/(N-1) ??
  p3 <- - (mu1 / (r*mu1)^2)
  p4 <- - (mu1 / (r*mu1)^2) * 1/(N-1) * sum((post-mean(post))^2)
  return(1/(N-1)*p1*p2 + 1/(N-1)*p3*p4)
}  # with N=1000 and tr=0.2 and any k/mu


# Work off this??
sig2 <- function(pre, post, k){
  stopifnot(length(pre)==length(post))
  N <- length(pre)
  r <- r_mle(pre, post)
  mu1 <- mu1_mle(pre)
  p1 <- 1/(r*mu1)
  p2 <- 1/(N-1) * sum((pre-mean(pre))^2) * 1/(r*mu1)
  p3 <- - (mu1 / (r*mu1)^2)
  p4 <- - (mu1 / (r*mu1)^2) * 1/(N-1) * sum((post-mean(post))^2)
  
  # mu1+k*mu1^2; (1/(N-1))*sum((pre-mean(pre))^2)
  # r*mu1+k*r^2*mu1^2; (1/(N-1))*sum((post-mean(post))^2)
  
  p1 <- 1/(N-1) * p1
  p3 <- 1/(N-1) * p3
  return(p1*p2 + p3*p4)
}


sig2 <- function(pre, post, k){
  stopifnot(length(pre)==length(post))
  N <- length(pre)
  r <- r_mle(pre, post)
  mu1 <- mu1_mle(pre)
  p1 <- 1/(r*mu1)
  p2 <- 1/(N-1) * sum((pre-mean(pre))^2) * 1/(r*mu1)
  p3 <- - (mu1 / (r*mu1)^2)
  p4 <- - (mu1 / (r*mu1)^2) * 1/(N-1) * sum((post-mean(post))^2)
  
  p1 <- 1/(N-1) * p1
  p3 <- 1/(N-1) * p3
  return(p1*p2 + p3*p4)
}


sig2 <- function(pre, post, k){
  stopifnot(length(pre)==length(post))
  N <- length(pre)
  r <- r_mle(pre, post)
  mu1 <- mu1_mle(pre)
  p1 <- 1/(r*mu1)
  p2 <- sum((pre-mean(pre))^2) * 1/(N-1) * 1/(r*mu1)
  p3 <- - (mu1 / (r*mu1)^2)
  p4 <- - (mu1 / (r*mu1)^2) * 1/(N-1) * sum((post-mean(post))^2)
  
  # mu1 + (mu1^2)/k; (1/(N-1))*sum((pre-mean(pre))^2)
  # r*mu1+r^2*mu1^2/k; (1/(N-1))*sum((post-mean(post))^2)
  
#  p1 <- 1/(N-1) * p1
#  p3 <- 1/(N-1) * p3
  return(p1*p2 + p3*p4)
}


mu <- 1000
k <- 10
r <- 0.5
(mu + 1/k * mu^2) / (r*mu)^2 + mu^2/(r*mu)^4 * (r*mu + 1/k*r^2*mu^2)

s1 <- rnbinom(10000, k, mu=mu)
s2 <- rnbinom(10000, k, mu=mu*r)
ratio <- s2/s1
mean(ratio[ratio < Inf])
median(ratio[ratio < Inf])
sd(ratio[ratio < Inf])


r=a1=a2 <- 1
mu1 <- sample(250:1000,1)
delta2a = (mu1 + a1 * mu1^2) / (r*mu1)^2 + mu1^2/(r*mu1)^4 * (r*mu1 + a2*r^2*mu1^2)
delta2a
delta2b

# Fix r at 1:
r <- 0.8
# Randomly sample other parameters:
mu1 <- exp(runif(1, log(1), log(100)))
a1 <- runif(1, 0.1, 5)
a2 <- runif(1, 0.1, 5)
N <- sample(250:1000, 1)

# Monte Carlo approximation to the distribution of the ML estimator for r:
ratios <- sapply(1:1000, function(x) sum(rnbinom(N,size=1/a2,mu=mu1*r))/sum(rnbinom(N,size=1/a1,mu=mu1)))
# Generally good approximation to normal for the parameter ranges chosen:
hist(ratios, breaks='fd')

# Monte Carlo approximation to the variance of the distribution of ML estimator:
var(ratios)
# Direct calculation of the same quantity:
delta2 = (mu1 + a1 * mu1^2) / (r*mu1)^2 + mu1^2/(r*mu1)^4 * (r*mu1 + a2*r^2*mu1^2)
1/N * delta2
# Relative error:
(1/N * delta2 - var(ratios)) / var(ratios)
# Good agreement
# ECDF plots of both curves for good measure:
plot.ecdf(ratios, col='blue')
plot.ecdf(rnorm(10000, r, sd=sqrt(1/N * delta2)), col='red', add=TRUE)


### Best attempt:

# Obs vs calculated variances:
sum((pre-mean(pre))^2) / (N-1)
var(pre)
mu + mu^2 / k

sum((post-mean(post))^2) / (N-1)
var(post)
tr*mu + tr^2*mu^2 / k

# Assuming known mu, r, k:
sig2k <- function(mu1, r, k){
  a1 <- 1/k
  a2 <- 1/k
  return((mu1+a1+mu1^2)/(r*mu1)^2 + (mu1^2/(r*mu1)^4 * (r*mu1 + a2*r^2*mu1^2)))
}


N <- 1000 #sample(250:1000,1)
k <- 10  #round(runif(1,0.1,5), 1)
mu <- 10000 #round(runif(1,5,100))
tr <- 0.2 #round(runif(1,0,0.25), 2)

iters <- 1000
reds <- numeric(iters)
for(i in 1:iters){
  post <- rnbinom(N, size=k, mu=mum*rm)
  pre <- rnbinom(N, k, mu=mum)
  reds[i] <- sum(post) / sum(pre)
}
hist(reds, freq=FALSE, main=paste('k:', k, ', mu:', mu, ', r:', tr, ', N:', N), breaks='fd')
curve(dnorm(x, tr, sd=sqrt(1/N * sig2k(mu, tr, k))), col='red', add=TRUE)
sqrt(1/N * sig2k(mu, tr, k))

# If r=k=1:
a1=a2=r <- 1
m <- mu1
(mu1+a1+mu1^2)/(r*mu1)^2 + (mu1^2/(r*mu1)^4 * (r*mu1 + a2*r^2*mu1^2))
(m+1+m^2) / m^2 + (m+m^2)/m^2
m <- 10:100
var <- (m+1+m^2) / m^2 + (m+m^2)/m^2
plot(m, sqrt(var/50), type='l')


sig2d(mu,tr,k)

trues2 <- (mu1+k*mu1^2)/(r*mu1)^2 + (mu1^2/(r*mu1)^4 * (r*mu1 + k*r^2*mu1^2))

curve(dnorm(x,mu-(tr*mu),sqrt(trues2/N)), from=0, to=mu)


# There is some interaction between tr and N, but k and mu about right
# But this is just coincidence for correct sd

N <- 1000
pre <- rnbinom(N, size=k, mu=mu)
post <- rnbinom(N, size=k, mu=mu*tr)
s2 <- sig2(pre, post, 1/k)



trues2
s2

curve(dnorm(x, r_mle(pre, post), sd=sqrt(1/N * s2)), main='Sampling distribution of r')

###### Compare this to fecrt_precision

rm <- r_mle(pre, post)
mum <- mu1_mle(pre)
mdiff <- mean(pre) - mean(post)

iters <- 1000
diffs=reds <- numeric(iters)
for(i in 1:iters){
  post <- rnbinom(N, k, mu=mum*rm)
  pre <- rnbinom(N, k, mu=mum)
  reds[i] <- sum(post) / sum(pre)
  diffs[i] <- mean(pre) - mean(post)
}
hist(reds, freq=FALSE, main=paste('k:', k, ', mu:', mu, ', r:', tr, ', N:', N))
curve(dnorm(x, rm, sd=sqrt(1/N * s2)), col='red', add=TRUE)

hist(diffs, freq=FALSE, main=paste('k:', k, ', mu:', mu, ', r:', tr, ', N:', N))
curve(dnorm(x, mdiff, sd=sqrt(sig2d(mu,tr,k)/N)), col='red', add=TRUE)

ar <- rnorm(iters, rm, sd=sqrt(1/N * s2))
qqplot(reds, ar)
abline(0,1)



sums <- sapply(1:1000, function(x) return(sum(rnbinom(N,k,prob))))
inds <- rnbinom(1000,N*k,prob)

qqplot(sums, inds); abline(0,1)


N <- 20
k <- 1.5
pre <- rnbinom(N, k, mu=20)
post <- rnbinom(N, k, mu=20*0.1)

fecrt_pval <- bayescount:::reduction_pval

bayescount:::reduction_pval(pre, post, delta=FALSE)
sum(post)

post <- rnbinom(N, k, mu=20*0.2)
pre <- rnbinom(N, k, mu=20)
fecrt_pval(pre, post, true_k=1)
fecrt_pval(pre, post, true_k=NA)
post <- rnbinom(N, k, mu=20*0.075)
fecrt_pval(pre, post)
post <- rnbinom(N, k, mu=20*0.01)
fecrt_pval(pre, post)

post <- rnbinom(N, k, mu=20*0.1)
fecrt_pval(pre, post)

bayescount:::pbeta_nbinom(12, 1.5, 27, 2800)
meanp <- 27/(27+2800)
meanp
sim <- rnbinom(1000, 1.5, meanp)

bayescount:::pbeta_nbinom(12, 1.5, 26.617355, 2794.419795)
bayescount:::pghyper(12, -2794.419795, -27.967380, 25.617355)

pre

it <- 10^7

s <- Sys.time()
cbt <- .C('conjbeta_ci_wrap', preN=as.integer(N), presum=as.integer(sum(pre)), preK=as.numeric(k), postN=as.integer(N), postsum=as.integer(sum(post)), postK=as.numeric(k), iters=as.integer(it), prob_priors=as.numeric(c(1,1)), lci_c=as.numeric(0.025), uci_c=as.numeric(0.975), ci_l=as.numeric(0), ci_u=as.numeric(0), PACKAGE='bayescount')
cbt
difftime(Sys.time(), s)

s <- Sys.time()
b1 <- rbeta(it, N*k +1, sum(pre) +1)
b2 <- rbeta(it, N*k +1, sum(post) +1)
m1 <- k / b1 - k;
m2 <- k / b2 - k;
quantile(1- m2/m1, prob=c(0.025,0.975))
difftime(Sys.time(), s)

# Why doesn't reserve and push_back work on the laptop?




ps <- c(0.1,0.9)
qbeta(ps, 10, 5) / qbeta(ps[2:1], 10, 7)



get_waavp <- function(pre, post){
  arith.pre <- mean(pre)
  arith.post <- mean(post)
  var.pre <- var(pre)
  var.post <- var(post)
  var.red <- ifelse(var.post==0, 0, var.post/(length(post) * arith.post^2)) + var.pre/(length(pre) * arith.pre^2)
  
  return((1 - (arith.post/arith.pre * exp(c(-1,1)*2.048 * sqrt(var.red)))))
}

N <- 10
k <- 0.5
pre <- rnbinom(N, k, mu=20)
post <- rnbinom(N, k, mu=20*0.1)

unlist(fecrt_pval(pre, post))
fecrt_pval(pre, post, delta=FALSE)


round(get_pee(pre, post), 5)

tval <- 2.048
wvp <- .C('waavp_ci_wrap', presum=as.integer(sum(pre)), predata=as.integer(pre), preN=as.integer(length(pre)), postsum=as.integer(sum(post)), postdata=as.integer(post), postN=as.integer(length(post)), tval=as.numeric(tval), ci_l=as.numeric(0), ci_u=as.numeric(0))
wvp
get_waavp(pre, post)

res <- testpees(1000)
summary(res[,1]-res[,5])
summary(res[,2]-res[,6])
# Nothing wrong with underlying pval function

res2 <- fecrt_pvals(iterations=10000, N=10, pre_mean=50, reduction=0.01, pre_k=1)

table(Sus=res[,1]<0.025, Res=res[,2]<0.025)
table(Sus=res[,5]<0.025, Res=res[,6]<0.025)
table(Sus=res2[,1]<0.025, Res=res2[,2]<0.025)


# Can provoke a high error rate only with big discrepancy between pre and post k
# This can be offset with k_change=0.8 though
prek=1.5
postk=1.5
res2 <- fecrt_pvals(N=10, pre_mean=100, reduction=0.11, pre_k=prek, post_k=postk, k_change=1)
table(Sus=res2[,1]<0.025, Res=res2[,2]<0.025)

res <- testpees(N=20, premean=50, reduction=0.04, k=prek, postk=postk)
table(Sus=res[,1]<0.025, Res=res[,2]<0.025)
table(Sus=res[,5]<0.025, Res=res[,6]<0.025)


lci_c <- 0.025
uci_c <- 0.975
dobson_priors <- c(1,1)
dobson <- .C('dobson_ci_wrap', presum=as.integer(sum(pre)), postsum=as.integer(sum(post)), lci_c=as.numeric(lci_c), uci_c=as.numeric(uci_c), dobson_priors=as.numeric(dobson_priors), ci_l=as.numeric(0), ci_u=as.numeric(0))
wvp
dobson

# Matches publication:
dobson <- .C('dobson_ci_wrap', pre=as.integer(100), preN=as.integer(1), post=as.integer(10), postN=as.integer(1), lci_c=as.numeric(lci_c), uci_c=as.numeric(uci_c), dobson_priors=as.numeric(dobson_priors), ci_l=as.numeric(0), ci_u=as.numeric(0))


var(rgamma(1000, 1, scale=10))

res <- fecrt_pvals(10000, N=10, pre_mean=10, reduction=0.01, pre_k=0.2, H0_1=0.5, H0_2=0.5, delta=FAL)
table(Sus = res$p1 < 0.025, Res = res$p2 < 0.025)

testpees(1000)

hist(res$p2)

res <- fecrt_power(N=10, pre_mean=10, reduction=0.04, pre_k=0.4, psig=0.025, beta_iters=1000)
apply(res, c(3,4), sum)
apply(res, c(1,4), sum)

apply(res, c(1,4,5), sum)
apply(res, c(2,4,5), sum)


res <- fecrt_power(10000, N=10, pre_mean=20, reduction=0.11, pre_k=0.8, post_k=0.8, psig=0.025)
apply(res, c(1,3,4), sum)
apply(res, c(2,3,4), sum)

# Dobson method has variable false positive/negative rates

387/10000

res1 <- .C('fecrt_power_comparison', iters=as.integer(iters), preN=as.integer(length(pre)), postN=as.integer(length(post)), maxN=as.integer(length(maxN)), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), edt_change=as.numeric(edt_change), prob_priors=as.numeric(prob_priors), kchange=as.numeric(kchange), truek=as.numeric(truek), usetruek=as.integer(usetruek), delta=as.integer(TRUE), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0), package='bayescount')

pre=as.integer(pre), 

fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, double *premean, double *reduction, double *edt_change, double *animalk, double *prek, double *postk, double *H0_1, double *H0_2, double *lci_c, double *uci_c, double *dobson_priors, double *tval, double *psig, double *prob_priors, double *kchange, double *truek, int *usetruek, int *delta, int *beta_iters, int *predata, int *postdata, int *classifications)


N <- 10
k <- 1
pre <- rnbinom(N, k, mu=20)
post <- rnbinom(N, k, mu=20*0.1)
H0_1 <- 0.95
H0_2 <- 0.95
edt_change <- 1
prob_priors <- c(1,1)
kchange <- 1
truek <- 1
usetruek <- FALSE
delta <- TRUE
beta_iters <- 10^6
# Max 10^8 - 10^6 is pretty instant and gives the same results

res1 <- .C('fecrt_pee_wrap', pre=as.integer(pre), preN=as.integer(length(pre)), post=as.integer(post), postN=as.integer(length(post)), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), edt_change=as.numeric(edt_change), prob_priors=as.numeric(prob_priors), kchange=as.numeric(kchange), truek=as.numeric(truek), usetruek=as.integer(usetruek), delta=as.integer(TRUE), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0), package='bayescount')
res2 <- .C('fecrt_pee_wrap', pre=as.integer(pre), preN=as.integer(length(pre)), post=as.integer(post), postN=as.integer(length(post)), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), edt_change=as.numeric(edt_change), prob_priors=as.numeric(prob_priors), kchange=as.numeric(kchange), truek=as.numeric(truek), usetruek=as.integer(usetruek), delta=as.integer(FALSE), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0), package='bayescount')
res1$p_1
res2$p_1
res1$p_2
res2$p_2


iters <- 1000
preN=postN <- 10
premean <- 100
reduction <- 0.1
edt_change <- 1
animalk <- 1  # Can be >100 -> means unpaired
prek <- 1
postk <- 1

H0_1 <- c(0.95, 0.9)
H0_2 <- c(0.9, 0.95)

stopifnot(length(H0_1)==length(H0_2))
H0_N <- length(H0_1)
maxN <- max(preN, postN)

p_1=p_2 <- rep(0, iters*H0_N)

powres <- .C('fecrt_power', iters=as.integer(iters), preN=as.integer(length(pre)), postN=as.integer(length(post)), maxN=as.integer(maxN), premean=as.numeric(premean), reduction=as.numeric(reduction), edt_change=as.numeric(edt_change), animalk=as.numeric(animalk), prek=as.numeric(prek), postk=as.numeric(postk), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), H0_N=as.integer(H0_N), prob_priors=as.numeric(prob_priors), kchange=as.numeric(kchange), truek=as.numeric(truek), usetruek=as.integer(usetruek), delta=as.integer(delta), beta_iters=as.integer(beta_iters), predata=as.integer(rep(0,maxN)), postdata=as.integer(rep(0,maxN)), p_1=as.numeric(p_1), p_2=as.numeric(p_2), package='bayescount')

p_1 <- powres$p_1
p_2 <- powres$p_2
dim(p_1) <- c(iters, H0_N)
dim(p_2) <- c(iters, H0_N)

head(p_1)
head(p_2)


iters <- 1000
res <- expand.grid(iter=1:iters, r=NA, pal=NA, pau=NA, ual=NA, uau=NA, wal=NA, wau=NA, pv1=NA, pv2=NA)

pb <- txtProgressBar(style=3)
for(i in 1:iters){
  pre <- rnbinom(N, 1, mu=100)
  post <- rnbinom(N, 1, mu=5)
  
  as <- bayescount:::asymptotic_ci_wrap(pre, post)
  wa <- bayescount:::waavp_ci_wrap(pre, post)
  pv <- bayescount:::fecrt_pee_wrap(pre, post)
  
  res$r[i] <- 1 - mean(post) / mean(pre)
  res$pal[i] <- as$p_ci_l
  res$pau[i] <- as$p_ci_u
  res$ual[i] <- as$u_ci_l
  res$uau[i] <- as$u_ci_u
  res$wal[i] <- wa$ci_l
  res$wau[i] <- wa$ci_u
  res$pv1[i] <- pv$p_1
  res$pv2[i] <- pv$p_2
  
  setTxtProgressBar(pb, i/iters)
}
close(pb)



get_type_ci <- function(lci, uci, m, t){
  if(uci < m){
    tp <- 1
  }else if(uci < t && lci < m){
    tp <- 2
  }else if(lci < m && uci >= t){
    tp <- 3
  }else if(lci >= m && uci < t){
    tp <- 4
  }else if(lci >= m && uci >= t){
    tp <- 5
  }else if(lci >= t){
    tp <- 6
  }else{
    stop('Error assigning type')
  }
  return(tp)
}
get_type_pv <- function(p1, p2, s=0.025){
  if(p2 <= s && p1 > s){
    tp <- 1.5
  }else if(p2 > s && p1 > s){
    tp <- 3
  }else if(p2 <= s && p1 <= s){
    tp <- 4
  }else if(p2 > s && p1 <= s){
    tp <- 5.5
  }else{
    stop('Error assigning type')
  }
  return(tp)
}

iters <- 1000
res <- expand.grid(iter=1:iters, N=50, k=1, mu=10, true=0.85, Target=0.95, tol=0.05, o1=NA, r=NA, ap=NA, au=NA, wa=NA, pv=NA)

pb <- txtProgressBar(style=3)
for(i in 1:iters){
  pre <- rnbinom(res$N[i], 1/res$k[i], mu=res$mu[i])
  post <- rnbinom(res$N[i], 1/res$k[i], mu=res$mu[i]*(1-res$true[i]))
  
  as <- bayescount:::asymptotic_ci_wrap(pre, post)
  wa <- bayescount:::waavp_ci_wrap(pre, post)
  pv <- bayescount:::fecrt_pee_wrap(pre, post, res$Target[i] - res$tol[i], res$Target[i])
  
  res$r[i] <- 1 - mean(post) / mean(pre)
  
  res$ap[i] <- get_type_ci(as$p_ci_l, as$p_ci_u, res$Target[i] - res$tol[i], res$Target[i])
  res$au[i] <- get_type_ci(as$u_ci_l, as$u_ci_u, res$Target[i] - res$tol[i], res$Target[i])
  res$wa[i] <- get_type_ci(wa$ci_l, wa$ci_u, res$Target[i] - res$tol[i], res$Target[i])
  res$pv[i] <- get_type_pv(pv$p_1, pv$p_2, 0.025)
  
  res$o1[i] <- sum(post)==0
  
  setTxtProgressBar(pb, i/iters)
}
close(pb)

library('tidyverse')
lres <- res %>% gather(Method, Type, -iter, -N, -k, -mu, -true, -Target, -tol, -o1, -r)
#lres$Type <- factor(lres$Type, levels=c(1,1.5,2,3,4,5,5.5,6))
#levels(lres$Type) <- c(1.5,1.5,1.5,3,4,5.5,5.5,5.5)
#levels(lres$Type) <- c(1.5,1.5,1.5,3,4,5.5,5.5,5.5)

lres %>% group_by(Method, Type, N, k, mu, true, Target, tol, o1) %>% summarise(num=n()) %>% group_by(Method, N, k, mu, true, Target, tol, o1) %>% mutate(rate=num/sum(num))
lres %>% group_by(Method, Type, N, k, mu, true, Target, tol) %>% summarise(num=n()) %>% group_by(Method, N, k, mu, true, Target, tol) %>% mutate(rate=num/sum(num))

