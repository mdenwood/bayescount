library('bayescount')


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

