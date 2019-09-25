# Ensure results are replicable:
set.seed(2017-08-29)

library('bayescount')

if(FALSE){
	

combinek <- function(a, b) ((a*b)/(a+b+1))
combinekk <- function(a, b, c) combinek(a, combinek(b, c))

get_type_pv <- function(p1, p2, s=0.025){
  if(p2 <= s && p1 > s){
    tp <- "Reduced Efficacy"
	nr <- 1
  }else if(p2 > s && p1 > s){
    tp <- "Inconclusive"
	nr <- 2
  }else if(p2 <= s && p1 <= s){
    tp <- "Marginal Efficacy"
	nr <- 3
  }else if(p2 > s && p1 <= s){
    tp <- "Adequate Efficacy"
	nr <- 4
  }else{
    stop('Error assigning type')
  }
  return(nr)
}


# Test a simple unpaired model:
cit <- 100000
rit <- 10000

anik <- 2
effk <- 2
prek <- 3
postk <- 3

preN <- 20
postN <- 20
premean <- 10

rep_pre <- 1
rep_post <- 1
edt_pre <- 0.5
edt_post <- 0.25

reduction <- 0.88


#### Simple unpaired test

kpre <- combinek(anik, prek)
kpost <- combinekk(anik, effk, postk)

c1 <- bayescount:::fecrt_power_wrap(reduction, pair_type=0, preN=preN, postN=postN, animalk=anik, efficacyk=effk, prek=prek, postk=postk, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)
c2 <- bayescount:::fecrt_power_wrap(reduction, pair_type=1, preN=preN, postN=postN, animalk=anik, efficacyk=effk, prek=prek, postk=postk, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)
c3 <- bayescount:::fecrt_power_wrap(reduction, pair_type=0, preN=preN, postN=postN, animalk=Inf, efficacyk=Inf, prek=kpre, postk=kpost, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)

100*round(c1$classifications[1:4]/cit, 3)
100*round(c3$classifications[1:4]/cit, 3)
stopifnot(chisq.test(rbind(c1$classifications[c(1,2,4)], c3$classifications[c(1,2,4)]))$p.value > 0.05)
# Expected to be slightly different:
100*round(c2$classifications[1:4]/cit, 3)

c_reds <- c3$obsred
c_class <- c3$classifications[1:4]
stopifnot(abs(mean(c_reds-reduction)) < 0.01)

r_reds <- numeric(rit)
r_class <- numeric(4)

for(i in 1:rit){
	pre <- rnbinom(preN, kpre, mu=premean/edt_pre)
	post <- rnbinom(postN, kpost, mu=(premean*(1-reduction))/edt_post)
	
	r_reds[i] <- 1 - ((mean(post) * edt_post) / (mean(pre) * edt_pre))
	
	suppressWarnings(pv <- bayescount:::fecrt_pee_wrap(pre, post, edt_pre=edt_pre, edt_post=edt_post, rep_pre=rep_pre, rep_post=rep_post))
	r_class[get_type_pv(pv$p_1, pv$p_2)] <- r_class[get_type_pv(pv$p_1, pv$p_2)]+1	
}
summary(r_reds)
summary(c_reds)
reduction
# qqplot(r_reds, c_reds); abline(0,1)
stopifnot(t.test(r_reds, c_reds)$p.value > 0.05)

r_class/rit
c_class/cit
stopifnot(chisq.test(rbind(r_class[c(1,2,4)], c_class[c(1,2,4)]))$p.value > 0.05)



##### Unpaired with replicates

rep_pre <- 2
rep_post <- 5
edt_pre <- 0.25
edt_post <- 0.5

te <- try(c1 <- bayescount:::fecrt_power_wrap(reduction, pair_type=0, preN=preN, postN=postN, animalk=anik, efficacyk=effk, prek=prek, postk=postk, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit), silent=TRUE)
stopifnot(inherits(te, 'try-error'))
c2 <- bayescount:::fecrt_power_wrap(reduction, pair_type=1, preN=preN, postN=postN, animalk=anik, efficacyk=effk, prek=prek, postk=postk, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)

100*round(c2$classifications[1:4]/cit, 3)

c_reds <- c2$obsred
c_class <- c2$classifications[1:4]
stopifnot(abs(mean(c_reds-reduction)) < 0.01)

r_reds <- numeric(rit)
r_class <- numeric(4)

txk <- combinek(anik, effk)
premu <- premean/edt_pre
postmu <- (premean * (1-reduction)) / edt_post
for(i in 1:rit){
	ctmeans <- rgamma(preN, anik, rate=anik/premu)
	txmeans <- rgamma(postN, txk, rate=txk/postmu)
	
	pre <- sapply(1:preN, function(x) return(sum(rnbinom(rep_pre, prek, mu=ctmeans[x]))))
	post <- sapply(1:postN, function(x) return(sum(rnbinom(rep_post, postk, mu=txmeans[x]))))
	
	r_reds[i] <- 1 - ((mean(post)/rep_post * edt_post) / (mean(pre)/rep_pre * edt_pre))
	
	suppressWarnings(pv <- bayescount:::fecrt_pee_wrap(pre, post, edt_pre=edt_pre, edt_post=edt_post, rep_pre=rep_pre, rep_post=rep_post))
	r_class[get_type_pv(pv$p_1, pv$p_2)] <- r_class[get_type_pv(pv$p_1, pv$p_2)]+1	
}
summary(r_reds)
summary(c_reds)
reduction
# qqplot(r_reds, c_reds); abline(0,1)
stopifnot(t.test(r_reds, c_reds)$p.value > 0.05)

r_class/rit
c_class/cit
stopifnot(chisq.test(rbind(r_class[c(1,2,4)], c_class[c(1,2,4)]))$p.value > 0.05)




##### Paired with replicates

rep_pre <- 3
rep_post <- 4
edt_pre <- 1
edt_post <- 0.1

c2 <- bayescount:::fecrt_power_wrap(reduction, pair_type=2, preN=preN, postN=postN, animalk=anik, efficacyk=effk, prek=prek, postk=postk, premean=premean, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)

c_reds <- c2$obsred
c_class <- c2$classifications[1:4]
# Bigger discrepancy between population and sample mean reductions:
stopifnot(abs(mean(c_reds-reduction)) < 0.01)

r_reds <- numeric(rit)
r_class <- numeric(4)

stopifnot(preN==postN)
for(i in 1:rit){
	premeans <- rgamma(preN, anik, rate=anik/premean)
	postmeans <- premeans * rgamma(postN, effk, rate=effk/(1-reduction))
	
	pre <- sapply(1:preN, function(x) return(sum(rnbinom(rep_pre, prek, mu= premeans[x] / edt_pre))))
	post <- sapply(1:postN, function(x) return(sum(rnbinom(rep_post, postk, mu= postmeans[x] / edt_post))))
	
	r_reds[i] <- 1 - ((mean(post)/rep_post * edt_post) / (mean(pre)/rep_pre * edt_pre))
	
	suppressWarnings(pv <- bayescount:::fecrt_pee_wrap(pre, post, edt_pre=edt_pre, edt_post=edt_post, rep_pre=rep_pre, rep_post=rep_post))
	r_class[get_type_pv(pv$p_1, pv$p_2)] <- r_class[get_type_pv(pv$p_1, pv$p_2)]+1	
}
summary(r_reds)
summary(c_reds)
reduction
# qqplot(r_reds, c_reds); abline(0,1)
stopifnot(t.test(r_reds, c_reds)$p.value > 0.05)

r_class/rit
c_class/cit
stopifnot(chisq.test(rbind(r_class[c(1,2,4)], c_class[c(1,2,4)]))$p.value > 0.05)





##### Unpaired with pooling

rep_pre <- 1
rep_post <- 1
edt_pre <- 1
edt_post <- 0.1
poolsize_pre <- 5
poolsize_post <- 2

c2 <- bayescount:::fecrt_power_wrap(reduction, pair_type=0, preN=preN, postN=postN, premean=premean, poolsize_pre=poolsize_pre, poolsize_post=poolsize_post, animalk=anik, efficacyk=effk, prek=prek, postk=postk, rep_pre=rep_pre, rep_post=rep_post, edt_pre=edt_pre, edt_post=edt_post, iterations=cit)

c_reds <- c2$obsred
c_class <- c2$classifications[1:4]
# Bigger discrepancy between population and sample mean reductions:
stopifnot(abs(mean(c_reds-reduction)) < 0.01)

r_reds <- numeric(rit)
r_class <- numeric(4)

kpre <- combinek(anik, prek)
kpost <- combinekk(anik, effk, postk)

stopifnot(preN==postN)
for(i in 1:rit){	
	premeans <- matrix(rgamma(preN*poolsize_pre, kpre, rate= kpre/premean), ncol=poolsize_pre, nrow=preN)
	postmeans <- matrix(rgamma(postN*poolsize_post, kpost, rate= kpost/((1-reduction)*premean)), ncol=poolsize_post, nrow=postN)
	
	pre <- rpois(preN, apply(premeans, 1, mean) / edt_pre)
	post <- rpois(postN, apply(postmeans, 1, mean) / edt_post)
	
	r_reds[i] <- 1 - ((mean(post)/rep_post * edt_post) / (mean(pre)/rep_pre * edt_pre))
	
	suppressWarnings(pv <- bayescount:::fecrt_pee_wrap(pre, post, edt_pre=edt_pre, edt_post=edt_post, rep_pre=rep_pre, rep_post=rep_post))
	r_class[get_type_pv(pv$p_1, pv$p_2)] <- r_class[get_type_pv(pv$p_1, pv$p_2)]+1	
}
summary(r_reds)
summary(c_reds)
reduction
# qqplot(r_reds, c_reds); abline(0,1)
stopifnot(t.test(r_reds, c_reds)$p.value > 0.05)

r_class/rit
c_class/cit
stopifnot(chisq.test(rbind(r_class[c(1,2,4)], c_class[c(1,2,4)]))$p.value > 0.05)


### All OK to here


### Still need to add and check complex/paired pooling


}
