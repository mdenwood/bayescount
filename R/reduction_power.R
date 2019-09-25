#' @name reduction_power
#' @aliases reduction.power fecrt.power FECRT.power
#' @title Count reduction (eg FECRT) data power calculations
#'
#' @description
#' Power calculation based on BNB code and Christian approximation.
#'
#' @seealso \code{\link{count_analysis}}, \code{\link{reduction_power}}, \code{\link{bayescount}}


reduction_pval <- function(pre_data, post_data, H0_1=0.9, H0_2=0.95, edt_pre=1, edt_post=1, prob_priors=c(1,1), k_change=1, true_k=NA, delta_method=TRUE, beta_iters=10^6){
	
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results
	
	stopifnot(all(abs(round(pre_data)-pre_data) < 10^-6, na.rm=TRUE))
	stopifnot(all(abs(round(post_data)-post_data) < 10^-6, na.rm=TRUE))
	stopifnot(length(prob_priors)==2 && all(prob_priors > 0))
	
	edt_change <- edt_pre / edt_post
	stopifnot(length(edt_change)==1 && !is.na(edt_change) && edt_change > 0)

	stopifnot(length(k_change)==1 && !is.na(k_change) && k_change > 0)
	
	usetruek <- TRUE
	if(all(is.na(true_k))){
		true_k <- 1
		usetruek <- FALSE
	}
	stopifnot(length(true_k)==1 && !is.na(true_k) && true_k > 0)
	
	results <- .C_OLD(C_fecrt_pee_wrap, pre=as.integer(pre_data), preN=as.integer(length(pre_data)), post=as.integer(post_data), postN=as.integer(length(post_data)), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), edt_change=as.numeric(edt_change), prob_priors=as.numeric(prob_priors), kchange=as.numeric(k_change), truek=as.numeric(true_k), usetruek=as.integer(usetruek), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0))
	
	return(list(p_1=results$p_1, p_2=results$p_2))
	
}


reduction_pvals <- function(N, pre_mean, reduction, pre_k=1, post_k=pre_k, paired=FALSE, iterations=10000, animal_k=1, H0_1=0.9, H0_2=0.95, preN=N, postN=N, edt_pre=1, edt_post=1, prob_priors=c(1,1), k_change=1, true_k=NA, delta_method=TRUE, beta_iters=10^6){
	
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results
	stopifnot(iterations <= 10^8)

	usetruek <- TRUE
	if(all(is.na(true_k))){
		true_k <- 1
		usetruek <- FALSE
	}
	stopifnot(length(true_k)==1 && !is.na(true_k) && true_k > 0)

	if(paired){
		stopifnot(preN==postN)
	}else{
		animal_k <- 1000  # Can be >100 -> means unpaired
	}
	edt_change <- edt_pre / edt_post
	
	stopifnot(length(H0_1)==length(H0_2))
	H0_N <- length(H0_1)
	maxN <- max(preN, postN)

	p_1=p_2 <- rep(0, iterations*H0_N)
	
	pre <- rep(0, maxN)
	post <- rep(0, maxN)
	
	powres <- .C_OLD(C_fecrt_pvals, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(maxN), premean=as.numeric(pre_mean), reduction=as.numeric(reduction), edt_change=as.numeric(edt_change), animalk=as.numeric(animal_k), prek=as.numeric(pre_k), postk=as.numeric(post_k), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), H0_N=as.integer(H0_N), prob_priors=as.numeric(prob_priors), kchange=as.numeric(k_change), truek=as.numeric(true_k), usetruek=as.integer(usetruek), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), predata=as.integer(pre), postdata=as.integer(post), p_1=as.numeric(p_1), p_2=as.numeric(p_2))

	p_1 <- powres$p_1
	p_2 <- powres$p_2
	dim(p_1) <- c(iterations, H0_N)
	dim(p_2) <- c(iterations, H0_N)

	return(data.frame(p1=p_1, p2=p_2))
	
}


reduction_power <- function(N, pre_mean, reduction, pre_k=1, post_k=pre_k, paired=FALSE, iterations=10000, animal_k=1, H0_1=0.9, H0_2=0.95, preN=N, postN=N, edt_pre=1, edt_post=1, prob_priors=c(1,1), k_change=1, true_k=NA, delta_method=TRUE, beta_iters=10, lci_c=0.005, uci_c=0.995, lci_b=0.025, uci_b=0.975, dobson_priors=c(1,1), tval=2.048, psig=0.025){
	
	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results
	stopifnot(iterations <= 10^8)
	
	edt_change <- edt_pre / edt_post
	stopifnot(length(edt_change)==1 && !is.na(edt_change) && edt_change > 0)

	stopifnot(length(k_change)==1 && !is.na(k_change) && k_change > 0)
	
	usetruek <- TRUE
	if(all(is.na(true_k))){
		true_k <- 1
		usetruek <- FALSE
	}
	stopifnot(length(true_k)==1 && !is.na(true_k) && true_k > 0)

	if(paired){
		stopifnot(preN==postN)
	}else{
		animal_k <- 1000  # Can be >100 -> means unpaired
	}
	
	# Needs other checks

	maxN <- max(preN, postN)
	predata=postdata <- rep(0, maxN)
	classifications <- rep(0, 3*3*3*3*2)	
	
	results <- .C_OLD(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(maxN), premean=as.numeric(pre_mean), reduction=as.numeric(reduction), edt_change=as.numeric(edt_change), animalk=as.numeric(animal_k), prek=as.numeric(pre_k), postk=as.numeric(post_k), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), lci_c=as.numeric(lci_c), uci_c=as.numeric(uci_c), lci_b=as.numeric(lci_b), uci_b=as.numeric(uci_b), dobson_priors=as.numeric(dobson_priors), tval=as.numeric(tval), psig=as.numeric(psig), prob_priors=as.numeric(prob_priors), kchange=as.numeric(k_change), truek=as.numeric(true_k), usetruek=as.integer(usetruek), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), predata=as.integer(predata), postdata=as.integer(postdata), classifications=as.integer(classifications))

	classifications <- results$classifications
	dim(classifications) <- c(3,3,3,3,2)
	dimnames(classifications) <- list(paste('WAAVP', c('Susceptible','Inconclusive','Resistant'), sep='_'), paste('Dobson', c('Susceptible','Inconclusive','Resistant'), sep='_'), paste('Beta', c('Susceptible','Inconclusive','Resistant'), sep='_'), paste('Delta', c('Susceptible','Inconclusive','Resistant'), sep='_'), c('100pcObs', 'Lt100pcObs'))
	
	return(classifications)	
}


reduction.power <- reduction_power
fecrt.power <- reduction_power
FECRT.power <- reduction_power
