waavp_ci_wrap <- function(pre, post){
	
	stopifnot(length(pre)>=1)
	stopifnot(length(post)>=1)

	rv <- .C(C_waavp_ci_wrap, presum=as.integer(sum(pre)), predata=as.integer(pre), preN=as.integer(length(pre)), postsum=as.integer(sum(post)), postdata=as.integer(post), postN=as.integer(length(post)), ci_l=as.double(0), ci_u=as.double(0))
	
	return(rv)
}

dobson_ci_wrap <- function(pre, post){
	
	stopifnot(length(pre)>=1)
	stopifnot(length(post)>=1)
	
	rv <- .C(C_dobson_ci_wrap, presum=as.integer(sum(pre)), postsum=as.integer(sum(post)), lci_c=as.double(lci_c), uci_c=as.double(uci_c), dobson_priors=as.double(dobson_priors), ci_l=as.double(0), ci_u=as.double(0))
	
	return(rv)
}

conjbeta_ci_wrap <- function(pre, post){
	
	stopifnot(length(pre)>=1)
	stopifnot(length(post)>=1)
	
	rv <- .C(C_conjbeta_ci_wrap, preN=as.integer(length(pre)), presum=as.integer(sum(pre)), preK=as.double(preK), postN=as.integer(length(post)), postsum=as.integer(sum(post)), postK=as.double(postK), iters=as.integer(iters), prob_priors=as.double(prob_priors), lci_b=as.double(lci_b), uci_b=as.double(uci_b), u_ci_l=as.double(0), u_ci_u=as.double(0))
	
	return(rv)
}

asymptotic_ci_wrap <- function(pre, post){
	
	stopifnot(length(pre)==length(post))
	stopifnot(length(pre)>1)
	N <- length(pre)
	
	rv <- .C(C_asymptotic_ci_wrap, N=as.integer(N), presum=as.integer(sum(pre)), predata=as.integer(pre), postsum=as.integer(sum(post)), postdata=as.integer(post), u_ci_l=as.double(0), u_ci_u=as.double(0), p_ci_l=as.double(0), p_ci_u=as.double(0))
	
	return(rv)
}

fecrt_pee_wrap <- function(pre_data, post_data, H0_1=0.9, H0_2=0.95, edt_pre=1, edt_post=1, rep_pre=1, rep_post=1, prob_priors=c(1,1), k_change=NA, true_k=NA, delta_method=TRUE, beta_iters=10^6){
	
	stopifnot(length(pre_data)>=1)
	stopifnot(length(post_data)>=1)

	stopifnot(beta_iters <= 10^8)
	# Max 10^8 - 10^6 is pretty instant and gives the same results

	stopifnot(all(abs(round(pre_data)-pre_data) < 10^-6, na.rm=TRUE))
	stopifnot(all(abs(round(post_data)-post_data) < 10^-6, na.rm=TRUE))
	stopifnot(length(prob_priors)==2 && all(prob_priors > 0))
	
	# NB edt is a threshold, i.e opposite of replicates!
	edt_change <- (rep_post*edt_pre) / (rep_pre*edt_post)
	stopifnot(length(edt_change)==1 && !is.na(edt_change) && edt_change > 0)
	
	stopifnot(length(true_k)==1)
	stopifnot(length(k_change)==1)
	stopifnot(is.na(k_change) || k_change > 0)
	stopifnot(sum(pre_data) > 0)
	if(!is.na(true_k)){
		preK <- true_k
		if(is.na(k_change)){
			postK <- (true_k/rep_pre) * rep_post
		}else{
			postK <- (true_k/rep_pre) * rep_post * k_change
		}
	}else{
		preK <- mean(pre_data)^2 / (var(pre_data) - mean(pre_data))
		if(!is.na(k_change)){
			postK <- (preK/rep_pre) * rep_post * k_change
		}else if(sum(post_data)==0){
			postK <- (preK/rep_pre) * rep_post
		}else{
			postK <- mean(post_data)^2 / (var(post_data) - mean(post_data))
			if(postK==Inf || postK <= 0)
				postK <- (preK/rep_pre) * rep_post
		}
	}
	stopifnot(preK > 0 & postK > 0)
	
	results <- .C(C_fecrt_pee_wrap, predata=as.integer(pre_data), preN=as.integer(length(pre_data)), presum=as.integer(sum(pre_data)), preK=as.numeric(preK), postdata=as.integer(post_data), postN=as.integer(length(post_data)), postsum=as.integer(sum(post_data)), postK=as.numeric(postK), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), edt_change=as.numeric(edt_change), prob_priors=as.numeric(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0))

#	if(is.na(results$p_1) || is.na(results$p_2)) browser()
#	return(list(p_1=results$p_1, p_2=results$p_2))
	return(results)
	
}
	
