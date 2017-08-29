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

asymptotic_ci_wrap <- function(pre, post, tail=0.025){
	
	stopifnot(length(pre)==length(post))
	stopifnot(length(pre)>1)
	N <- length(pre)
	stopifnot(tail < 1)
	stopifnot(tail > 0)
	
	rv <- .C(C_asymptotic_ci_wrap, N=as.integer(N), presum=as.integer(sum(pre)), predata=as.integer(pre), postsum=as.integer(sum(post)), postdata=as.integer(post), tail=as.double(tail), u_ci_l=as.double(0), u_ci_u=as.double(0), p_ci_l=as.double(0), p_ci_u=as.double(0))
	
	return(rv)
}

fecrt_pee_wrap <- function(pre_data, post_data, H0_1=0.9, H0_2=0.95, edt_pre=1, edt_post=1, rep_pre=1, rep_post=1, pool_pre=1, pool_post=1, prob_priors=c(1,1), k_change=NA, true_k=NA, delta_method=TRUE, beta_iters=10^6){
	
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
	
	# Need to adjust H0 thresholds according to the change in mean from replicates and edt:
	H0_1 <- 1 - ((1 - H0_1) * edt_change)
	H0_2 <- 1 - ((1 - H0_2) * edt_change)
	
	stopifnot(length(true_k)==1)
	stopifnot(length(k_change)==1)
	stopifnot(is.na(k_change) || k_change > 0)
	stopifnot(sum(pre_data) > 0)
	
	## Providing fixed k with different replicates and pools needs thoroughly testing:
	if(!is.na(true_k)){
		preK <- true_k
		if(is.na(k_change)){
			postK <- (true_k/(rep_pre * pool_pre)) * rep_post * pool_post
		}else{
			postK <- (true_k/(rep_pre * pool_pre)) * rep_post * pool_post * k_change
		}
	}else{
		preK <- mean(pre_data)^2 / (var(pre_data) - mean(pre_data))
		
		# If preK can't be calculated then set to 10 with a warning:
		if(var(pre_data) < mean(pre_data)){
			warning('The mean of the pre-treatment/control data is higher than its variance: using an estimate of k=10')
			preK <- 10
		}
		if(preK > 10){
			warning('The variance of the pre-treatment/control data is very close to its variance: using an estimate of k=10')
			preK <- 10
		}
		
		postK <- mean(post_data)^2 / (var(post_data) - mean(post_data))

		# If k_change is set always use it and ignore above calculation:
		if(!is.na(k_change)){
			postK <- (preK/(rep_pre * pool_pre)) * rep_post * pool_post * k_change
		
		# Otherwise if postK can't be calculated:
		}else if(sum(post_data)==0 || var(post_data) < mean(post_data) || postK > 10){
			postK <- (preK/(rep_pre * pool_pre)) * rep_post * pool_post
		}
	}
	stopifnot(preK > 0 & postK > 0)
	
	results <- .C(C_fecrt_pee_wrap, presum=as.integer(sum(pre_data)), preN=as.integer(length(pre_data)), preK=as.numeric(preK), postsum=as.integer(sum(post_data)), postN=as.integer(length(post_data)),  postK=as.numeric(postK), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), prob_priors=as.numeric(prob_priors), delta=as.integer(delta_method), beta_iters=as.integer(beta_iters), p_1=as.numeric(0), p_2=as.numeric(0))
	
#	if(is.na(results$p_1) || is.na(results$p_2)) browser()
#	return(list(p_1=results$p_1, p_2=results$p_2))
	return(results)
	
}


fecrt_pee_direct <- function(presum, preN, preK, postsum, postN, postK, H0_1, H0_2, prob_priors=c(1,1)){
	
	results <- .C(C_fecrt_pee_wrap, presum=as.integer(presum), preN=as.integer(preN), preK=as.numeric(preK), postsum=as.integer(postsum), postN=as.integer(postN),  postK=as.numeric(postK), H0_1=as.numeric(H0_1), H0_2=as.numeric(H0_2), prob_priors=as.numeric(prob_priors), delta=as.integer(1), beta_iters=as.integer(0), p_1=as.numeric(0), p_2=as.numeric(0))
	return(results)
	
}	

# There are 4 versions of this function:
	# Main wrapper - allows anything, pair type specified manually, arguments change meaning depending on pair type
	# Others call main wrapper:
		# Simple unpaired - 1 rep only, no pooling, pair_type 0
		# Complex unpaired - either reps or pooling (not both?), pair_type 1
		# Paired - 1+ reps, no pooling (but add pooling??), pair_type 2 or 3 depending on post reps (or always 3??)
	
# TODO: add pooling (just adjust k sources as mean unchanged)
	
fecrt_power_wrap <- function(reduction, pair_type=0, preN=20, postN=20, poolsize_pre=1, poolsize_post=1, premean=10, animalk=2, efficacyk=2, prek=3, postk=prek, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, target=0.95, tol=0.05, tail=0.025, iterations=10000, prob_priors=c(1,1)){

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(pair_type)==1 && pair_type >= 0)
	stopifnot(length(preN)==1 && preN > 0 && round(preN)==preN)
	stopifnot(length(postN)==1 && postN > 0 && round(postN)==postN)
	stopifnot(length(poolsize_pre)==1 && poolsize_pre > 0 && round(poolsize_pre)==poolsize_pre)
	stopifnot(length(poolsize_post)==1 && poolsize_post > 0 && round(poolsize_post)==poolsize_post)
	stopifnot(length(premean)==1 && premean > 0)
	stopifnot(length(animalk)==1 && animalk > 0)
	stopifnot(length(prek)==1 && prek > 0)
	stopifnot(length(postk)==1 && postk > 0)
	stopifnot(length(efficacyk)==1 && efficacyk > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_pre)==1 && edt_pre > 0)
	stopifnot(length(edt_post)==1 && edt_post > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))
	
	H0_1 <- target-tol
	H0_2 <- target
	reduction 

	# The paired int is either:
	# 0:  no pairing, no animalk, no efficacyk, can't use replicates or pooling
	# 1:  no pairing, use animalk and efficacyk so we can use replicates or pooling
	# 2:  with pairing, use animalk and efficacyk so we can use replicates or pooling
	# suspended:  with pairing, use animalk but no efficacyk, can't use replicates or pooling - disadvantage is distribution change from combining gammas
	
	if(pair_type == 0){

		stopifnot(rep_pre==1 && rep_post==1)
		
		# C++ function does not combine variances so need to do this if animalk AND efficacyk not supplied as Inf:
		if(animalk == Inf || efficacyk == Inf){
			stopifnot(animalk==Inf && efficacyk==Inf)

			# Avoid NA/Inf in foreign call:
			animalk <- 0
			efficacyk <- 0			
		}else{
			prek <- (animalk * prek) / (animalk + prek + 1)
			intk <- (animalk * efficacyk) / (animalk + efficacyk + 1)
			postk <- (intk * postk) / (intk + postk + 1)
			# Probably need a combine_k function!
		}
		
	}else if(pair_type == 1){
		
		# C++ function combines animal and efficacy variability sources for treatment group

	}else if(pair_type == 2){
		
		# All variability sources are separate
		
	}else if(pair_type == 3){
		
		stopifnot(rep_pre==1 && rep_post==1)
		
		# C++ function combines efficacyk and postk
		stop('pair_type is disabled')
		
	}else{
		stop('Unrecognised pair_type')
	}
	
	# Pooling works by reducing k sources:
	if(poolsize_pre > 1 || poolsize_post > 1){
		
		# NOTE that in the case of pooling, preN and postN becomes number of pools!
		
		# Currently only available for simple unpaired:
		stopifnot(pair_type == 0)
		
		prek <- prek * poolsize_pre
		postk <- postk * poolsize_post
		
		# When taking into account rep_pre and rep_post will need to adjust k and then set rep to 1
		# (assumes that repeat samples from the same animal are always pooled together)
		# Alternatively, rep_pre and rep_post could represent repeat samples from the same pool...
		
	}
	
	# Currently only 5*2 but will be 5*4:
	classifications <- 20

	results <- .C(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre),
	rep_post=as.integer(rep_post), edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), pair_type=as.integer(pair_type),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail),
	prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(preN, postN))), postdata=as.integer(rep(0, max(preN, postN))), classifications=as.integer(rep(0, classifications)),
	obsred=as.double(rep(0, iterations)))

	return(results)

}

fecrt_power_paired <- function(reduction, preN=20, postN=20, animalk=2, efficacyk=2, prek=3, postk=prek, premean=10, rep_pre=1, rep_post=1, edt_pre=1, edt_post=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){
	
	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(preN)==1 && preN > 0 && round(preN)==preN)
	stopifnot(length(postN)==1 && postN > 0 && round(postN)==postN)
	stopifnot(length(animalk)==1 && animalk > 0)
	stopifnot(length(prek)==1 && prek > 0)
	stopifnot(length(postk)==1 && postk > 0)
	stopifnot(length(efficacyk)==1 && efficacyk > 0)
	stopifnot(length(premean)==1 && premean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_pre)==1 && edt_pre > 0)
	stopifnot(length(edt_post)==1 && edt_post > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))
	
	H0_1 <- target-tol
	H0_2 <- target
	reduction 

	# The paired int is either:
	# 0:  no pairing, no animalk, no efficacyk, can't use replicates or pooling
	# 1:  no pairing, use animalk and efficacyk so we can use replicates or pooling
	# 2:  with pairing, use animalk but no efficacyk, can't use replicates or pooling
	# 3:  with pairing, use animalk and efficacyk so we can use replicates or pooling

	if(rep_pre > 1 || rep_post > 1){
		# The C++ function will combine efficacy and postk
		paired <- 2
	}else{
		paired <- 3
	}

	
	# Not using but keep for info:
	if(FALSE && paired==0){
		prek <- (animalk * prek) / (animalk + prek + 1)
		intk <- (animalk * efficacyk) / (animalk + efficacyk + 1)
		postk <- (intk * postk) / (intk + postk + 1)
	}
	
	# Currently only 5*2 but will be 5*4:
	classifications <- 20
	
	results <- .C(C_fecrt_power_comparison, iters=as.integer(iterations), preN=as.integer(preN), postN=as.integer(postN), maxN=as.integer(max(preN, postN)), rep_pre=as.integer(rep_pre),
	rep_post=as.integer(rep_post), edt_pre=as.double(edt_pre), edt_post=as.double(edt_post), premean=as.double(premean), reduction=as.double(reduction), paired=as.integer(paired),
	animalk=as.double(animalk), efficacyk=as.double(efficacyk), prek=as.double(prek), postk=as.double(postk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail),
	prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(preN, postN))), postdata=as.integer(rep(0, max(preN, postN))), classifications=as.integer(rep(0, classifications)),
	obsred=as.double(rep(0, iterations)))
	
	return(results)
	
}

fecrt_power_unpaired <- function(reduction, controlN=20, treatmentN=20, controlk=1, treatmentk=0.5, controlmean=10, edt_control=1, edt_treatment=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){
	
	paired <- FALSE
	# No replicates for simple unpaired:

	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(paired)==1 && is.logical(paired))
	stopifnot(length(controlN)==1 && controlN > 0 && round(controlN)==controlN)
	stopifnot(length(treatmentN)==1 && treatmentN > 0 && round(treatmentN)==treatmentN)
	stopifnot(length(controlk)==1 && controlk > 0)
	stopifnot(length(treatmentk)==1 && treatmentk > 0)
	stopifnot(length(controlmean)==1 && controlmean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_control)==1 && edt_control > 0)
	stopifnot(length(edt_treatment)==1 && edt_treatment > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))
	
	H0_1 <- target-tol
	H0_2 <- target
	
	# Currently only 5*2 but will be 5*4:
	classifications <- 20
	
	results <- .C(C_fecrt_power_comparison, iters=as.integer(iterations), controlN=as.integer(controlN), treatmentN=as.integer(treatmentN), maxN=as.integer(max(controlN, treatmentN)), rep_pre=as.integer(rep_pre), rep_post=as.integer(rep_post), edt_control=as.double(edt_control), edt_treatment=as.double(edt_treatment), controlmean=as.double(controlmean), reduction=as.double(reduction), paired=as.integer(paired), animalk=as.double(Inf), efficacyk=as.double(Inf), prek=as.double(controlk), postk=as.double(treatmentk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail), prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(controlN, treatmentN))), postdata=as.integer(rep(0, max(controlN, treatmentN))), classifications=as.integer(rep(0, classifications)), obsred=as.double(rep(0, iterations)))
	
	return(results)
	
}

fecrt_power_pooled <- function(reduction, pools_control=5, pools_treatment=5, Nperpool_control=4, Nperpool_treatment=4, controlk=1, treatmentk=0.5, controlmean=10, edt_control=1, edt_treatment=1, target=0.95, tol=0.05, tail=0.025, iterations=1000, prob_priors=c(1,1)){
	
	# Pools is for unpaired data only currently
	# (but can be extended to unpaired just by adjusting all k values ???  Check by simulation)
	paired <- FALSE
	
	# No replicates for pools currently
	# (although theoretically justifiable as long as all replicates from the same animal in)
	rep_pre <- 1
	rep_post <- 1
	
	# 
	
	stopifnot(length(reduction)==1 && reduction <= 1 && reduction > -Inf)
	stopifnot(length(paired)==1 && is.logical(paired))
	stopifnot(length(controlN)==1 && controlN > 0 && round(controlN)==controlN)
	stopifnot(length(treatmentN)==1 && treatmentN > 0 && round(treatmentN)==treatmentN)
	stopifnot(length(controlk)==1 && controlk > 0)
	stopifnot(length(treatmentk)==1 && treatmentk > 0)
	stopifnot(length(controlmean)==1 && controlmean > 0)
	stopifnot(length(rep_pre)==1 && rep_pre > 0 && round(rep_pre)==rep_pre)
	stopifnot(length(rep_post)==1 && rep_post > 0 && round(rep_post)==rep_post)
	stopifnot(length(edt_control)==1 && edt_control > 0)
	stopifnot(length(edt_treatment)==1 && edt_treatment > 0)
	stopifnot(length(target)==1 && target <= 1 && target > -Inf)
	stopifnot(length(tol)==1 && tol < 1)
	stopifnot(length(tail)==1 && tail > 0 && tail < 1)
	stopifnot(length(iterations)==1 && iterations > 0 && round(iterations)==iterations)
	stopifnot(length(prob_priors)==2 && is.numeric(prob_priors) && all(prob_priors > 0))
	
	H0_1 <- target-tol
	H0_2 <- target
	
	# Currently only 5*2 but will be 5*4:
	classifications <- 20
	
	results <- .C(C_fecrt_power_comparison, iters=as.integer(iterations), controlN=as.integer(controlN), treatmentN=as.integer(treatmentN), maxN=as.integer(max(controlN, treatmentN)), rep_pre=as.integer(rep_pre), rep_post=as.integer(rep_post), edt_control=as.double(edt_control), edt_treatment=as.double(edt_treatment), controlmean=as.double(controlmean), reduction=as.double(reduction), paired=as.integer(paired), animalk=as.double(Inf), efficacyk=as.double(Inf), prek=as.double(controlk), postk=as.double(treatmentk), H0_1=as.double(H0_1), H0_2=as.double(H0_2), tail=as.double(tail), prob_priors=as.double(prob_priors), predata=as.integer(rep(0, max(controlN, treatmentN))), postdata=as.integer(rep(0, max(controlN, treatmentN))), classifications=as.integer(rep(0, classifications)), obsred=as.double(rep(0, iterations)))
	
	return(results)
	
}
