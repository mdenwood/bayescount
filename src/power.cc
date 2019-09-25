/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for power calculations for FECRT studies
Last updated by MJD June 2018
*/

/*****************************************/
/*   C wrappers to be exported to R      */
/*****************************************/

#include <math.h>
#include <R.h>
#include <Rmath.h>

#include "fecrt.h"

extern "C" {
		
	#include "power.h"
	
	/*************************************************/
	/*   Power analysis wrapper (obtain p-values)    */
	/*************************************************/
	
	/*
	Note to self:
	fecrt_power and _comparison need to save the two means and variances and pass them to WAAVP, asymptotic and BNB methods
	None of these underlying functions will then need the raw data so the relevant summary stats can be calculated in the calling function/R and arrays of data removed from the argument list
	*/
	
	void fecrt_power(int *iters, int *paired_analysis, int *use_truek, double *kchange, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *conjugate_priors, int *delta, int *beta_iters, int *approximation, double *p_1, double *p_2, double *obsred, int *ndeltafail){

		// Note: pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// If not paired model then rep_pre and rep_post must be 1
		// Pooling must be accounted for in R
		// Note in the paper for asymptotic methods that p-values or CI can be extracted
		
		GetRNGstate();
		
		int a, h;
		int ndfail = 0;
		long in;
		double anmean1, anmean2, est_prek, est_postk, estpremu, estprevar, estpostmu, estpostvar, estcov;
		
		// Resolve pointers (and some casts) once:
		int it = *iters;
		int delta_method = *delta;  // Note: delta_method can take 3 values:  0=never, 1=unless_fails, 2=always but should probably be 1 or 2 here!
		int beta_it = *beta_iters;  // Probably don't want too many iterations if the delta method fails!
		int approx = *approximation;  	// Note: approx can take 3 values:  0=never, 1=if_necessary, 2=always
		bool pan = (bool) *paired_analysis;
		bool usetruek = *use_truek;
		double k_change = *kchange;
		int pair_type = *pairtype;
		int pre_N = *preN;
		int post_N = *postN;
		int max_N = *maxN;
		int H0N = *H0_N;
		double pre_Nd = (double) pre_N;
		double post_Nd = (double) post_N;
		double pre_rep_d = (double) *rep_pre;
		double post_rep_d = (double) *rep_post;
		double red_val = 1.0 - *reduction;
		double animal_k = *animalk;
		double efficacy_k = *efficacyk;
		double pre_k = *prek;
		double post_k = *postk;
		// Only needed for unpaired model with replicates:
		double treatment_k = (animal_k * efficacy_k) / (animal_k + efficacy_k + 1);   // See solve 1/a = 1/b+1/c+1/(b*c) for a in WA
		// Only needed for paired model without replicates:
		// double postcombined_k = (efficacy_k * post_k) / (efficacy_k + post_k + 1);
		// double eff_postcombined_k = postcombined_k * pre_rep_d;
		
		// Take account of edt_change parameters:
		double eff_pre_mu = *premean / (double) *edt_pre;
		double eff_post_mu = (*premean * red_val) / (double) *edt_post;
		double eff_red_val = red_val * ((double) *edt_pre / (double) *edt_post);
		double eff_pre_k = pre_k * pre_rep_d;
		double eff_post_k = post_k * post_rep_d;
		double mean_ratio = (post_rep_d * *edt_pre) / (pre_rep_d * *edt_post);
		
		// This correction is only used for the hypothesis test thresholds and calculating the obs reduction:
		double lts[H0N];
		double uts[H0N];
		double correction = (pre_rep_d * (double) *edt_pre) / (post_rep_d * (double) *edt_post);
		for(h=0; h<H0N; h++){
//			lts[h] = 1.0 - ((1.0 - H0_1[h]) * correction);
//			uts[h] = 1.0 - ((1.0 - H0_2[h]) * correction);
		}
		double pre_post_ratio = (pre_Nd / post_Nd) * (1.0 / correction);
		
		// Use long types for simulated counts:
		int predata[max_N];
		int postdata[max_N];
		long long presum, postsum;
		
		// If using the paired analysis then preN must equal postN:
		if(pan && (pre_N!=post_N)){
			// Error handling:  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Error-signaling
			error("preN must equal postN when using the paired analysis");
		}
		if(pair_type==2 && (pre_N!=post_N)){
			// Error handling:  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Error-signaling
			error("preN must equal postN when using pair_type 2");
		}
		
		// Calculate true pre and post tx ks:
		double prek_true, postk_true, cov_true;		
		// Simple unpaired model ignores animalk and efficacyk:			
		if(!usetruek || pair_type == 0){
			prek_true = pre_k;
			postk_true = post_k;
			cov_true = 0.0;
		// Unpaired model but with replicates needs to use animalk and efficacyk:
		}else if(pair_type == 1){
			// Combination of ks:
			double repk = eff_pre_k * pre_rep_d;
			prek_true = (animal_k * repk) / (animal_k + repk + 1.0);
			repk = eff_post_k * post_rep_d;
			postk_true = (treatment_k * repk) / (treatment_k + repk + 1.0);
			cov_true = 0.0;
		// If it is a paired model postk_true and cov_true are more complex:
		}else if(pair_type == 2){
			double repk = eff_pre_k * pre_rep_d;
			prek_true = (animal_k * repk) / (animal_k + repk + 1.0);
			repk = eff_post_k * post_rep_d;
			postk_true = (treatment_k * repk) / (treatment_k + repk + 1.0);
			
			// Estimate true covariance empirically:
			cov_true = 0.0;
			long long bigN = 1000000;
			double tpre, tpost;
		    for(a=0; a<bigN; a++){
				anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
				anmean2 = anmean1 * rgamma(efficacy_k, (eff_red_val / efficacy_k));
				tpre = rgamma(eff_pre_k, (anmean1*pre_rep_d) / eff_pre_k);
				tpost = rgamma(eff_post_k, (anmean2*post_rep_d) / eff_post_k);
				cov_true += (tpre - eff_pre_mu) * (tpost - eff_post_mu);
			}
		    cov_true = cov_true / (bigN - 1.0);
		
		// pair_type should match one of the above but have this for safety:
		}else{
			// Error handling:  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Error-signaling
			error("Unmatched pair_type");				
		}
		
		// First loop over iterations:
		int i=0;
		while(i < it){
			// Then simulate data using array allocated in R:
			
			presum = 0;
			postsum = 0;
			
			// Simple unpaired model ignores animalk and efficacyk:			
			if(pair_type == 0){
				// We only need pre and post data where required:
				for(a=0; a<pre_N; a++){
					 // temp = rgamma(pre_k, (eff_pre_mu / pre_k));
					 // predata[a] = rpois(temp);
					predata[a] = rnbinom_mu(pre_k, eff_pre_mu);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					// temp = rgamma(post_k, (eff_post_mu / post_k));
					// postdata[a] = rpois(temp);
					postdata[a] = rnbinom_mu(post_k, eff_post_mu);
					postsum += postdata[a];
				}

			// Unpaired model but with replicates needs to use animalk and efficacyk:
			}else if(pair_type == 1){
				for(a=0; a<pre_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					// temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					// predata[a] = rpois(temp);
					predata[a] = rnbinom_mu(eff_pre_k, anmean1*pre_rep_d);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					anmean2 = rgamma(treatment_k, (eff_post_mu / treatment_k));
					// temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					// postdata[a] = rpois(temp);
					postdata[a] = rnbinom_mu(eff_post_k, anmean2*post_rep_d);
					postsum += postdata[a];
				}
			
			// If it is a paired model we need pre and post means for all animals:
			// (don't need all poisson variates but do them to avoid another loop
			// or conditional, as preN nearly always = postN for paired anyway):
			}else if(pair_type == 2){
				for(a=0; a<max_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					anmean2 = anmean1 * rgamma(efficacy_k, (eff_red_val / efficacy_k));
					// temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					// predata[a] = rpois(temp);
					// temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					// postdata[a] = rpois(temp);
					predata[a] = rnbinom_mu(eff_pre_k, anmean1*pre_rep_d);
					postdata[a] = rnbinom_mu(eff_post_k, anmean2*post_rep_d);
				}
				for(a=0; a<pre_N; a++){
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					postsum += postdata[a];
				}				
			
			/* Suspend this	code: little is gained by avoiding a gamma draw, and power change from
				1->2 replicates is unpredictable because the distribution changes
				
			// If not using replicates then use combined post and efficacy k:
			}else if(pair_type == 3){
				for(a=0; a<max_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					anmean2 = eff_red_val * anmean1;
					temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					predata[a] = rpois(temp);
					temp = rgamma(eff_postcombined_k, (anmean2 / postcombined_k));  // replicates cancels out of scale parameter
					postdata[a] = rpois(temp);
				}
				for(a=0; a<pre_N; a++){
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					postsum += postdata[a];
				}				

			*/

			// pair_type should match one of the above but have this for safety:
			}else{
				// Error handling:  https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Error-signaling
				error("Unmatched pair_type");				
			}

			// Rprintf("Pre sum: %i, Post sum: %i, red_val: %f\n", presum, postsum, red_val);
			
			// If the presum is 0 try to simulate data again:
			if(presum==0){
				continue;
			}
			
			// Save observed reduction:
			obsred[i] = 1.0 - (((double) postsum) / ((double) presum) * pre_post_ratio);
			
			estpremu = (double) presum / pre_Nd;
			estpostmu = (double) postsum / post_Nd;

			// Estimate pre and post k
			estprevar = 0.0;
			for(a=0; a<pre_N; a++){
				estprevar += pow(((double) predata[a] - estpremu), 2);
			}
			estprevar = estprevar / (pre_Nd - 1.0);
			// If the estimated variance is greater than the mean:
			if(estprevar > estpremu){
				est_prek = pow(estpremu, 2) / (estprevar - estpremu);
				// If either variance < mean or close then use default of 10:
				if(est_prek > 10.0){
					est_prek = 10.0;
				}
			}else{
				est_prek = 10.0;
			}
		
			estpostvar = 0.0;
			for(a=0; a<post_N; a++){
				estpostvar += pow(((double) postdata[a] - estpostmu), 2);
			}
			estpostvar = estpostvar / (post_Nd - 1.0);
			// If the estimated variance is greater than the mean:
			if(estpostvar > estpostmu){
				est_postk = pow(estpostmu, 2) / (estpostvar - estpostmu);
				// If either variance < mean or close then use default of prek:
				if(est_postk > 10.0){
					est_postk = est_prek * k_change;
				}
			}else{
				est_postk = est_prek * k_change;
			}
		
			// If using the paired analysis (and pre and post variance both >0), estimate the covariance:
			estcov = 0.0;
			if(pan && estprevar > 0 && estpostvar > 0){
				// Note: pre_N and post_N are the same for paired:
			    for(a=0; a<max_N; a++){
					estcov += ((double) predata[a] - estpremu) * ((double) postdata[a] - estpostmu);
				}
			    estcov = estcov / (max_N - 1.0);
			}
			
			if(usetruek){
				est_prek = prek_true;
				est_postk = postk_true;
				estcov = cov_true;
				estprevar = estpremu + pow(estpremu, 2)/prek_true;
				if(postsum==0){
					estpostvar = 0.0;
				}else{
					estpostvar = estpostmu + pow(estpostmu, 2)/postk_true;
				}
			}
			
			/*
			printf("pre = c(");
			for(a=0; a<max_N; a++){
				printf("%i", predata[a]);
				if(a != (max_N-1)){
					printf(", ");
				}
			}
			printf("; post = c(");
			for(a=0; a<max_N; a++){
				printf("%i", postdata[a]);
				if(a != (max_N-1)){
					printf(", ");
				}
			}
			printf(")\nestpremu = %f, estpostmu = %f, est_prek = %f, est_postk = %f\n", estpremu, estpostmu, est_prek, est_postk);
			*/
			
			// Then loop over H0 values:
			for(h=0; h<H0N; h++){
				// The array index to save to:
				in = h*it + i;
				
				// Pass pointers to correct index of p_1 and p_2:
				ndfail += bnb_pval(presum, pre_N, est_prek, estpremu, estprevar, postsum, post_N, est_postk, estpostmu, estpostvar, estcov, mean_ratio, H0_1[h], H0_2[h], conjugate_priors, delta_method, beta_it, approx, &p_1[in], &p_2[in]);		
						
				// This is for compensated H0/H1 but I think not necessary - test reps before deleting:
//				ndfail += bnb_pval(presum, pre_N, est_prek, estpremu, estprevar, postsum, post_N, est_postk, estpostmu, estpostvar, estcov, mean_ratio, lts[h], uts[h], conjugate_priors, delta_method, beta_it, approx, &p_1[in], &p_2[in]);				
				
				// int bnb_pval(long long sum1, int N1, double K1, double mu1, double var1, long long sum2, int N2, double K2, double mu2, double var2, double cov12, 
				// double mean_ratio, double H0_1, double H0_2, double *conjugate_priors, int delta, int beta_iters, int approx, double *p_1, double *p_2);
				// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
				
			}
			
			i++;
		}
		
		ndeltafail[0] = ndfail;
		PutRNGstate();
	}
	

	/*************************************************/
	/*   Power comparison for different tests        */
	/*************************************************/
	
	void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, double *tail, double *conjugate_priors, int *delta, int *beta_iters, int poscov, double *lci_c, double *uci_c, double *dobson_priors, int *classifications, double *obsred, int *ndeltafail){
		
		// Note: pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// If not paired model then rep_pre and rep_post must be 1
		// Pooling must be accounted for in R
		
		GetRNGstate();
		
		int a, r;
		int ndfail = 0;
		double anmean1, anmean2, est_prek, est_postk, estmu, estvar;
		double p1, p2, w1, w2, w3, w4;

		// Resolve pointers (and some casts) once:
		int delta_method = delta[0];  // Note: delta_method can take 3 values:  0=never, 1=unless_fails, 2=always but should probably be 1 or 2 here!
		int beta_it = beta_iters[0];  // Probably don't want too many iterations if the delta method fails!
		int pair_type = pairtype[0];
		int pre_N = preN[0];
		int post_N = postN[0];
		int max_N = maxN[0];
		double pre_Nd = (double) pre_N;
		double post_Nd = (double) post_N;
		double pre_rep_d = (double) rep_pre[0];
		double post_rep_d = (double) rep_post[0];
		double pt = tail[0];
		double red_val = 1.0 - reduction[0];
		double animal_k = animalk[0];
		double efficacy_k = efficacyk[0];
		double pre_k = prek[0];
		double post_k = postk[0];
		// Only needed for unpaired model with replicates:
		double treatment_k = (animal_k * efficacy_k) / (animal_k + efficacy_k + 1);
		// Only needed for paired model without replicates:
		// double postcombined_k = (efficacy_k * post_k) / (efficacy_k + post_k + 1);
		// double eff_postcombined_k = postcombined_k * pre_rep_d;
		
		// Take account of edt_change parameters:
		double eff_pre_mu = premean[0] / (double) edt_pre[0];
		double eff_post_mu = (premean[0] * red_val) / (double) edt_post[0];
		double eff_red_val = red_val * ((double) edt_pre[0] / (double) edt_post[0]);
		double eff_pre_k = pre_k * pre_rep_d;
		double eff_post_k = post_k * post_rep_d;
		
		// This correction is only used for the hypothesis test thresholds and calculating the obs reduction:
		double correction = (pre_rep_d * (double) edt_pre[0]) / (post_rep_d * (double) edt_post[0]);
		double lt = 1.0 - ((1.0 - H0_1[0]) * correction);
		double ut = 1.0 - ((1.0 - H0_2[0]) * correction);
		double pre_post_ratio = (pre_Nd / post_Nd) * (1.0 / correction);
		
		// Use long types for sums of simulated counts:
		long predata[max_N];
		long postdata[max_N];
		long long presum, postsum;
		
		/*
	double premu = (double) presum / (double) preN;
	double prevarsum = 0;
	for(int a=0; a<preN; a++){
		prevarsum += std::pow((predata[a]-premu), 2);
	}
	double postmu = (double) postsum / (double) postN;
	double postvarsum = 0;
	for(int a=0; a<postN; a++){
		postvarsum += std::pow((postdata[a]-postmu), 2);
	}
	double varpre = 0;
	if(prevarsum > 0.00000001){
		varpre = prevarsum / (double (preN - 1));
	}
	double varpost = 0;
	if(postvarsum > 0.00000001){
		varpost = postvarsum / (double (postN - 1));
	}
	
	double prevar = 0.0;
	double postvar = 0.0;
	double covar = 0.0;
	for(int a=0; a<N; a++){
		prevar += pow((predata[a]-premu), 2);
		postvar += pow((postdata[a]-postmu), 2);
		covar += (predata[a]-premu) * (postdata[a]-postmu);
	}
	prevar = prevar / (Nd-1.0);
	postvar = postvar / (Nd-1.0);
	covar = covar / (Nd-1.0);
	// Option to limit the covariance to biologically sensible values:
	if(poscov && covar < 0.0){
		covar = 0.0;
	}
	
	
		*/
				
		// First loop over iterations:
		int i = 0;
		while(i < iters[0]){
			// Then simulate data using array allocated in R:
			
			presum = 0;
			postsum = 0;
			
			// Simple unpaired model ignores animalk and efficacyk:			
			if(pair_type == 0){
				// We only need pre and post data where required:
				for(a=0; a<pre_N; a++){
					// temp = rgamma(pre_k, (eff_pre_mu / pre_k));
					// predata[a] = rpois(temp);
					predata[a] = rnbinom_mu(pre_k, eff_pre_mu);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					// temp = rgamma(post_k, (eff_post_mu / post_k));
					// postdata[a] = rpois(temp);
					postdata[a] = rnbinom_mu(post_k, eff_post_mu);
					postsum += postdata[a];
				}

			// Unpaired model but with replicates needs to use animalk and efficacyk:
			}else if(pair_type == 1){
				for(a=0; a<pre_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					// temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					// predata[a] = rpois(temp);
					predata[a] = rnbinom_mu(eff_pre_k, anmean1*pre_rep_d);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					anmean2 = rgamma(treatment_k, (eff_post_mu / treatment_k));
					// temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					// postdata[a] = rpois(temp);
					postdata[a] = rnbinom_mu(eff_post_k, anmean2*post_rep_d);
					postsum += postdata[a];
				}
			
			// If it is a paired model we need pre and post means for all animals:
			// (don't need all poisson variates but do them to avoid another loop
			// or conditional, as preN nearly always = postN for paired anyway):
			}else if(pair_type == 2){
				for(a=0; a<max_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					anmean2 = anmean1 * rgamma(efficacy_k, (eff_red_val / efficacy_k));
					// temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					// predata[a] = rpois(temp);
					// temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					// postdata[a] = rpois(temp);
					predata[a] = rnbinom_mu(eff_pre_k, anmean1*pre_rep_d);
					postdata[a] = rnbinom_mu(eff_post_k, anmean2*post_rep_d);
				}
				for(a=0; a<pre_N; a++){
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					postsum += postdata[a];
				}				
			
			/* Suspend this	code: little is gained by avoiding a gamma draw, and power change from
				1->2 replicates is unpredictable because the distribution changes
				
			// If not using replicates then use combined post and efficacy k:
			}else if(pair_type == 3){
				for(a=0; a<max_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					anmean2 = eff_red_val * anmean1;
					temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					predata[a] = rpois(temp);
					temp = rgamma(eff_postcombined_k, (anmean2 / postcombined_k));  // replicates cancels out of scale parameter
					postdata[a] = rpois(temp);
				}
				for(a=0; a<pre_N; a++){
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					postsum += postdata[a];
				}				

			*/

			// pair_type should match one of the above but have this for safety:
			}else{
				for(i=0; i<iters[0]; i++){
					obsred[i] = -100;
				}
				break;				
			}

			// Rprintf("Pre sum: %i, Post sum: %i, red_val: %f\n", presum, postsum, red_val);
			
			// If the presum is 0 try to simulate data again:
			if(presum==0){
				continue;
			}
			
			// Save observed reduction:
			obsred[i] = 1.0 - (((double) postsum) / ((double) presum) * pre_post_ratio);
			
			// Estimate pre and post k

			estmu = (double) presum / pre_Nd;
			estvar = 0.0;
			for(a=0; a<pre_N; a++){
				estvar += pow((predata[a]-estmu), 2);
			}
			estvar = estvar / (pre_Nd - 1.0);
			// If the estimated variance is greater than the mean:
			if(estvar > estmu){
				est_prek = pow(estmu, 2) / (estvar - estmu);
				// If either variance < mean or close then use default of 10:
				if(est_prek > 10.0){
					est_prek = 10.0;
				}
			}else{
				est_prek = 10.0;
			}
			
			estmu = (double) postsum / post_Nd;
			estvar = 0.0;
			for(a=0; a<post_N; a++){
				estvar += pow((postdata[a]-estmu), 2);
			}
			estvar = estvar / (post_Nd - 1.0);
			// If the estimated variance is greater than the mean:
			if(estvar > estmu){
				est_postk = pow(estmu, 2) / (estvar - estmu);
				// If either variance < mean or close then use default of prek:
				if(est_postk > 10.0){
					est_postk = (est_prek / rep_pre[0]) * rep_post[0];
				}
			}else{
				est_postk = (est_prek / rep_pre[0]) * rep_post[0];
			}
			
			
			/***** Get BNB results *****/
//			ndfail += fecrt_bnb(presum, pre_N, est_prek, postsum, post_N, est_postk, lt, ut, conjugate_priors, delta_method, beta_it, &p1, &p2);
			// int fecrt_bnb(long long presum, int preN, double preK, double premean, double prevar, long long postsum, int postN, double postK, double postmean, 
			// double postvar, double H0_1, double H0_2, double *conjugate_priors, int delta, int beta_iters, double *p_1, double *p_2);
			
			// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
			if(p2 <= pt && p1 > pt){
				r = 0;
			}else if(p2 > pt && p1 > pt){
				r = 1;
			}else if(p2 <= pt && p1 <= pt){
				r = 2;
			}else if(p2 > pt && p1 <= pt){
				r = 3;
			}else{
				Rprintf("Unexpected failure for BNB method (with delta_method=%i) - values:  p1: %f, p2: %f  (presum: %i, preN: %i, prek: %f, postsum: %i, postN: %i, postk: %f, lt: %f, ut: %f)\n", delta_method, p1, p2, presum, pre_N, est_prek, postsum, post_N, est_postk, lt, ut);
				r = 4;
			}
			// Use indices 0-4 for BNB:
			classifications[r]++;
			

			/***** Get WAAVP results *****/
			if(postsum == 0){
				r = 4;
			}else{
				//waavp_ci(presum, predata, pre_N, postsum, postdata, post_N, pt, &w1, &w2);
				// void waavp_ci(long long presum, long *predata, int preN, long long postsum, long *postdata, int postN, double tail, double *ci_l, double *ci_u);

				// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
				if(w2 < ut && w1 <= lt){
					r = 0;
				}else if(w2 >= ut && w1 <= lt){
					r = 1;
				}else if(w2 < ut && w1 > lt){
					r = 2;
				}else if(w2 >= ut && w1 > lt){
					r = 3;
				}else{
					Rprintf("Unexpected failure for WAAVP method\n");
					r = 4;
				}
			}
			// Use indices 5-9 for WAAVP:
			classifications[r+5]++;
				
			// Note: not using Beta conjugate as too computationally expensive
			
			/***** Get Binomial conjugate results *****/
			if(postsum > presum){
				r = 4;
			}else{
				//dobson_ci(presum, postsum, lci_c[0], uci_c[0], dobson_priors, &w1, &w2);
				// void dobson_ci(long long presum, long long postsum, double lci_c, double uci_c, double *dobson_priors, double *ci_l, double *ci_u);
			
				// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
				if(w2 < ut && w1 <= lt){
					r = 0;
				}else if(w2 >= ut && w1 <= lt){
					r = 1;
				}else if(w2 < ut && w1 > lt){
					r = 2;
				}else if(w2 >= ut && w1 > lt){
					r = 3;
				}else{
					Rprintf("Unexpected failure for Dobson method with: presum=%i, postsum=%i\n", presum, postsum);
					r = 4;
				}
			}
			// Use indices 10-14 for Dobson:
			classifications[r+5+5]++;
				
			
			/***** Get Asymptotic method (paired and unpaired) results *****/
			if(pre_N != post_N){
				r = 4;
				// Use indices 15-19 for Asymptotic unpaired:
				classifications[r+5+5+5]++;
				// Use indices 20-24 for Asymptotic paired:
				classifications[r+5+5+5+5]++;
			}else{
				//asymptotic_ci(pre_N, presum, predata, postsum, postdata, poscov, pt, &w1, &w2, &w3, &w4);
				// void asymptotic_ci(int N, long long presum, long *predata, long long postsum, long *postdata, int poscov, double tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u);
				
				// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
				if(w2 < ut && w1 <= lt){
					r = 0;
				}else if(w2 >= ut && w1 <= lt){
					r = 1;
				}else if(w2 < ut && w1 > lt){
					r = 2;
				}else if(w2 >= ut && w1 > lt){
					r = 3;
				}else{
					Rprintf("Unexpected failure for Asymptotic unpaired method\n");
					r = 4;
				}
				// Use indices 15-19 for Asymptotic unpaired:
				classifications[r+5+5+5]++;
				
				// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure
				if(w4 < ut && w3 <= lt){
					r = 0;
				}else if(w4 >= ut && w3 <= lt){
					r = 1;
				}else if(w4 < ut && w3 > lt){
					r = 2;
				}else if(w4 >= ut && w3 > lt){
					r = 3;
				}else{
					Rprintf("Unexpected failure for Asymptotic paired method\n");
					r = 4;
				}
				// Use indices 20-24 for Asymptotic paired:
				classifications[r+5+5+5+5]++;
			}	
			
			i++;	
		}
		
		ndeltafail[0] = ndfail;
		PutRNGstate(); 
			
	}

}