/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for analysis and power calculations for FECRT studies
Last updated by MJD February 2017
*/

/*****************************************/
/*   C wrappers to be exported to R      */
/*****************************************/

#include <math.h>
#include <R.h>
#include <Rmath.h>

#include "power.h"

/// TODO:  convert int presum and int postsum to long

extern "C" {
		
	#include "power_wrappers.h"
	
	/*****************************************/
	/*   Simple C wrappers                   */
	/*****************************************/
		
	void pbnb_lower_wrap(int *q, double *bnb_k, double *bnb_alpha, double *bnb_beta, double *p){
		p[0] = pbnb_lower(q[0], bnb_k[0], bnb_alpha[0], bnb_beta[0]);
	}
	void pbnb_upper_wrap(int *q, double *bnb_k, double *bnb_alpha, double *bnb_beta, double *p){
		p[0] = pbnb_upper(q[0], bnb_k[0], bnb_alpha[0], bnb_beta[0]);
	}
	
	void waavp_ci_wrap(int *presum, int *predata, int *preN, int *postsum, int *postdata, int *postN, double *tail, double *ci_l, double *ci_u){
		waavp_ci(presum[0], predata, preN[0], postsum[0], postdata, postN[0], tail[0], ci_l, ci_u);
	}
	void dobson_ci_wrap(int *presum, int *postsum, double *lci_c, double *uci_c, double *dobson_priors, double *ci_l, double *ci_u){
		dobson_ci(presum[0], postsum[0], lci_c[0], uci_c[0], dobson_priors, ci_l, ci_u);
	}
	void conjbeta_ci_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, int *iters, double *prob_priors, double *tail, double *ci_l, double *ci_u){
		conjbeta_ci(preN[0], presum[0], preK[0], postN[0], postsum[0], postK[0], iters[0], prob_priors, tail[0], ci_l, ci_u);
	}
	void asymptotic_ci_wrap(int *N, int *presum, int *predata, int *postsum, int *postdata, double *tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u){
		// It is NOT checked that preN==postN here - that must be done in R //
		asymptotic_ci(N[0], presum[0], predata, postsum[0], postdata, tail[0], u_ci_l, u_ci_u, p_ci_l, p_ci_u);
	}

	void fecrt_pee_wrap(int *presum, int *preN, double *preK, int *postsum, int *postN, double *postK, double *H0_1, double *H0_2, double *prob_priors, int *delta, int *beta_iters, double *p_1, double *p_2){
		if(delta[0]==0){
			GetRNGstate();
		}
		fecrt_pee(presum[0], preN[0], preK[0], postsum[0], postN[0], postK[0], H0_1[0], H0_2[0], prob_priors, delta[0], beta_iters[0], p_1, p_2);
		if(delta[0]==0){
			PutRNGstate();
		}
	}

	/*************************************************/
	/*   Power analysis wrapper (obtain p-values)    */
	/*************************************************/

	void fecrt_pvals(int *iters, int *preN, int *postN, int *maxN, double *premean, double *reduction, double *edt_change, double *animalk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *prob_priors, double *kchange, double *truek, int *usetruek, int *delta, int *beta_iters, int *predata, int *postdata, double *p_1, double *p_2){
		
		// Note:  predata and postdata must have length >= max(preN, postN)
		// Pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// reduction parameter should take account of edt_change
		
		// TODO:  BROKEN
		
		// TODO:  calculate pre and post k in function and clean up arguments
		// TODO:  vectorised variables (including 4 x k)
		
		GetRNGstate();
		
		int i, h, a, in;
		int presum, postsum;		
		double anmean, temp, k, mu, varsum;
		
		double simred = reduction[0] * edt_change[0];
		
		// First loop over iterations:
		for(i=0; i<iters[0]; i++){
			// Then simulate data using array allocated in R:
			for(a=0; a<maxN[0]; a++){
				if(animalk[0] >= 100.0){
					anmean = premean[0];
				}else{
					anmean = rgamma(animalk[0], (premean[0] / animalk[0]));
				}
				anmean = premean[0];
				temp = rgamma(prek[0], (anmean / prek[0]));
				predata[a] = rpois(temp);
				temp = rgamma(postk[0], ((anmean*simred) / postk[0]));
				postdata[a] = rpois(temp);
			}
			
			presum=0;
			for(int a=0; a<preN[0]; a++){
				presum += predata[a];
			}
			postsum=0;
			for(int a=0; a<postN[0]; a++){
				postsum += postdata[a];
			}
			if(usetruek[0] != 0){
				k = truek[0];
			}else{
				mu = (double) presum / (double) preN[0];
				varsum = 0;
				for(a=0; a<preN[0]; a++){
					varsum += pow((predata[a]-mu), 2);
				}
				k = pow(mu, 2) / ((varsum / (double) (preN[0] - 1)) - mu);
				// TODO:  pre and post k - and what happens when Inf or < 0??
			}
			
			// Then loop over H0 values:
			for(h=0; h<H0_N[0]; h++){
				// The array index to save to:
				in = i*H0_N[0] + h;
				// Pass pointers to correct index of p_1 and p_2:
				fecrt_pee(presum, preN[0], k, postsum, postN[0], k*kchange[0], H0_1[h], H0_2[h], prob_priors, delta[0], beta_iters[0], &p_1[in], &p_2[in]);
			}
		}
		
		PutRNGstate(); 
			
	}
	

	/*************************************************/
	/*   Power comparison for different tests        */
	/*************************************************/
	
	void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, double *tail, double *prob_priors, int *delta, int *beta_iters, int *predata, int *postdata, int *classifications, double *obsred){
		
		// Note:  predata and postdata must have length >= max(preN, postN)
		// Pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// If not paired model then rep_pre and rep_post must be 1
		// TODO:  Add pooling (and replicates directly here??)
		
		GetRNGstate();
		
		int i, a, in, r;
		long presum, postsum;
		double anmean1, anmean2, temp, est_prek, est_postk, estmu, estvar;
		double p1, p2, b1, b2, w1, w2, d1, d2;
		
		// Resolve pointers (and some casts) once:
		int delta_method = delta[0];  // Note: delta_method can take 3 values:  0=never, 1=unless_fails, 2=always but should probably be 1 or 2 here!
		int beta_it = beta_iters[0];  // Probably don't want too many iterations if the delta method fails!
		int pair_type = pairtype[0];
		int pre_N = preN[0];
		int post_N = postN[0];
		int max_N = maxN[0];
		double pre_Nd = (double) pre_N;
		double post_Nd = (double) post_N;
		double pt = tail[0];
		double red_val = 1.0 - reduction[0];
		double animal_k = animalk[0];
		double efficacy_k = efficacyk[0];
		double pre_k = prek[0];
		double post_k = postk[0];
		// Only needed for unpaired model with replicates:
		double treatment_k = (animal_k * efficacy_k) / (animal_k + efficacy_k + 1);
		// Only needed for paired model without replicates:
		double postcombined_k = (efficacy_k * post_k) / (efficacy_k + post_k + 1);
		double eff_postcombined_k = postcombined_k * (double) rep_post[0];
		
		// Take account of edt_change parameters:
		double eff_pre_mu = premean[0] / edt_pre[0];
		double eff_post_mu = (premean[0] * red_val) / edt_post[0];
		double eff_red_val = red_val * (edt_pre[0] / edt_post[0]);
		double eff_pre_k = pre_k * (double) rep_pre[0];
		double eff_post_k = post_k * (double) rep_post[0];
		
		// This correction is only used for the hypothesis test thresholds and calculating the obs reduction:
		double correction = ((double) rep_post[0] * (double) edt_pre[0]) / ((double) rep_pre[0] * (double) edt_post[0]);
		double lt = 1.0 - ((1.0 - H0_1[0]) * correction);
		double ut = 1.0 - ((1.0 - H0_2[0]) * correction);
		double pre_post_ratio = (pre_Nd / post_Nd) * (1.0 / correction);
				
		// First loop over iterations:
		for(i=0; i<iters[0]; i++){
			// Then simulate data using array allocated in R:
			
			presum = 0;
			postsum = 0;
			
			// Simple unpaired model ignores animalk and efficacyk:			
			if(pair_type == 0){
				// We only need pre and post data where required:
				for(a=0; a<pre_N; a++){
					temp = rgamma(pre_k, (eff_pre_mu / pre_k));
					predata[a] = rpois(temp);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					temp = rgamma(post_k, (eff_post_mu / post_k));
					postdata[a] = rpois(temp);
					postsum += postdata[a];
				}

			// Unpaired model but with replicates needs to use animalk and efficacyk:
			}else if(pair_type == 1){
				for(a=0; a<pre_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					predata[a] = rpois(temp);
					presum += predata[a];
				}
				for(a=0; a<post_N; a++){
					anmean2 = rgamma(treatment_k, (eff_post_mu / treatment_k));
					temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					postdata[a] = rpois(temp);
					postsum += postdata[a];
				}
			
			// If it is a paired model we need pre and post means for all animals:
			// (don't need all poisson variates but do them to avoid another loop
			// or conditional, as preN nearly always = postN for paired anyway):
			}else if(pair_type == 2){
				for(a=0; a<max_N; a++){
					anmean1 = rgamma(animal_k, (eff_pre_mu / animal_k));
					anmean2 = anmean1 * rgamma(efficacy_k, (eff_red_val / efficacy_k));
					temp = rgamma(eff_pre_k, (anmean1 / pre_k));  // replicates cancels out of scale parameter
					predata[a] = rpois(temp);
					temp = rgamma(eff_post_k, (anmean2 / post_k));  // replicates cancels out of scale parameter
					postdata[a] = rpois(temp);
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
			

			// Save observed reduction:
			obsred[i] = 1.0 - (((double) postsum) / ((double) presum) * pre_post_ratio);
			
			// Estimate pre and post k
			estmu = (double) presum / pre_Nd;
			
			estvar = 0.0;
			for(a=0; a<pre_N; a++){
				estvar += pow((predata[a]-estmu), 2);
			}
			estvar = estvar / (pre_Nd - 1.0);
			
			est_prek = pow(estmu, 2) / (estvar - estmu);

			/* 	Not sure if I should do this or base prek on postk in this situation??
				Or simply generate another dataset??
			*/
			// If we have data with less variation than Poisson (or close to), make prek quite big:
			if(estvar <= estmu || est_prek > 10.0){
				est_prek = 10.0;
			}
				
			estmu = (double) postsum / post_Nd;
			estvar = 0.0;
			for(a=0; a<post_N; a++){
				estvar += pow((postdata[a]-estmu), 2);
			}
			estvar = estvar / (post_Nd - 1.0);
			est_postk = pow(estmu, 2) / (estvar - estmu);
			
			// If postk can't be calculated then base on prek:
			if(postsum==0 || estvar <= estmu || est_postk > 10.0){
				est_postk = (est_prek / rep_pre[0]) * rep_post[0];
			}
			
			// Get pee results:
			fecrt_pee(presum, pre_N, est_prek, postsum, post_N, est_postk, lt, ut, prob_priors, delta_method, beta_it, &p1, &p2);
			
			// Results are 1: reduced, 2: inconclusive, 3: marginal, 4: adequate, 5: method failure - but indexing starts at 1
			if(est_prek < 0 || est_postk < 0){
				r = 4;
			}else if(p2 <= pt && p1 > pt){
				r = 0;
			}else if(p2 > pt && p1 > pt){
				r = 1;
			}else if(p2 <= pt && p1 <= pt){
				r = 2;
			}else if(p2 > pt && p1 <= pt){
				r = 3;
			}else{
				Rprintf("Failure values (with delta_method=%i):  p1: %f, p2: %f  (presum: %i, preN: %i, prek: %f, postsum: %i, postN: %i, postk: %f, lt: %f, ut: %f)\n", delta_method, p1, p2, presum, pre_N, est_prek, postsum, post_N, est_postk, lt, ut);
				r = 4;
			}
			classifications[r]++;
			
			// Get WAAVP results:
			waavp_ci(presum, predata, pre_N, postsum, postdata, post_N, pt, &w1, &w2);
				
			// Results are 0: method failure, 1: reduced, 2: inconclusive, 3: marginal, 4: adequate
			if(postsum == 0){
				r = 4;
			}else if(w2 < ut && w1 <= lt){
				r = 0;
			}else if(w2 >= ut && w1 <= lt){
				r = 1;
			}else if(w2 < ut && w1 > lt){
				r = 2;
			}else if(w2 >= ut && w1 > lt){
				r = 3;
			}else{
				r = 4;
			}
			// Use indices 5-9 for WAAVP:
			classifications[r+5]++;
				
			
			/*
			if(pre_N==post_N){
				
			TODO:  Add asymptotoc results
				asymptotic_ci(preN[0], int presum, int *predata, int postsum, int *postdata, double tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u);
			}
			
			// Then classify results according to 3 categories for each (3*3*3*3) and observed reduction = 1 / < 1 (*2) - so 162 categories (index 0-161)
			
			in = 0;
			
			// If reduction > 0 start at 81st element:
			if(postsum > 0){
				in = 81;
			}
			
			r = 1;
			if(p1 <= pt){
				r = 0;
			}else if(p2 <= pt){
				r = 2;
			}
			in += r*3*3*3;

			r = 1;
			if(b1 >= lt){
				r = 0;
			}else if(b2 <= ut){
				r = 2;
			}
			in += r*3*3;
			
			r = 1;
			if(d1 >= lt){
				r = 0;
			}else if(d2 <= ut){
				r = 2;
			}
			in += r*3;
			
			r = 1;
			if(w1 >= lt){
				r = 0;
			}else if(w2 <= ut){
				r = 2;
			}
			in += r;
			
			classifications[in]++;
			
			*/
				
		}
		
		PutRNGstate(); 
			
	}

}