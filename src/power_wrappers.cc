/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for analysis and power calculations for FECRT studies
Last updated by MJD February 2017
*/

/*****************************************/
/*   C wrappers to be exported to R)     */
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
	
	void waavp_ci_wrap(int *presum, int *predata, int *preN, int *postsum, int *postdata, int *postN, double *ci_l, double *ci_u){
		waavp_ci(presum[0], predata, preN[0], postsum[0], postdata, postN[0], ci_l, ci_u);
	}
	void dobson_ci_wrap(int *presum, int *postsum, double *lci_c, double *uci_c, double *dobson_priors, double *ci_l, double *ci_u){
		dobson_ci(presum[0], postsum[0], lci_c[0], uci_c[0], dobson_priors, ci_l, ci_u);
	}
	void conjbeta_ci_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, int *iters, double *prob_priors, double *lci_b, double *uci_b, double *ci_l, double *ci_u){
		conjbeta_ci(preN[0], presum[0], preK[0], postN[0], postsum[0], postK[0], iters[0], prob_priors, lci_b[0], uci_b[0], ci_l, ci_u);
	}
	void asymptotic_ci_wrap(int *N, int *presum, int *predata, int *postsum, int *postdata, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u){
		// It is NOT checked that preN==postN here - that must be done in R //
		asymptotic_ci(N[0], presum[0], predata, postsum[0], postdata, u_ci_l, u_ci_u, p_ci_l, p_ci_u);
	}

		
	void fecrt_pee_wrap(int *predata, int *preN, int *presum, double *preK, int *postdata, int *postN, int *postsum, double *postK, double *H0_1, double *H0_2, double *edt_change, double *prob_priors, int *delta, int *beta_iters, double *p_1, double *p_2){
		fecrt_pee(presum[0], predata, preN[0], preK[0], postsum[0], postdata, postN[0], postK[0], H0_1[0], H0_2[0], edt_change[0], prob_priors, delta[0], beta_iters[0], p_1, p_2);
	}

	/*************************************************/
	/*   Power analysis wrapper (obtain p-values)    */
	/*************************************************/

	void fecrt_pvals(int *iters, int *preN, int *postN, int *maxN, double *premean, double *reduction, double *edt_change, double *animalk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *prob_priors, double *kchange, double *truek, int *usetruek, int *delta, int *beta_iters, int *predata, int *postdata, double *p_1, double *p_2){
		
		// Note:  predata and postdata must have length >= max(preN, postN)
		// Pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// reduction parameter should take account of edt_change
		
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
				fecrt_pee(presum, predata, preN[0], k, postsum, postdata, postN[0], k*kchange[0], H0_1[h], H0_2[h], edt_change[0], prob_priors, delta[0], beta_iters[0], &p_1[in], &p_2[in]);
			}
		}
		
		PutRNGstate(); 
			
	}
	
	/*************************************************/
	/*   Power comparison for different tests        */
	/*************************************************/
	
	void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, double *premean, double *reduction, double *edt_change, double *animalk, double *prek, double *postk, double *H0_1, double *H0_2, double *lci_c, double *uci_c, double *lci_b, double *uci_b, double *dobson_priors, double *tval, double *psig, double *prob_priors, double *kchange, double *truek, int *usetruek, int *delta, int *beta_iters, int *predata, int *postdata, int *classifications){
		
		// Note:  predata and postdata must have length >= max(preN, postN)
		// Pre and post data will be simulated for max(preN, postN) but only used for preN / postN respectively
		// reduction parameter should take account of edt_change
		
		// TODO:  remove tval
		// TODO:  same pre-data can be used for each reduction
		
		GetRNGstate();
		
		int i, a, presum, postsum, in, r;
		double anmean, temp, k, mu, varsum;
		double p1, p2, b1, b2, w1, w2, d1, d2;
		
		double lt = H0_1[0];
		double ut = H0_2[0];
		double pt = psig[0];
		
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
				temp = rgamma(prek[0], (anmean / prek[0]));
				predata[a] = rpois(temp);
				temp = rgamma(postk[0], ((anmean*simred) / postk[0]));
				postdata[a] = rpois(temp);
			}
			
			presum = 0;
			for(a=0; a<preN[0]; a++){
				presum += predata[a];
			}
			postsum = 0;
			for(a=0; a<postN[0]; a++){
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
				k = pow(mu, 2) / (varsum / (double) (preN[0] - 1));
			}
			
			// Get results:
			waavp_ci(presum, predata, preN[0], postsum, postdata, postN[0], &w1, &w2);
			dobson_ci(presum, postsum, lci_c[0], uci_c[0], dobson_priors, &d1, &d2);
			conjbeta_ci(preN[0], presum, k, postN[0], postsum, k*kchange[0], beta_iters[0], prob_priors, lci_b[0], uci_b[0], &b1, &b2);
			fecrt_pee(presum, predata, preN[0], k, postsum, postdata, postN[0], k*kchange[0], H0_1[0], H0_2[0], edt_change[0], prob_priors, delta[0], beta_iters[0], &p1, &p2);
			
			// Then classify results according to 3 categories for each (3*3*3*3) and observed reduction = 1 / < 1 (*2) - so 162 categories (index 0-161)
			
			in = 0;
			
			// If reduction > 0 start at 81st element:
			if(postsum > 0){
				in = 81;
			}
			
			// Results are 0: susceptible, 1: inconclusive, 2: resistant
			
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
				
		}
		
		PutRNGstate(); 
			
	}

}