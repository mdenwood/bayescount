/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for analysis and power calculations for FECRT studies
Last updated by MJD February 2017
*/

#include <math.h>
#include <R.h>
#include <Rmath.h>

// For std::nth_element used in ratio of betas MC approximation:
#include <vector>
#include <algorithm>
#include <functional>

#include "dists.h"

#include "power.h"

/**************************************************************/
/*   Beta negative binomial functions (not exported to R)     */
/**************************************************************/

double pbnb_lower(int q, double bnb_k, double bnb_alpha, double bnb_beta){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1;
	double p;

	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = phypergeometric(q, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = pgenhypergeometric(q, ghg_a, ghg_k, ghg_N, variety);
	
	return(p);
}
double pbnb_upper(int q, double bnb_k, double bnb_alpha, double bnb_beta){
	
	// Convert from beta-NB to generalised hypergeometric:
	double ghg_a = -bnb_beta;
	double ghg_k = -bnb_k;
	double ghg_N = bnb_alpha - 1;
	double p;
	
	// We redefine upper tail as inclusive, so if q=0:
	if(q==0){
		return(1.0);
	}
	
	hyperType variety = typeHyper(ghg_a, ghg_k, ghg_N);
	if (! checkHyperArgument(q, ghg_a, ghg_k, ghg_N, variety))
		p = NA_REAL;
	else if (variety==classic)
		p = 1.0 - phypergeometric(q-1, (int) ghg_a, (int) ghg_k, (int) ghg_N);
	else
		p = 1.0 - pgenhypergeometric(q-1, ghg_a, ghg_k, ghg_N, variety);
	// Note q-1 here to make p inclusive
	
	return(p);
}


/**************************************************************/
/*   Derivative and utility functions (not exported to R)     */
/**************************************************************/

// Simple function to convert mean and variance to alpha and beta:
void beta_params(double mu, double s2, double &alpha, double &beta){
	double ab = mu * (1.0 - mu) / s2 -1;
	alpha = mu * ab;
	beta = (1-mu)*ab;
}

// Non-linear transformation function:
double g_fun(double p, double r, double s, double t){
	return(p*r*t / (p*r*t - p*s + s));
}

// First 4 derivatives of g:
double gp1_fun(double p, double r, double s, double t){
	return(r*s*t / std::pow((p-1.0)*s - p*r*t, 2));
}
double gp2_fun(double p, double r, double s, double t){
	return(-2.0*r*s*t*(s - r*t) / std::pow((p-1.0)*s - p*r*t, 3));
}
double gp3_fun(double p, double r, double s, double t){
	return(6.0*r*s*t* std::pow(s - r*t, 2) / std::pow(p*r*t - p*s + s, 4));
}
double gp4_fun(double p, double r, double s, double t){
	return(24.0*r*s*t* std::pow(r*t - s, 3) / std::pow((p-1.0)*s - p*r*t, 5));
}

// Functions for expectations of powers 2-5 of a Beta distribution with parameters a(lpha), b(eta):
double epf2(double a, double b){
	return(((a+1.0)*a)  / ((a+b+1.0)*(a+b)));
}
double epf3(double a, double b){
	return(((a+2.0)*(a+1.0)*a)  / ((a+b+2.0)*(a+b+1.0)*(a+b)));
}
double epf4(double a, double b){
	return(((a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
}
double epf5(double a, double b){
	return(((a+4.0)*(a+3.0)*(a+2.0)*(a+1.0)*a)  / ((a+b+4.0)*(a+b+3.0)*(a+b+2.0)*(a+b+1.0)*(a+b)));
}

// The binomial expansions for E((x - mu)^t) for t in 2:5:
double expow2(double a, double b){
	double m = a / (a+b);
	return(epf2(a,b) - std::pow(m,2));
}
double expow3(double a, double b){
	double m = a / (a+b);
	return(epf3(a,b) - 3.0*m*epf2(a,b) + 2*std::pow(m,3));
}
double expow4(double a, double b){
	double m = a / (a+b);
	return(epf4(a,b) - 4.0*m*epf3(a,b) + 6.0*std::pow(m,2)*epf2(a,b) - 3.0*std::pow(m,4));
}
double expow5(double a, double b){
	double m = a / (a+b);
	return(epf5(a,b) - 5.0*m*epf4(a,b) + 10.0*std::pow(m,2)*epf3(a,b) - 10.0*std::pow(m,3)*epf2(a,b) + 4.0*std::pow(m,5));
}


// Delta method approximation to the change in mean (first 2 terms of Taylor series):
double delta_mean(double p_mu, double p_var, double r, double s, double t){
	
	// p_mu and p_var are calculated by calling function (for efficiency) as:
	// double p_mu = alpha / (alpha + beta);
	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));

	double rm = g_fun(p_mu, r, s, t) + 0.5 * gp2_fun(p_mu, r, s, t) * p_var;
	return(rm);
}

// Delta method approximation to the change in variance (higher order Taylor series):
double delta_var(double alpha, double beta, double p_mu, double p_var, double r, double s, double t){
	
	// p_mu and p_var are calculated by calling function (for efficiency) as:
	// double p_mu = alpha / (alpha + beta);
	// double p_var = (alpha * beta) / (std::pow(alpha + beta, 2) * (alpha + beta + 1));
	
	double g1 = gp1_fun(p_mu, r, s, t);
	double g2 = gp2_fun(p_mu, r, s, t);
	double g3 = gp3_fun(p_mu, r, s, t);
	double g4 = gp4_fun(p_mu, r, s, t);
	
	double rv = std::pow(g1,2) * p_var +
				2.0 * g1 * g2/2.0 * expow3(alpha, beta) +
				(std::pow(g2,2)/4.0 + 2.0*g1 * g3/6.0) * expow4(alpha, beta) +
				(2.0*g1 * g4/24.0 + g2 * g3/6.0) * expow5(alpha, beta);
	return(rv);
}


/*********************************************************************/
/*   Underlying C++ hypothesis test function (not exported to R)     */
/*********************************************************************/

void fecrt_pee(int presum, int preN, double preK, int postsum, int postN, double postK, double H0_1, double H0_2, double *prob_priors, int delta, int beta_iters, double *p_1, double *p_2){
	
	// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
	
	double alpha1 = presum + prob_priors[1];
	double beta1 = preN * preK + prob_priors[0];
	
	double tmu = alpha1 / (alpha1 + beta1);
	double tvar = (alpha1 * beta1) / (std::pow(alpha1 + beta1, 2) * (alpha1 + beta1 + 1));
	
	double effK = postK * (double) postN;
	
	/*  This is now done in the wrapper as no other C function calls this with delta=1	
	if(delta==0){
		GetRNGstate();
	}
	*/
	
	// Calculate the observed reduction
	// Note: this does not take account of replicates and edt but neither does H0_1 and H0_2 so all fine:
	double obsred = 1.0 - (((double) postsum) / ((double) presum));

	// First hypothesis:

	// Only calculate the pvalue when it makes sense to do so:
	if(obsred < H0_1){

		*p_1 = 1.0;

	}else{
		
		double meanchange = (1.0 - H0_1);
		if(meanchange < 0.0){
			meanchange = 0.0;
		}
		
		double newEprob = 0.0;
		double newVprob = 0.0;
		
		bool intensive=true;
		// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
		if(delta > 0){
			
			newEprob = delta_mean(tmu, tvar, preK, postK, meanchange);
			newVprob = delta_var(alpha1, beta1, tmu, tvar, preK, postK, meanchange);
			
			// If we get a negative or zero variance then fall back to not using the delta method:
			intensive = delta > 0 && (newEprob <= 0 || newVprob <= 0);
		}

		if(intensive){
			double rmean=0;
			double rvar=0;
			double delta=0;
			for(int i=1; i<=beta_iters; i++){
				delta = g_fun(rbeta(alpha1, beta1), preK, postK, meanchange) - rmean;
				rmean += delta / (double) i;
				rvar += delta*delta;
			}
			newEprob = rmean;
			newVprob = rvar / (double) (beta_iters-1);
		}
		
		double alpha2 = 0.0;
		double beta2 = 0.0;		
		beta_params(newEprob, newVprob, beta2, alpha2);
	
		// Note alpha and beta swapped as BNB is failures before successes and we want vice versa!
		*p_1 = pbnb_lower(postsum, effK, alpha2, beta2);
	
	}
	
	// Second hypothesis

	// Only calculate the pvalue when it makes sense to do so:
	if(obsred >= H0_2 || postsum == 0){

		*p_2 = 1.0;

	}else{
		
		double meanchange = (1.0 - H0_2);
		if(meanchange < 0.0){
			meanchange = 0.0;
		}
		
		double newEprob = 0.0;
		double newVprob = 0.0;
		
		bool intensive=true;
		// Note: delta can take 3 values:  0=never, 1=unless_fails, 2=always
		if(delta > 0){
			
			newEprob = delta_mean(tmu, tvar, preK, postK, meanchange);
			newVprob = delta_var(alpha1, beta1, tmu, tvar, preK, postK, meanchange);
			
			// If we get a negative or zero variance then fall back to not using the delta method:
			intensive = delta > 0 && (newEprob <= 0 || newVprob <= 0);
		}

		if(intensive){
			double rmean=0;
			double rvar=0;
			double delta=0;
			for(int i=1; i<=beta_iters; i++){
				delta = g_fun(rbeta(alpha1, beta1), preK, postK, meanchange) - rmean;
				rmean += delta / (double) i;
				rvar += delta*delta;
			}
			newEprob = rmean;
			newVprob = rvar / (double) (beta_iters-1);
		}
	
		double alpha2 = 0.0;
		double beta2 = 0.0;
		beta_params(newEprob, newVprob, alpha2, beta2);
		
		// Note alpha and beta swapped as BNB is failures before successes and we want vice versa!
		*p_2 = pbnb_upper(postsum, effK, beta2, alpha2);
		
	}
	
	/*  This is now done in the wrapper as no other C function calls this with delta=1	
	if(delta==0){
		PutRNGstate();
	}
	*/
	
}

/**********************************************************************/
/*   Underlying C++ function for WAAVP method (not exported to R)     */
/**********************************************************************/

void waavp_ci(int presum, int *predata, int preN, int postsum, int *postdata, int postN, double tail, double *ci_l, double *ci_u){
	
	if(postsum==0){

		*ci_l = 1.0;
		*ci_u = 1.0;

	}else{
		
		int df = preN + postN - 2;
		// Signature of qt is:  double  qt(double, double, int, int);
		double tval = qt(1.0 - tail, (double) df, 1, 0);
		
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
		
		double varred = varpre / (preN * premu * premu) + varpost / (postN * postmu * postmu);
		
		*ci_u = 1 - (postmu / premu * std::exp(-tval * std::sqrt(varred)));
		*ci_l = 1 - (postmu / premu * std::exp(tval * std::sqrt(varred)));
		
	}	
}

/***********************************************************************/
/*   Underlying C++ function for Dobson method (not exported to R)     */
/***********************************************************************/

void dobson_ci(int presum, int postsum, double lci_c, double uci_c, double *dobson_priors, double *ci_l, double *ci_u){
	
	*ci_l = 1.0 - qbeta(lci_c, postsum+dobson_priors[0], presum-postsum+dobson_priors[1], 0, 0);
	*ci_u = 1.0 - qbeta(uci_c, postsum+dobson_priors[0], presum-postsum+dobson_priors[1], 0, 0);
	
}

/*******************************************************************************/
/*   Underlying C++ function for conjugate Beta method (not exported to R)     */
/*******************************************************************************/

void conjbeta_ci(int preN, int presum, double preK, int postN, int postsum, double postK, int iters, double *prob_priors, double tail, double *ci_l, double *ci_u){
	
	std::vector<double> ratio(iters, 0);

//	Theoretically faster but makes no difference in practice:
//	std::vector<double> ratio;
//	ratio.reserve(iters);
	
	double alpha1 = preN * preK + prob_priors[0];
	double beta1 = presum + prob_priors[1];
	double alpha2 = postN * postK + prob_priors[0];
	double beta2 = postsum + prob_priors[1];
	
	for(int i=0; i<iters; i++){
		double p1 = rbeta(alpha1, beta1);
		double m1 = preK / p1 - preK;
		double p2 = rbeta(alpha2, beta2);
		double m2 = postK / p2 - postK;

		ratio[i] = 1.0 - m2/m1;
//		ratio.push_back(1.0 - m2/m1);
	}
	int in = (int) ((double) iters * tail);
	std::nth_element(ratio.begin(), ratio.begin() + in, ratio.end());
	*ci_l = ratio[in];

	in = (int) ((double) iters * tail);
	std::nth_element(ratio.begin(), ratio.begin() + in, ratio.end(), std::greater<double>());
	*ci_u = ratio[in];

}

/*******************************************************************************/
/*   Underlying C++ function for asymptotic methods (not exported to R)        */
/*******************************************************************************/

void asymptotic_ci(int N, int presum, int *predata, int postsum, int *postdata, double tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u){
	
	double Nd = (double) N;
	double premu = (double) presum / Nd;
	double postmu = (double) postsum / Nd;
	double r = postmu / premu;

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
	
	// Signature of qt is:  double  qt(double, double, int, int);
	// double critval = qt(0.975, Nd-1.0, 1, 0);
	// But it is MLE so normal is appropriate:
	// Signature of qnorm is:  qnorm(double, double, double, int, int)
	double critval = qnorm((double) 1.0 - tail, 0.0, 1.0, 1, 0);
	// double critval = 1.96;  
	
	double delta2u = std::pow(postmu, 2) / std::pow(premu, 4) * prevar + std::pow(premu, -2) * postvar;
	double delta2p = delta2u + 2 * postmu / std::pow(premu, 3) * covar;
	
	delta2u = std::sqrt(delta2u / Nd);
	delta2p = std::sqrt(delta2p / Nd);
	
	*u_ci_l = 1.0 - (r + critval * delta2u);
	*u_ci_u = 1.0 - (r - critval * delta2u);
	*p_ci_l = 1.0 - (r + critval * delta2p);
	*p_ci_u = 1.0 - (r - critval * delta2p);
	
}

