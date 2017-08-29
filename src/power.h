/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for older precision-type functions which are being phased out
Last updated by MJD February 2017
*/


double pbnb_lower(int q, double bnb_k, double bnb_alpha, double bnb_beta);

double pbnb_upper(int q, double bnb_k, double bnb_alpha, double bnb_beta);

void waavp_ci(int presum, int *predata, int preN, int postsum, int *postdata, int postN, double tail, double *ci_l, double *ci_u);

void dobson_ci(int presum, int postsum, double lci_c, double uci_c, double *dobson_priors, double *ci_l, double *ci_u);

void conjbeta_ci(int preN, int presum, double preK, int postN, int postsum, double postK, int iters, double *prob_priors, double tail, double *ci_l, double *ci_u);

void fecrt_pee(int presum, int preN, double preK, int postsum, int postN, double postK, double H0_1, double H0_2, double *prob_priors, int delta, int beta_iters, double *p_1, double *p_2);

void asymptotic_ci(int N, int presum, int *predata, int postsum, int *postdata, double tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u);
