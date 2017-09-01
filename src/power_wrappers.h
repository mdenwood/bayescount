/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for older precision-type functions which are being phased out
Last updated by MJD February 2017
*/


void pbnb_lower_wrap(int *q, double *bnb_k, double *bnb_alpha, double *bnb_beta, double *p);

void pbnb_upper_wrap(int *q, double *bnb_k, double *bnb_alpha, double *bnb_beta, double *p);

void waavp_ci_wrap(int *presum, int *predata, int *preN, int *postsum, int *postdata, int *postN, double *tail, double *ci_l, double *ci_u);

void dobson_ci_wrap(int *presum, int *postsum, double *lci_c, double *uci_c, double *dobson_priors, double *ci_l, double *ci_u);

void conjbeta_ci_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, int *iters, double *prob_priors, double *tail, double *ci_l, double *ci_u);

void asymptotic_ci_wrap(int *N, int *presum, int *predata, int *postsum, int *postdata, double *tail, double *u_ci_l, double *u_ci_u, double *p_ci_l, double *p_ci_u);

void fecrt_pee_wrap(int *preN, int *presum, double *preK, int *postN, int *postsum, double *postK, double *H0_1, double *H0_2, double *prob_priors, int *delta, int *beta_iters, double *p_1, double *p_2);

void fecrt_pvals(int *iters, int *preN, int *postN, int *maxN, double *premean, double *reduction, double *edt_change, double *animalk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *prob_priors, double *kchange, double *truek, int *usetruek, int *delta, int *beta_iters, int *predata, int *postdata, double *p_1, double *p_2);

void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, double *tail, double *prob_priors, int *delta, int *beta_iters, int *predata, int *postdata, int *classifications, double *obsred);


