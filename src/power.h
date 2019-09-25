/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for power calculations for FECRT studies
Last updated by MJD June 2018
*/


void fecrt_power(int *iters, int *paired_analysis, int *use_truek, double *kchange, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, int *H0_N, double *prob_priors, int *delta, int *beta_iters, int *approximation, double *p_1, double *p_2, double *obsred, int *ndeltafail);

void fecrt_power_comparison(int *iters, int *preN, int *postN, int *maxN, int *rep_pre, int *rep_post, double *edt_pre, double *edt_post, double *premean, double *reduction, int *pairtype, double *animalk, double *efficacyk, double *prek, double *postk, double *H0_1, double *H0_2, double *tail, double *prob_priors, int *delta, int *beta_iters, int poscov, double *lci_c, double *uci_c, double *dobson_priors, int *classifications, double *obsred, int *ndeltafail);

/*
#### TODO
# Modify and check fecrt_power_comparison
# Website interface and fecrt_power etc could have paired type with 2 k values - ani+pre then pre+efficacy+post - use 2 new k parameters to avoid confusion and set ones not using to negative to be sure
# Modify and check fecrt_power
# Remove functions below here or leave as shims for the function above
# Analysis method could quickly run sims and advise which has optimal rates based on data characteristics from simplified optim (unless 100% reduction in which case guess true reduction and post k)
 */