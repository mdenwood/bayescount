/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for older precision-type functions which are being phased out
Last updated by MJD February 2017
*/


#include <R.h>

double myround(double number, int decimals);

// Precision for count data:
void precision_count(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarmcm, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *iterations, int *print, double *meancounts);

// Precision for reduction data:
void precision_reduction(double *premean, double *reduction, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, int *iterations, int *print, double *obsreds);