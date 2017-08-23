/*
This file is part of the bayescount package:  https://cran.r-project.org/package=bayescount
Copyright of Matthew James Denwood, licensed GPL>=2
The code is for older precision-type functions which are being phased out
Last updated by MJD February 2017
*/

#include "precision.h"

#include <math.h>
#include <Rmath.h>
#include <stdio.h>


double myround(double number, int decimals){
	number = number * pow(10,decimals);
	return ((double)((number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5))) / pow(10,decimals);
}


void precision_count(double *meanepg, double *gfaeces, double *sensitivity, int *replicates, int *animals, double *coeffvarmcm, double *coeffvarrep, double *coeffvarind, double *coeffvargroup, int *iterations, int *print, double *meancounts){

	double shapegp;
	double shapeindt[animals[0]];

	int i;

	for(i=0; i<animals[0]; i++){
		shapeindt[i] = pow(coeffvarmcm[i],2) + pow(coeffvarrep[i]/sqrt(gfaeces[i]),2) + (pow(coeffvarmcm[i],2)*pow(coeffvarrep[i]/sqrt(gfaeces[i]),2));
		shapeindt[i] = 1.0 / (shapeindt[i] + pow(coeffvarrep[i],2) + (shapeindt[i]*pow(coeffvarrep[i],2)));
	}

	shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

	//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

	double indmeans;
	double replicatemeans;

	int set, a, skipset;
	double lci, uci, meancount, sumcount;

	if(print[0]){
		Rprintf("< Running simulation >\n");
	}

	GetRNGstate();

	for(set=0; set<iterations[0]; set++){
		sumcount = 0.;
	
		for(a=0; a<animals[0]; a++){
			indmeans = rgamma(shapegp, (meanepg[0] / shapegp));
			replicatemeans = rgamma(shapeindt[a]*replicates[a], indmeans/(shapeindt[a]*replicates[a]));
			sumcount = sumcount + (((double)rpois(replicatemeans*replicates[a]*sensitivity[a]))*(1/sensitivity[a]))/replicates[a];
		}

		meancount = sumcount/animals[0];
	
		meancounts[set] = meancount;
	
		if(print[0]){
			Rprintf("\r%i%% complete", (int)(((double)set/(double)iterations[0])*100));
		}
	}

	PutRNGstate(); 

	if(print[0]){
		Rprintf("\r< Finished >                  \n");
	}

}

// Needs editing!
void precision_reduction(double *premean, double *reduction, int *replicates, int *animals, double *precoeffvarrep, double *precoeffvarind, double *precoeffvargroup, double *postcoeffvarrep, double *postcoeffvarind, double *postcoeffvargroup, int *iterations, int *print, double *obsreds){

	/*
	double shapegp;
	double shapeindt[animals[0]];

	int i;

	for(i=0; i<animals[0]; i++){
		shapeindt[i] = pow(coeffvarmcm[i],2) + pow(coeffvarrep[i]/sqrt(gfaeces[i]),2) + (pow(coeffvarmcm[i],2)*pow(coeffvarrep[i]/sqrt(gfaeces[i]),2));
		shapeindt[i] = 1.0 / (shapeindt[i] + pow(coeffvarrep[i],2) + (shapeindt[i]*pow(coeffvarrep[i],2)));
	}

	shapegp = 1/(coeffvargroup[0]*coeffvargroup[0]);

	//Rprintf("CVS:  %f, %f\n", shapeindt, shapegp);

	double indmeans;
	double replicatemeans;

	int set, a, skipset;
	double lci, uci, meancount, sumcount;

	if(print[0]){
		Rprintf("< Running simulation >\n");
	}

	GetRNGstate();

	for(set=0; set<maxiterations[0]; set++){
		sumcount = 0.;
	
		for(a=0; a<animals[0]; a++){
			indmeans = rgamma(shapegp, (meanepg[0] / shapegp));
			replicatemeans = rgamma(shapeindt[a]*replicates[a], indmeans/(shapeindt[a]*replicates[a]));
			sumcount = sumcount + (((double)rpois(replicatemeans*replicates[a]*sensitivity[a]))*(1/sensitivity[a]))/replicates[a];
		}

		meancount = sumcount/animals[0];
	
		meancounts[set] = meancount;
	
		if(print[0]){
			Rprintf("\r%i%% complete", (int)(((double)set/(double)maxiterations[0])*100));
		}
	}

	PutRNGstate(); 

	if(print[0]){
		Rprintf("\r< Finished >                  \n");
	}
	*/
}

