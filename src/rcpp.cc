#include "Rcpp.h"

#include "fecrt.h"


// [[Rcpp::export]]
double pbnb(int q, double k, double alpha, double beta, bool lower, bool inclusive){
	
	double p = NA_REAL;
	
	if(lower){
		if(!inclusive){
			if(q==0L){
				return 0.0;
			}
			q--;
		}
		p = pbnb_lower(q, k, alpha, beta);
	}else{
		if(inclusive){
			if(q==0L){
				return 1.0;
			}
			q--;
		}
		p = 1.0 - pbnb_lower(q, k, alpha, beta);
	}
	
	return p;
}


// [[Rcpp::export]]
Rcpp::NumericVector estimate_k_useless(const double mean_1, const double var_1, const double mean_2, const double var_2, const double cov_12, const bool paired){
	
	Rcpp::NumericVector retval(2);
	
	// Remove Poisson variation:
	double var_eff_1 = var_1 - mean_1;
	double var_eff_2 = var_2 - mean_2;
	
	double k1 = (mean_1 * mean_1) / var_eff_1;
	double k2 = (mean_2 * mean_2) / var_eff_2;
	
	// Adjust for the covariance if necessary and possible (cov_12 should be ignored if not positive):
	if(paired && var_eff_1 > 0.0 && var_eff_2 > 0.0 && cov_12 > 0.0){
		
		// Mean-adjusted covariance (Poisson variance cancels out at this stage):
		double adj_cov = cov_12 / (mean_1 * mean_2);
	
		// Adjust to the uncorrelated k components:
		k1 = (adj_cov*mean_1*mean_1 + k1*var_eff_1) / (var_eff_1 - adj_cov*mean_1*mean_1);
		k2 = (adj_cov*mean_2*mean_2 + k2*var_eff_2) / (var_eff_2 - adj_cov*mean_2*mean_2);
		
		// Or should k1 be the correlated component?
		
		// If this generates a negative value then revert to the uncorrelated k estimates:
		if(k1 <= 0.0){
			k1 = (mean_1 * mean_1) / var_eff_1;
		}
		if(k2 <= 0.0){
			k2 = (mean_2 * mean_2) / var_eff_2;
		}
	}
	
	// TODO: this assumes that mean_2 is lower than mean_1 (so k1 is better estimated than k2) but that may not be true - detect and switch for below here
	
	// Can happen if the control variance is less than the control mean,
	// in which case re-calculate ignoring Poisson and correlation:
	if( var_eff_1 <= 0.0 ){
		retval[0] = (mean_1 * mean_1) / var_1;
	}else{
		retval[0] = k1;
	}

	// Can happen if the treatment variance is either 0 or less than the control mean,
	// in which case re-calculate ignoring Poisson and correlation:
	if( var_eff_2 <= 0.0 ){
		if(var_2 > 0.0){
			retval[1] = (mean_2 * mean_2) / var_2;
		}else{
			// If there is no variance at all assume k2=k1:
			retval[1] = retval[0];
		}
		// TODO: more sophisticated heuristic of when to use k1 and when to use observed data_2 when var_eff_2 < 0??
		retval[1] = retval[0];
	}else{
		retval[1] = k2;
	}
	
	// Test:
	if(retval[0] <= 0.0){
		Rcpp::stop("k1 below 0.0 generated");
	}
	if(retval[1] <= 0.0){
		Rcpp::stop("k2 below 0.0 generated");
	}
	
	return retval;
	
}

// [[Rcpp::export]]
Rcpp::NumericVector estimate_k_experiment(const double mean_1, const double var_1, const double mean_2, const double var_2, const double cov_12, const bool paired){
	
	Rcpp::NumericVector retval(2);
	
	// Remove Poisson variation:
	double var_eff_1 = var_1 - mean_1;
	double var_eff_2 = var_2 - mean_2;
	
	double k1 = (mean_1 * mean_1) / var_eff_1;
	double k2 = (mean_2 * mean_2) / var_eff_2;
	
	// Adjust for the covariance if necessary and possible (cov_12 should be ignored if not positive):
	if(paired){
		
		// Mean-adjusted covariance (Poisson variance cancels out at this stage):
		double adj_cov = cov_12 / (mean_1 * mean_2);
	
		// Adjust to the uncorrelated k components:
		double uk1 = (adj_cov*mean_1*mean_1 + k1*var_eff_1) / (var_eff_1 - adj_cov*mean_1*mean_1);
		double uk2 = (adj_cov*mean_2*mean_2 + k2*var_eff_2) / (var_eff_2 - adj_cov*mean_2*mean_2);
		
		double kc = 1 / adj_cov;
		
		k1 = k1 / (1- std::sqrt(k1)*std::sqrt(k2)/kc);
		k2 = k2 / (1- std::sqrt(k1)*std::sqrt(k2)/kc);
		
	}
	
	// TODO: this assumes that mean_2 is lower than mean_1 (so k1 is better estimated than k2) but that may not be true - detect and switch for below here
	
	// Can happen if the control variance is less than the control mean,
	// in which case re-calculate ignoring Poisson and correlation:
	if( k1 <= 0.0 ){
		retval[0] = (mean_1 * mean_1) / var_1;
	}else{
		retval[0] = k1;
	}

	// Can happen if the treatment variance is either 0 or less than the control mean,
	// in which case re-calculate ignoring Poisson and correlation:
	if( k2 <= 0.0 ){
		if(var_2 > 0.0){
			retval[1] = (mean_2 * mean_2) / var_2;
		}else{
			// If there is no variance at all assume k2=k1:
			retval[1] = retval[0];
		}
		// TODO: more sophisticated heuristic of when to use k1 and when to use observed data_2 when var_eff_2 < 0??
		retval[1] = retval[0];
	}else{
		retval[1] = k2;
	}
	
	// Test:
	if(retval[0] <= 0.0){
		Rcpp::stop("k1 below 0.0 generated");
	}
	if(retval[1] <= 0.0){
		Rcpp::stop("k2 below 0.0 generated");
	}
	
	return retval;
	
}


// [[Rcpp::export]]
Rcpp::NumericVector estimate_k(const double mean_1, const double var_1, const double mean_2, const double var_2, const double cov_12, const bool paired){
	
	Rcpp::NumericVector retval(2);
	
	// Remove Poisson variation:
	double var_eff_1 = var_1 - mean_1;
	double var_eff_2 = var_2 - mean_2;
	
	// If using the paired model and both adjusted variances are positive:
	double correlation = 0.0;
	if(paired && var_eff_1 > 0.0 && var_eff_2 > 0.0){	
		// Assume that none of the covariance is due to Poisson variation:
		// correlation = cov_12 / (std::sqrt(var_eff_1) * std::sqrt(var_eff_2));
	}
	// Empirically this seems to give the best results:
	if(paired && var_1 > 0.0 && var_2 > 0.0){	
		// Assume that none of the covariance is due to Poisson variation:
		correlation = cov_12 / (std::sqrt(var_1) * std::sqrt(var_2));
	}
	// TODO: correlation can come out negative here - make min 0.0??
	// TODO: ks are over-estimated when covariance is huge
	
	if(correlation < 0.0){
		//Rcpp::warning("Negative correlation generated");
	}
	if(correlation >= 1.0){
		Rcpp::warning("correlation >= 1.0 corrected");
		correlation = 0.99;
	}
	
	
	// Can happen if the control variance is less than the control mean:
	if( var_eff_1 <= 0.0 ){
		retval[0] = (mean_1 * mean_1) / (var_1 * (1.0 - correlation));
	}else{
		retval[0] = (mean_1 * mean_1) / (var_eff_1 * (1.0 - correlation));
	}
	
	// Can happen if the treatment variance is either 0 or less than the control mean:
	if( var_eff_2 <= 0.0 ){
		retval[1] = retval[0];
	}else{
		retval[1] = (mean_2 * mean_2) / (var_eff_2 * (1.0 - correlation));
	}
	
	return retval;
	
}


// Forward declaration:
double find_theta(const Rcpp::IntegerVector data, const double mu, const double ax, const double bx, const double tol);

// [[Rcpp::export]]
Rcpp::NumericVector estimate_k_ml(const Rcpp::IntegerVector data_1, const double mean_1, const double var_1, const Rcpp::IntegerVector data_2, const double mean_2, const double var_2, const double cov_12, const bool paired){
	
	Rcpp::NumericVector retval(2);
	
	// Get MLE estimates of two k values:
	double k1 = find_theta(data_1, mean_1, 0.001, 20.0, 0.01);
	double k2 = find_theta(data_2, mean_2, 0.001, 20.0, 0.01);
	
	if(paired && var_1 > 0.0 && var_2 > 0.0 && cov_12 > 0.0){
		double correlation = 0.0;

		// Assuming that none of the covariance is due to Poisson variation:
		// correlation = cov_12 / (std::sqrt(var_1 - mean_1) * std::sqrt(var_2 - mean_2));

		// Empirically this seems to give the best results:
		correlation = cov_12 / (std::sqrt(var_1) * std::sqrt(var_2));
		
		if(correlation >= 1.0){
			Rcpp::warning("correlation >= 1.0");
			correlation = 0.99;
		}
		
		// Adjust k:
		k1 /= (1.0 - correlation);
		k2 /= (1.0 - correlation);		
	}
	
	// TODO: might be that var_2 > 0 but var_1 = 0??
	if(var_2 <= 0.0){
		k2 = k1;
	}
	
	if(!R_finite(k2) || k2 <= 0.0){
		Rcpp::Rcout << k2 << " - " << mean_2 << " - " << var_2 << " - " << cov_12 << "\n";
	}
	
	retval[0] = k1;
	retval[1] = k2;
	
	return retval;
	
}



// [[Rcpp::export]]
Rcpp::String get_type_ci(const double eff, const double lci, const double uci, const double ta, const double ti){

	Rcpp::String tp;
	
	if(uci < ta){
		tp = "1a";
	}else if(uci < ti && eff < ta){
		tp = "1b";
	}else if(uci < ti && lci < ta){
		tp = "1c";
	}else if(lci < ta && uci >= ti && eff < ta){
		tp = "2a";
	}else if(lci < ta && uci >= ti && eff >= ta && eff < ti){
		tp = "2b";
	}else if(lci < ta && uci >= ti && eff >= ti){
		tp = "2c";
	}else if(lci >= ta && uci < ti){
		tp = "3";
	}else if(lci >= ta && uci >= ti && eff < ti){
		tp = "4a";
	}else if(lci >= ta && uci >= ti && eff >= ti){
		tp = "4b";
	}else if(lci >= ti){
		tp = "4c";
	}else{
//		Rcpp::Rcout << eff << " - " << lci << " - " << uci << " - " << ta << " - " << ti << "\n";
		tp = "error";
//		Rcpp::stop("Error assigning ci type");
	}

	return tp;
}

// [[Rcpp::export]]
Rcpp::String get_type_pv(const double eff, const double p1, const double p2, const double tail, const double ta, const double ti){
	
	Rcpp::String tp;
	
	if(p2 <= tail && p1 > tail && eff < ta){
		tp = "1ab";
	}else if(p2 <= tail && p1 > tail && eff >= ta){
		tp = "1c";
	}else if(p2 > tail && p1 > tail && eff < ta){
		tp = "2a";
	}else if(p2 > tail && p1 > tail && eff >= ta && eff < ti){
		tp = "2b";
	}else if(p2 > tail && p1 > tail && eff >= ti){
		tp = "2c";
	}else if(p2 <= tail && p1 <= tail){
		tp = "3";
	}else if(p2 > tail && p1 <= tail && eff < ti){
		tp = "4a";
	}else if(p2 > tail && p1 <= tail && eff >= ti){
		tp = "4bc";
	}else{
		tp = "error";
//		Rcpp::stop("Error assigning pv type");
	}

	return tp;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix methodcomp(long long sum1, int N1, double K1_p, double K1_u, double mu1, double var1, long long sum2, int N2, double K2_p, double K2_u, double mu2, double var2, double cov12, double mean_ratio, double H0_1, double H0_2, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, Rcpp::NumericVector dobson_cl, Rcpp::NumericVector dobson_priors)
{
	
	Rcpp::NumericMatrix output(9, 2);
	
	// TODO: add Dobson method and int switch for paired, unpaired, both and calculate sum/N/ks here from data_1 and data_2
	
	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };

	int i=0;

	bnb_pval(sum1, N1, K1_p, mu1, var1, sum2, N2, K2_p, mu2, var2, cov12, mean_ratio, H0_1, H0_2, conjugate_priors_db, delta, beta_iters, approx, &output(i,0), &output(i,1)); i++;
	bnb_pval(sum1, N1, K1_u, mu1, var1, sum2, N2, K2_u, mu2, var2, 0.0, mean_ratio, H0_1, H0_2, conjugate_priors_db, delta, beta_iters, approx, &output(i,0), &output(i,1)); i++;
	
	levecke_p_ci(mu1, mu2, var1, var2, cov12, N1, tail, &output(i,0), &output(i,1)); i++;
	levecke_u_ci(mu1, mu2, var1, var2, N1, N2, tail, &output(i,0), &output(i,1)); i++;
	
	waavp_p_ci(mu1, mu2, var1, var2, cov12, N1, tail, &output(i,0), &output(i,1)); i++;
	waavp_u_ci(mu1, mu2, var1, var2, N1, N2, tail, &output(i,0), &output(i,1)); i++;
	
	mle_p_ci(mu1, mu2, var1, var2, cov12, N1, tail, &output(i,0), &output(i,1)); i++;
	mle_u_ci(mu1, mu2, var1, var2, N1, N2, tail, &output(i,0), &output(i,1)); i++;
	
	dobson_ci(sum1, sum2, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output(i,0), &output(i,1)); i++;

	return output;
}


// [[Rcpp::export]]
Rcpp::DataFrame fecrt_sim_unpaired(int iters, Rcpp::NumericVector red, int N_tx, int N_ctl, double mu, double k_tx, double k_ctl, Rcpp::NumericMatrix thresholds, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, bool useml)
{

	// Check that the thresholds have 2 columns and >0 rows:
	if(thresholds.ncol() != 2){
		Rcpp::stop("There must be exactly 2 columns");
	}
	if(thresholds.nrow() < 1){
		Rcpp::stop("There must be 1 or more rows");
	}
	for(int i=0; i<thresholds.nrow(); i++){
		if(thresholds(i,1) < thresholds(i,0)){
			Rcpp::stop("Second thresholds must be higher than the first");
		}
		if(thresholds(i,0) < 0 || thresholds(i,0) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
		if(thresholds(i,1) < 0 || thresholds(i,1) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
	}
	
	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");
	
	Rcpp::StringVector methods = { "BNB", "Levecke", "WAAVP", "MLE" };
	Rcpp::IntegerVector tseq = Rcpp::seq(1L, thresholds.nrow());
	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);
	
	Rcpp::DataFrame output = expGrid(methods, red, tseq, td, td, iteration, td, td, td, ts, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Method", "Reduction", "ThresholdIndex", "ThresholdLower", "ThresholdUpper", "Iteration", "ObsReduction", "Stat1", "Stat2", "Classification");

	// Get references to list items:
	Rcpp::StringVector output_Method = output[0L];
	Rcpp::DoubleVector output_Reduction = output[1L];
	Rcpp::IntegerVector output_ThresholdIndex = output[2L];
	Rcpp::DoubleVector output_ThresholdLower = output[3L];
	Rcpp::DoubleVector output_ThresholdUpper = output[4L];
	Rcpp::IntegerVector output_Iteration = output[5L];
	Rcpp::DoubleVector output_ObsReduction = output[6L];
	Rcpp::DoubleVector output_Stat1 = output[7L];
	Rcpp::DoubleVector output_Stat2 = output[8L];
	Rcpp::StringVector output_Classification = output[9L];
	
	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	
	// To control the row index to write to:
	int row = 0L;
	
	// TODO: tidy up check and error:
	if((double(N_ctl)*double(mu)) > (0.1 * INT_MAX) || (double(N_ctl)*double(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}
	
	// First loop over control datasets to simulate:
	for(int i=0; i<iters; i++){		
		Rcpp::IntegerVector ctl_data;
		long long sum_ctl = 0L;
		do {
			sum_ctl = 0L;
			ctl_data = Rcpp::rnbinom_mu(N_ctl, k_ctl, mu);
			for(int i=0; i<N_tx; i++){
				sum_ctl += ctl_data[i];
			}
		} while (sum_ctl == 0L);

		double mean_ctl = Rcpp::mean(ctl_data);
		double var_ctl = Rcpp::var(ctl_data);
		
		// Then loop over reductions:
		for(int r=0; r<red.length(); r++){
			Rcpp::IntegerVector tx_data;
			tx_data = Rcpp::rnbinom_mu(N_tx, k_tx, mu*red[r]);
			
			// Avoid possible segfaults with int vs long long:
			long long sum_tx = 0L;
			for(int i=0; i<N_tx; i++){
				sum_tx += tx_data[i];
			}
			double mean_tx = Rcpp::mean(tx_data);		
			double var_tx = Rcpp::var(tx_data);
			double obsred = mean_tx / mean_ctl;

			Rcpp::NumericVector ks;
			if(useml){
				ks = estimate_k_ml(ctl_data, mean_ctl, var_ctl, tx_data, mean_tx, var_tx, 0.0, false);
			}else{
				ks = estimate_k(mean_ctl, var_ctl, mean_tx, var_tx, 0.0, false);
			}
			// Test:
			if(!R_finite(ks[0])){
				Rcpp::stop("Non-finite k1 generated");
			}
			if(ks[0] <= 0.0){
				Rcpp::stop("k1 below 0.0 generated");
			}
			if(!R_finite(ks[1])){
				Rcpp::stop("Non-finite k2 generated");
			}
			if(ks[1] <= 0.0){
				Rcpp::stop("k2 below 0.0 generated");
			}
			
			double estk_ctl = ks[0];
			double estk_tx = ks[1];
			
			// For cheating:
			//estk_ctl = k_ctl;
			//estk_tx = k_tx;
			
			// Then loop over thresholds:
			for(int t=0; t<thresholds.nrow(); t++){
				
				// This stays the same for all rows within this threshold:
				int tind = output_ThresholdIndex[row] - 1L;
				double th1 = thresholds(tind,0L);
				double th2 = thresholds(tind,1L);				
				
				// TODO: Method should not need to be set but do as a debug check?
				
				// Apply each test and save the results:
				output_Method[row] = "BNB";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				int fl = bnb_pval(sum_ctl, N_ctl, estk_ctl, mean_ctl, var_ctl, sum_tx, N_tx, estk_tx, mean_tx, var_tx, 0.0, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);
				row++;
				
				output_Method[row] = "Levecke";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				levecke_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "SumTx0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "WAAVP";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				waavp_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "SumTx0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "MLE";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				mle_u_ci(mean_ctl, mean_tx, var_ctl, var_tx, N_ctl, N_tx, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_tx==0) ? "SumTx0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;
				
				Rcpp::checkUserInterrupt();
				
			}
		}
	}
	
	return output;
}


// [[Rcpp::export]]
Rcpp::DataFrame fecrt_sim_paired(int iters, Rcpp::NumericVector red, int N, double mu, double k_pre, double k_post, double k_c, Rcpp::NumericMatrix thresholds, Rcpp::NumericVector conjugate_priors, int delta, int beta_iters, int approx, double tail, bool useml, Rcpp::NumericVector dobson_cl, Rcpp::NumericVector dobson_priors)
{

	// Check that the thresholds have 2 columns and >0 rows:
	if(thresholds.ncol() != 2){
		Rcpp::stop("There must be exactly 2 columns");
	}
	if(thresholds.nrow() < 1){
		Rcpp::stop("There must be 1 or more rows");
	}
	for(int i=0; i<thresholds.nrow(); i++){
		if(thresholds(i,1) < thresholds(i,0)){
			Rcpp::stop("Second thresholds must be higher than the first");
		}
		if(thresholds(i,0) < 0 || thresholds(i,0) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
		if(thresholds(i,1) < 0 || thresholds(i,1) > 1){
			Rcpp::stop("Thresholds must be between 0 and 1");
		}
	}
	
	// Check that k_c is greater than pre and post:
	if(k_c < k_pre || k_c < k_post){
		Rcpp::stop("k_c must be larger than k_pre and k_post");
	}
	
	// TODO: tidy up check and error:
	if((double(N)*double(mu)) > (0.1 * INT_MAX)){
		Rcpp::stop("Possible int overflow");
	}
	
	// Use the expand.grid R function:
	Rcpp::Function expGrid("expand.grid");
	
	Rcpp::StringVector methods = { "BNB", "Levecke", "WAAVP", "MLE", "Dobson" };
	Rcpp::IntegerVector tseq = Rcpp::seq(1L, thresholds.nrow());
	Rcpp::DoubleVector td = { NA_REAL };
	Rcpp::StringVector ts = { NA_STRING };
//	Rcpp::IntegerVector ti = { NA_INTEGER };
	Rcpp::IntegerVector iteration = Rcpp::seq(1L, iters);
	
	Rcpp::DataFrame output = expGrid(methods, red, tseq, td, td, iteration, td, td, td, ts, Rcpp::Named("stringsAsFactors")=false);
	output.names() = Rcpp::StringVector::create("Method", "Reduction", "ThresholdIndex", "ThresholdLower", "ThresholdUpper", "Iteration", "ObsReduction", "Stat1", "Stat2", "Classification");

	// Get references to list items:
	Rcpp::StringVector output_Method = output[0L];
	Rcpp::DoubleVector output_Reduction = output[1L];
	Rcpp::IntegerVector output_ThresholdIndex = output[2L];
	Rcpp::DoubleVector output_ThresholdLower = output[3L];
	Rcpp::DoubleVector output_ThresholdUpper = output[4L];
	Rcpp::IntegerVector output_Iteration = output[5L];
	Rcpp::DoubleVector output_ObsReduction = output[6L];
	Rcpp::DoubleVector output_Stat1 = output[7L];
	Rcpp::DoubleVector output_Stat2 = output[8L];
	Rcpp::StringVector output_Classification = output[9L];
	
	// Note:  Rcpp allows conversion from long long but not to double array directly
	double conjugate_priors_db[2] = { conjugate_priors[0], conjugate_priors[1] };
	double dobson_priors_db[2] = { dobson_priors[0], dobson_priors[1] };
	
	// To control the row index to write to:
	int row = 0L;
	
	// Precalculate some things:
	double muadj_pre = mu / k_pre;
	double muadj_post = mu / k_post;
	double b_pre = k_c - k_pre;
	double b_post = k_c - k_post;
	
	// First loop over correlation and pre-treatment datasets to simulate:
	for(int i=0; i<iters; i++){
		Rcpp::NumericVector gammas;
		gammas = Rcpp::rgamma(N, k_c, 1.0);
		
		Rcpp::IntegerVector pre_data(N);
		long long sum_pre = 0L;
		do {
			sum_pre = 0L;
			for(int i=0; i<N; i++){
				pre_data[i] = R::rpois(R::rbeta(k_pre, b_pre) * gammas[i] * muadj_pre);
				sum_pre += pre_data[i];
			}
		} while (sum_pre == 0L);

		double mean_pre = Rcpp::mean(pre_data);
		double var_pre = Rcpp::var(pre_data);
		
		// Then loop over reductions:
		for(int r=0; r<red.length(); r++){
			long long sum_post = 0L;
			Rcpp::IntegerVector post_data(N);
			for(int i=0; i<N; i++){
				post_data[i] = R::rpois(R::rbeta(k_post, b_post) * gammas[i] * muadj_post * red[r]);
				sum_post += post_data[i];
			}

			double mean_post = Rcpp::mean(post_data);		
			double var_post = Rcpp::var(post_data);
			double obsred = mean_post / mean_pre;
			
			// We need both scaled_cov for estimate_k and cov for the estimation functions (including bnb if using a large sample approximation)
			double scaled_cov = 0.0;
			double cov = 0.0;
		    for(int i=0; i<N; i++){
				scaled_cov += (double(pre_data[i])/mean_pre - 1.0) * (double(post_data[i])/mean_post - 1.0);
				cov += (double(pre_data[i]) - mean_pre) * (double(post_data[i]) - mean_post);
			}
		    scaled_cov = scaled_cov / double(N - 1L);
		    cov = cov / double(N - 1L);
			
			Rcpp::NumericVector ks;
			if(useml){
				ks = estimate_k_ml(pre_data, mean_pre, var_pre, post_data, mean_post, var_post, cov, true);
			}else{
				ks = estimate_k(mean_pre, var_pre, mean_post, var_post, cov, true);
			}
			// Test:
			if(!R_finite(ks[0])){
				Rcpp::stop("Non-finite k1 generated");
			}
			if(ks[0] <= 0.0){
				Rcpp::stop("k1 below 0.0 generated");
			}
			if(!R_finite(ks[1])){
				Rcpp::stop("Non-finite k2 generated");
			}
			if(ks[1] <= 0.0){
				Rcpp::stop("k2 below 0.0 generated");
			}
			double estk_pre = ks[0];
			double estk_post = ks[1];
			
			// For cheating:
			//estk_pre = k_pre / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));
			//estk_post = k_post / (1.0 - (std::sqrt(k_pre) * std::sqrt(k_post) / k_c));

			Rcpp::NumericVector uks = estimate_k(mean_pre, var_pre, mean_post, var_post, 0.0, false);
			double estk_pre_u = uks[0];
			double estk_post_u = uks[1];

//			Rcpp::Rcout << cov << " - " << estk_pre << " - " << estk_pre_u << " - " << estk_post << " - " << estk_post_u << "\n";
			
			// Then loop over thresholds:
			for(int t=0; t<thresholds.nrow(); t++){
				
				// This stays the same for all rows within this threshold:
				int tind = output_ThresholdIndex[row] - 1L;
				double th1 = thresholds(tind,0L);
				double th2 = thresholds(tind,1L);
				
				
				// TODO: Method should not need to be set but do as a debug check?
				
				// Apply each test and save the results:
				
				output_Method[row] = "BNB";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				int fl = bnb_pval(sum_pre, N, estk_pre, mean_pre, var_pre, sum_post, N, estk_post, mean_post, var_post, cov, 1.0, th1, th2, conjugate_priors_db, delta, beta_iters, approx, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = get_type_pv(1.0-obsred, output_Stat1[row], output_Stat2[row], tail, th1, th2);
				row++;
				
				output_Method[row] = "Levecke";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				levecke_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "WAAVP";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				waavp_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;

				output_Method[row] = "MLE";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				mle_p_ci(mean_pre, mean_post, var_pre, var_post, cov, N, tail, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post==0) ? "SumPost0" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;
	
				output_Method[row] = "Dobson";
				output_ObsReduction[row] = obsred;
				output_ThresholdLower[row] = th1;
				output_ThresholdUpper[row] = th2;
				dobson_ci(sum_pre, sum_post, dobson_cl[0], dobson_cl[1], dobson_priors_db, &output_Stat1[row], &output_Stat2[row]);
				output_Classification[row] = (sum_post>sum_pre) ? "post>pre" : get_type_ci(1.0-obsred, output_Stat1[row], output_Stat2[row], th1, th2);
				row++;
					
				Rcpp::checkUserInterrupt();
				
			}
		}
	}

	return output;
}

