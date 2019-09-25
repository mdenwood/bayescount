estimate_ks <- function(data_1, data_2){
	
	# Working that justifies the equations:
	if(FALSE){

		# Equating kc to correlation:
		# Three types of k:  ko = overall, kc=correlated fraction, ku=uncorrelated fraction, c=correlation
		# ku = m^2 / (v*(1-c)) = ko * 1/(1-c)
		# c = 1-(ko/ku)

		# bm = k1/(k1+kc-k1) = k1/kc
		# bv = (k1*(kc-k1)) / ( (k1+(kc-k1))^2 * (k1+(kc-k1)+1) )
		# ku = bm^2 / bv = (ko/kc)^2 / ( (ko*(kc-ko)) / ( (ko+(kc-ko))^2 * (ko+(kc-ko)+1) ) )
		# Agrees with:
		# ko = (kc * ku) / (kc + ku + 1)
		# ku = ((kc+1)*ko)/(kc-ko)

		# Total variance of product of gammas is:
		# mu^2 / (kc * ku) / (kc + ku + 1)
		ka <- runif(1,0.1,5); kb <- runif(1,0.1,5)
		var(rgamma(it,ka,scale=1/ka) * rgamma(it, kb, scale=1/kb))
		(ka+kb+1) / (ka*kb)

		# So correlation = proportion of correlated variance = (mu^2 / kc) / (mu^2 / (kc * ku) / (kc + ku + 1))
		# Or for known ko:  (mu^2 / kc) / (mu^2 / ko) = ko/kc

		mu <- rgamma(1,1,scale=10)
		red <- runif(1,0.1,2)
		# Types of ko:
		k1 <- runif(1,0.1,10)
		k2 <- runif(1,0.1,10)
		# kc:
		kc <- max(k1,k2)+runif(1,0,5)

		it <- 10^5
		gammas <- rgamma(it, kc, 1)
		d1 <- rbeta(it, k1, kc-k1)*gammas*mu/k1
		d2 <- rbeta(it, k2, kc-k2)*gammas*mu/k2*red
		cor(d1,d2)

		# We need to scale the respective d1 and d2 variances to a mean of 1:
		tcov <- cor(d1,d2) * (sd(d1)/mean(d1)) * (sd(d2)/mean(d2))
		# Simplify:
		# = cov(d1,d2)/(sd(d1)*sd(d2)) * (sd(d1)/mean(d1)) * (sd(d2)/mean(d2))
		# = cov(d1,d2) / (mean(d1) * mean(d2))

		tcov <- cov(d1,d2) / (mean(d1) * mean(d2))
		# Correlation:
		k1/kc
		tcov / (sd(d1)/mean(d1))^2
		k2/kc
		tcov / (sd(d2)/mean(d2))^2
		# Correlated k:
		kc
		k1 * (sd(d1)/mean(d1))^2 / tcov
		k2 * (sd(d2)/mean(d2))^2 / tcov
		# Uncorrelated k1:
		((kc+1)*k1)/(kc-k1)
		(tcov*mean(d1)^2 + k1*var(d1)) / (var(d1) - tcov*mean(d1)^2)
		# Uncorrelated k2:
		((kc+1)*k2)/(kc-k2)
		(tcov*mean(d2)^2 + k2*var(d2)) / (var(d2) - tcov*mean(d2)^2)

		
		# Reduces to in special cases:
		if(red==1){
			k1/kc
			cov(d1,d2) / var(d1)
			k2/kc
			cov(d1,d2) / var(d2)
	
			if(k1==k2){
				k1/kc
				cor(d1,d2)
			}
		}
		
		# This appears to be how to correct k for BNB though:
		mean(d1)^2 / (var(d1)*(1-cor(d1,d2)))
		k1 / (1-cor(d1,d2))
		# To estimate covariance:
		(mean(d1) * mean(d2)) * (sd(d1)/mean(d1))^2 * k1/kc
		(mean(d1) * mean(d2)) * (sd(d2)/mean(d2))^2 * k2/kc
		cov(d1,d2)
		# Then to estimate correlation:
		(mean(d1)/sd(d1) * mean(d2)/sd(d2)) * (sd(d1)/mean(d1))^2 * k1/kc
		(mean(d1)/sd(d1) * mean(d2)/sd(d2)) * (sd(d2)/mean(d2))^2 * k2/kc
		# Both simplify to:
		sqrt(k1)*sqrt(k2) * 1/kc
		cor(d1,d2)
		k1 / (1- sqrt(k1)*sqrt(k2)/kc)
		# Or to find kc from estimates of k1 and k2 and correlation:
		sqrt(k1*k2) / cor(d1,d2)
		kc
		
		

		# Adding Poisson difficulties actually just cancels out:
		d1 <- rpois(it, rbeta(it, k1, kc-k1)*gammas*mu/k1)
		d2 <- rpois(it, rbeta(it, k2, kc-k2)*gammas*mu/k2*red)

		# Simplifying proof:
		cor12 <- cov(d1,d2) / (sd1*sd2)
		cor12
		sd1 <- sqrt(var(d1)-mean(d1))
		sd2 <- sqrt(var(d2)-mean(d2))
		tcov <- cov(d1,d2) / (sd1*sd2) * (sd1/mean(d1)) * (sd2/mean(d2))
		# = cov(d1,d2) / (mean(d1) * mean(d2))
		# / simplifying

		tcov <- cov(d1,d2) / (mean(d1) * mean(d2))
		# Correlation:
		k1/kc
		tcov / (sqrt(var(d1)-mean(d1))/mean(d1))^2
		k2/kc
		tcov / (sqrt(var(d2)-mean(d2))/mean(d2))^2
		# Correlated k:
		kc
		k1 * (sqrt(var(d1)-mean(d1))/mean(d1))^2 / tcov
		k2 * (sqrt(var(d2)-mean(d2))/mean(d2))^2 / tcov
		# Uncorrelated k1:
		((kc+1)*k1)/(kc-k1)
		(tcov*mean(d1)^2 + k1*(var(d1)-mean(d1))) / ((var(d1)-mean(d1)) - tcov*mean(d1)^2)
		# Uncorrelated k2:
		((kc+1)*k2)/(kc-k2)
		(tcov*mean(d2)^2 + k2*(var(d2)-mean(d2))) / ((var(d2)-mean(d2)) - tcov*mean(d2)^2)

	}
	
	d1 <- data_1
	d2 <- data_2
	
	tcov <- cov(d1,d2) / (mean(d1) * mean(d2))
	# Correlation:
	correlation <- c(tcov / (sqrt(var(d1)-mean(d1))/mean(d1))^2, tcov / (sqrt(var(d2)-mean(d2))/mean(d2))^2)
	# Total k:
	effvar1 <- var(d1)-mean(d1)
	effvar2 <- var(d2)-mean(d2)
	k1 <- mean(d1)^2 / effvar1
	k2 <- mean(d2)^2 / effvar2
	total_k <- c(k1,k2)
	# Correlated k:
	correlated_k <- c(k1 * (sqrt(var(d1)-mean(d1))/mean(d1))^2 / tcov, k2 * (sqrt(var(d2)-mean(d2))/mean(d2))^2 / tcov)
	# Uncorrelated k:
	uncorrelated_k <- c((tcov*mean(d1)^2 + k1*(var(d1)-mean(d1))) / ((var(d1)-mean(d1)) - tcov*mean(d1)^2), (tcov*mean(d2)^2 + k2*(var(d2)-mean(d2))) / ((var(d2)-mean(d2)) - tcov*mean(d2)^2))
	
	return(list(correlation=correlation, total_k=total_k, correlated_k=correlated_k, uncorrelated_k=uncorrelated_k, cfun=estimate_k(as.double(mean(d1)), as.double(var(d1)), as.double(mean(d2)), as.double(var(d2)), as.double(cov(d1, d2)), as.logical(TRUE))))
	
}