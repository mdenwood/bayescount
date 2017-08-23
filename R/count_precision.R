#' @name count_precision
#' @aliases count.precision fec.precision FEC.precision fec.power.limits FEC.power.limits
#' @title Count data precision calculations
#'
#' @description
#' Precision analysis by Monte Carlo simulation and examination of distribution of mean
#'
#' @seealso \code{\link{count_analysis}}, \code{\link{count_power}}, \code{\link{bayescount}}


### Preserve old c code

# Thoroughly test everything against nb function.  g.faeces, coeffs (not group) can be vectorised if not using NB.  pre and post for faces, sens and reps


# From help file for version 0.99:
#count.precision(meanepg=200, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, coeffvarrep=0.4, coeffvarind=0.3, coeffvargroup=0.7, true.sample=FALSE, lower.limit=NA, upper.limit=NA, iterations=100000, power = 0.95, confidence = 0.99, feedback=FALSE, forcesim=FALSE)
   
count_precision <- function(meanepg=200, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, coeffvarmcm=0.1, coeffvarrep=0.4, coeffvarind=0.3, coeffvargroup=0.7, true.sample=FALSE, type = 'confidence', precision=if(type=='confidence') 0.1 else NA, lower.interval=meanepg*(1-precision), upper.interval=meanepg*(1+precision), confidence = 0.95, iterations=if(type=='confidence') 10^6 else 10^5, mc.decimals=2, mc.conf = 0.99, feedback=FALSE, forcesim=FALSE){
	
	if(true.sample)
		warning('Deprecated and ignored')
	
	type <- tolower(type)
	if(length(type)==2){
		if(any(type=='interval') & any(type=='mc')) type <- 'interval.mc'
	}
	if(length(type)!=1 | !any(type==c('interval', 'confidence', 'mc', 'interval.mc'))) stop('The specified type must be one of "interval", "confidence", or "mc" (or "interval.mc" for interval calculation with the Monte carlo sample returned)')

	# Make sure we don't set things that don't make any sense for the given type - give a warning if non-default options provided that wont be used and fail if things are set that need not to be
	if(type=='interval'){
		if(!is.na(lower.interval) & !is.na(upper.interval)) stop("One or both of lower.interval or upper.interval must be non-fixed (NA) when using the 'interval' type")
		fixedc <- TRUE
	}
	if(type=='interval.mc'){
		if(!is.na(lower.interval) & !is.na(upper.interval)) stop("One or both of lower.interval or upper.interval must be non-fixed (NA) when using the 'interval.mc' type")
		fixedc <- TRUE
		forcesim <- TRUE
	}
	if(type=='confidence'){
		if(is.na(lower.interval) | is.na(upper.interval)) stop("The precision (or lower.interval and upper.interval) must be supplied when using the 'confidence' type")
		if(formals(fec.precision)$confidence != confidence) warning("Any supplied value for 'confidence' is ignored when using the 'confidence' type")
		if(is.na(mc.decimals)) fixedc <- TRUE else fixedc <- FALSE
	}
	if(type=='mc'){
		if(formals(fec.precision)$confidence != confidence | !is.na(lower.interval) | !is.na(upper.interval) | formals(fec.precision)$mc.conf != mc.conf | formals(fec.precision)$mc.decimals != mc.decimals) warning("Any supplied value for 'confidence', 'precision', 'lower.interval', 'upper.interval', 'mc.conf' or 'mc.decimals' is ignored when using the 'mc' type")
		fixedc <- TRUE
		forcesim <- TRUE
	}

	if(is.na(as.integer(iterations))) stop('You must supply an integer for "iterations"')
	if(is.na(as.logical(feedback))) stop('You must supply a logical value for "feedback"')
	if(any(is.na(c(meanepg, g.faeces, sensitivity, replicates, animals, coeffvarmcm, coeffvarrep, coeffvarind, coeffvargroup)))) stop('You must supply numeric values for meanepg, g.faeces, sensitivity, replicates, animals, coeffvarmcm, coeffvarrep, coeffvarind, and coeffvargroup')
	if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")

	if(type!='mc'){		

		# confidence is the required confidence, mc.conf is the confidence in this confidence when using simulations
		if(confidence >= 1) stop("Required confidence must be < 1")

		if(mc.conf >= 1) stop("Confidence must be < 1")
		lci <- 0+((1-mc.conf)/2)
		uci <- 1-((1-mc.conf)/2)
	
	}

	# If we have 1 animal and true.sample=T then we can forget about cv.group and use nbconfidence for confidence type
	if(animals==1 & true.sample==TRUE) coeffvargroup <- 0#10^-10
	if(coeffvargroup > (coeffvarrep+coeffvarind)/10) approximate <- FALSE else approximate <- TRUE
	if(forcesim) approximate <- FALSE

	if(feedback & !approximate){
	if(any(.Platform$GUI == c("AQUA", "Rgui"))){
	#feedback <- FALSE
	warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
	}
	##### CHECK ON WINDOWS / UNIX
	}

	if(animals==1 & true.sample==FALSE) warning("NOTE:  The confidence calculated is for the population from which the animal is derived, to calculate the confidence for the true mean of this individual use true.sample=TRUE")


	if(approximate){
		# if any lengths > 1 then warn we can't approximate
	}else{
		if(length(g.faeces)==1) g.faeces <- rep(g.faeces, animals)
		if(length(sensitivity)==1) sensitivity <- rep(sensitivity, animals)
		if(length(replicates)==1) replicates <- rep(replicates, animals)
		if(length(coeffvarmcm)==1) coeffvarmcm <- rep(coeffvarmcm, animals)
		if(length(coeffvarrep)==1) coeffvarrep <- rep(coeffvarrep, animals)
		if(length(coeffvarind)==1) coeffvarind <- rep(coeffvarind, animals)
		if(length(coeffvargroup)==1) coeffvargroup <- rep(coeffvargroup, animals)

		print('make sure lengths are all OK')
	}

	fix.lower <- !is.na(lower.interval)
	fix.upper <- !is.na(upper.interval)

	fixed.lower <- lower.interval
	fixed.upper <- upper.interval

	target <- confidence


	start <- Sys.time()
	output <- list()




	# First run any simulations that need to be done:

	if(!approximate){

		if(fixedc){
			if(true.sample){
				stop()
		
			}else{
				mcout <- .C(C_precision_count, as.numeric(meanepg), as.numeric(g.faeces), as.numeric(sensitivity), as.integer(replicates), as.integer(animals), as.numeric(coeffvarmcm), as.numeric(coeffvarrep), as.numeric(coeffvarind), as.numeric(coeffvargroup), as.integer(iterations), as.integer(feedback), as.numeric(numeric(iterations)))

			}

			meancounts <- mcout[[length(mcout)]]
			nin <- sum(meancounts >= lower.interval & meancounts <= upper.interval)
			nout <- length(meancounts)-nin

		}else{
			if(true.sample){
				stop()
			}else{
				stop()
			}
		
			nin <- confout[[length(confout)-1]]
			nout <- confout[[length(confout)]] - nin
		}
	
		if(type=='confidence'){
			output=c(output, list(roundedci=round(qbeta(c(lci, 0.5, uci), nin+1, nout+1), mc.decimals), ci=qbeta(c(lci, 0.5, uci), nin+1, nout+1), within=nin, without=nout, total=nin+nout))
		
			names(output$roundedci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))
			names(output$ci) <- c(paste("lower.", confidence*100, "%", sep=""), "median", paste("upper.", confidence*100, "%", sep=""))
		}
	
		# And in case we need the optimising functions:
	
		f <- function(limit){
			prob <- limit
			#limits <- HPDinterval(mcmc(meancounts), prob=prob)[1,]
			limits <- quantile(meancounts, probs=c(0+((1-prob)/2), 1-((1-prob)/2)))
			nin <- sum(meancounts <= limits[2] & meancounts >= limits[1])
			nout <- length(meancounts)-nin
			med <- qbeta(0.5, nin+1, nout+1)
			return(med)
		}

		if(fix.lower){
		f <- function(limit){
			#limits <- HPDinterval(mcmc(meancounts), prob=prob)[1,]
			nin <- sum(meancounts <= limit & meancounts >= fixed.lower)
			nout <- length(meancounts)-nin
			med <- qbeta(0.5, nin+1, nout+1)
			return(med)
		}
		}
		if(fix.upper){
		f <- function(limit){
			#limits <- HPDinterval(mcmc(meancounts), prob=prob)[1,]
			nin <- sum(meancounts <= fixed.upper & meancounts >= -limit)
			nout <- length(meancounts)-nin
			med <- qbeta(0.5, nin+1, nout+1)
			return(med)
		}
		}
	
	}




	# Now do any approximate stuff that needs to be done:

	if(approximate){
	
		if(true.sample==FALSE | animals > 1){
			# If we're using this with animals > 1 then cvgroup must be so small compared to cvind and cvrep that we can just pretend that we're drawing all replicates from the same animal (with a slightly larger cv and animalmean=populationmean)
			# If we're using this with true.sample=FALSE then we need to increase the variability so that we're representing the group mean that the (probably single) animal is taken from
			coeffvarind <- sqrt(coeffvarind^2 + coeffvargroup^2 + coeffvarind^2*coeffvargroup^2)
			replicates <- replicates*animals
		}

		# Produces the same results as precisionanalysispopulation (only for animal=1 and/or cv.group<<cv.other)

		coeff.var <- sqrt(coeffvarmcm^2 + (coeffvarrep^2)/g.faeces + (coeffvarmcm^2 * (coeffvarrep^2)/g.faeces))
		coeff.var <- sqrt(coeff.var^2 + coeffvarind^2 + (coeff.var^2 * coeffvarind^2))
	
		eggs.counted <- replicates*meanepg*sensitivity

		shape <- replicates / coeff.var^2

		if(type=='confidence'){
			upper.tol <- replicates*sensitivity*(upper.interval)
			lower.tol <- replicates*sensitivity*(lower.interval)
	
			nbconfidence <- pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted))
		
			output = c(output, list(roundedci=replicate(3, round(nbconfidence, mc.decimals)), ci=replicate(3, nbconfidence)))
		
			names(output$roundedci) <- c("absolute.value", "absolute.value", "absolute.value")
			names(output$ci) <- c("absolute.value", "absolute.value", "absolute.value")
		
		}else{
		
			# The optimising functions:
	
			f <- function(limit){
				precision <- limit

				upper.tol <- replicates*sensitivity*(meanepg+meanepg*precision)
				lower.tol <- replicates*sensitivity*(meanepg-meanepg*precision)

				return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))

			}	
			if(fix.lower){
			f <- function(limit){
				upper.tol <- replicates*sensitivity*(limit)
				lower.tol <- replicates*sensitivity*(fixed.lower)

				return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))
			}
			}
			if(fix.upper){
			f <- function(limit){
				upper.tol <- replicates*sensitivity*(fixed.upper)
				lower.tol <- replicates*sensitivity*(-limit)

				return(pnbinom(upper.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) - pnbinom(lower.tol, size=shape, mu=eggs.counted, lower.tail=TRUE) + suppressWarnings(dnbinom(lower.tol, size=shape, mu=eggs.counted)))
			}
			}
		}
	
	}



	#  Now optimise if we need to (either MC output or approximate):

	if(any(type==c('interval.mc', 'interval'))){
		limits <- c(0,1)
		if(fix.upper) limits <- c(-Inf,0)
		if(fix.lower) limits <- c(0,Inf)

		bsres <- binary.search(f, target, limits)
		if(bsres$status!="OK"){
			if(bsres$status=="Absolute value not possible but this is closest"){
				status <- 'Closest possible'
			}else{
				status <- bsres$status
			}
		}else{
			status <- 'Exact'
		}
	
		if(is.na(bsres$value)){
			limits <- c(NA,NA)
			nin <- NA
			nout <- NA
		}else{

			if(fix.upper) limits <- c(-bsres$value, fixed.upper)
			if(fix.lower) limits <- c(fixed.lower, bsres$value)
			if(!fix.upper & !fix.lower) limits <- quantile(meancounts, prob=c(0+((1-bsres$value)/2), 1-((1-bsres$value)/2)))
	
			nin <- sum(meancounts <= limits[2] & meancounts >= limits[1])
			nout <- length(meancounts)-nin
		}

		if(!approximate){
			confidencelimits <- qbeta(c(lci,0.5,uci), nin+1, nout+1)
			names(limits) <- c("lower.interval", "upper.interval")
			names(confidencelimits) <- c(paste("lower.", confidence*100, "%.estimate", sep=""), "median.estimate", paste("upper.", confidence*100, "%.estimate", sep=""))
		}else{
			confidencelimits <- replicate(3, bsres$objective)
			names(limits) <- c("lower.interval", "upper.interval")
			names(confidencelimits) <- c("absolute.value", "absolute.value", "absolute.value")
		}

		output <- c(output, list(limits=limits, confidence=confidencelimits, status=status))
	}


	if(any(type==c('interval.mc', 'mc'))){
		output <- c(output, list(mc.output=meancounts))
	}

	time <- timestring(start, Sys.time(), units="s", show.units=FALSE)
	output <- c(output, list(time.taken=time))


	return(output)

}


count.precision <- count_precision
fec.precision <- count_precision
FEC.precision <- count_precision

fec.power.limits <- function(...){
  warning('Use of the fec.power.limits function is deprecated and will be removed in version 1.1 - please use the count_precision function instead')
  return(count_precision(...))
}
FEC.power.limits <- fec.power.limits
