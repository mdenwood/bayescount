#' @name reduction_precision
#' @aliases reduction.precision fecrt.precision FECRT.precision fecrt.power.limits FECRT.power.limits
#' @title Count data reduction precision calculations
#'
#' @description
#' Precision analysis by Monte Carlo simulation and examination of distribution of mean
#'
#' @seealso \code{\link{reduction_analysis}}, \code{\link{reduction_power}}, \code{\link{bayescount}}


### Preserve old c code

reduction_precision <- function(meanepg=200, reduction = 95, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, lower.interval=NA, upper.interval=NA, iterations=100000, mc.decimals = 0.95, mc.conf = 0.99, feedback=FALSE, mc.output=FALSE){
	
	# mc.decimals is the required mc.decimals, mc.conf is the mc.conf in this mc.decimals when using simulations
	if(mc.decimals >= 1) stop("Required mc.decimals must be < 1")

	if(mc.conf >= 1) stop("mc.conf must be < 1")
	conf <- (1-mc.conf)/2
	lci <- 0+conf
	uci <- 1-conf

	if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")

	if(!is.na(lower.interval) & !is.na(upper.interval)) stop("One or both of lower.interval or upper.interval must be non-fixed (NA)")

	fix.upper <- !is.na(lower.interval)
	fix.lower <- !is.na(upper.interval)

	fixed.lower = lowerl <- (100-upper.interval)/100
	fixed.upper = upperl <- (100-lower.interval)/100
	delta <- (100-reduction)/100

	tlowerl <- lowerl
	tupperl <- upperl
	tdelta <- delta
	if(is.na(tlowerl)) tlowerl <- 0.5
	if(is.na(tupperl)) tupperl <- 0.5
	if(is.na(tdelta)) tdelta <- 0.5

	if(tlowerl < 0 | tupperl < 0 | tdelta < 0 | tlowerl > 1 | tupperl > 1 | tdelta > 1) stop("reduction, lower.interval and upper.interval (if supplied) must all be between 0 and 100(%)")

	target <- confidence

	if(feedback){
	if(any(.Platform$GUI == c("AQUA", "Rgui"))){
	#feedback <- FALSE
	warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
	}
	##### CHECK ON WINDOWS / UNIX
	}


	if(animals==1 & true.sample==FALSE) warning("NOTE:  The confidence calculated is for the population from which the animal is derived, to calculate the confidence for the true mean of this individual use true.sample=TRUE")

	lowerl <- lower.interval
	upperl <- upper.interval

	delta <- (100-reduction)/100

	if(delta < 0 | delta > 1) stop("Supplied value for 'reduction' must be between 0 and 100(%)")

	start <- Sys.time()

	if(true.sample){
	
		out <- .C_OLD(C_precision_reduction, as.numeric(premean), as.numeric(reduction), as.integer(replicates), as.integer(animals), as.numeric(pre.coeffvarrep), as.numeric(pre.coeffvarind), as.numeric(pre.coeffvargroup), as.numeric(post.coeffvarrep), as.numeric(post.coeffvarind), as.numeric(post.coeffvargroup), as.integer(iterations), as.integer(feedback), numeric(iterations))

		lo <- length(out)
		meanreds <- out[[lo]]
	
	}else{
		stop()

		lo <- length(out)
		meanreds <- out[[lo]]
	
	}

	mcs <- meanreds

	f <- function(limit){
		prob <- limit
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		limits <- quantile(mcs, probs=c(0+((1-prob)/2), 1-((1-prob)/2)))
		nin <- sum(mcs <= limits[2] & mcs >= limits[1])
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}

	if(fix.lower){
	f <- function(limit){
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		nin <- sum(mcs <= limit & mcs >= fixed.lower)
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}
	}
	if(fix.upper){
	f <- function(limit){
		#limits <- HPDinterval(mcmc(mcs), prob=prob)[1,]
		nin <- sum(mcs <= fixed.upper & mcs >= -limit)
		nout <- length(mcs)-nin
		med <- qbeta(0.5, nin+1, nout+1)
		return(med)
	}
	}

	limits <- c(0,1)
	if(fix.upper) limits <- c(-1,0)
	if(fix.lower) limits <- c(0,1)

	bsres <- binary.search(f, target, limits)
	if(bsres$status!="OK"){
		if(bsres$status!="Absolute value not possible but this is closest"){
			limits <- c(0,1)
			if(fix.upper) limits[2] <- fixed.upper
			if(fix.lower) limits[1] <- fixed.lower
			nin <- sum(mcs <= limits[1] & mcs >= limits[2])
			nout <- length(mcs)-nin
			med <- round(qbeta(0.5, nin+1, nout+1), 2)
			stop(paste("It is not possible to achieve the desired confidence with any tolerance; maximum confidence =", med))
		}
	}

	if(fix.upper) limits <- c(-bsres$value, fixed.upper)
	if(fix.lower) limits <- c(fixed.lower, bsres$value)
	if(!fix.upper & !fix.lower) limits <- quantile(mcs, prob=c(0+((1-bsres$value)/2), 1-((1-bsres$value)/2)))

	nin <- sum(mcs <= limits[2] & mcs >= limits[1])
	nout <- length(mcs)-nin

	#limits <- paste(prettyround((100*(1-limits))[2:1], 2), "%", sep="")
	limits <- (100*(1-limits))[2:1]

	confidence <- qbeta(c(lci,0.5,uci), nin+1, nout+1)

	names(limits) <- c("lower.interval", "upper.interval")
	names(confidence) <- c(paste("lower.", mc.conf*100, "%.estimate", sep=""), "median.estimate", paste("upper.", mc.conf*100, "%.estimate", sep=""))

	if(mc.output) return(list(limits=limits, confidence=confidence, mc.output=meanreds)) else return(list(limits=limits, confidence=confidence))


}

reduction.precision <- reduction_precision
fecrt.precision <- reduction_precision
FECRT.precision <- reduction_precision

fecrt.power.limits <- function(...){
  warning('Use of the fecrt.power.limits function is deprecated and will be removed in version 1.1 - please use the count_precision function instead')
  return(reduction_precision(...))
}
FECRT.power.limits <- fecrt.power.limits
