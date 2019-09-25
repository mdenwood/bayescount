#' @name count_power
#' @aliases count.power fec.power FEC.power
#' @title Count data power calculations
#'
#' @description
#' Power calculation based on BNB code and christian approach yet to be written.  Do either test vs known single mean, or test for difference in mean to observed dataset (k can vary between datasets or be fixed).
#'
#' @seealso \code{\link{count_analysis}}, \code{\link{reduction_power}}, \code{\link{bayescount}}



### Replace with a new function based on new methods
### Preserve old arguments used in version 0.999 with warning

# STUFF TO BE REMOVED AT THE BOTTOM OF THIS FILE

count_power <- function(type='fec', cutoff=if('fec') 400 else 95, alternative='two.sided', conf.level=0.05, maxiterations=1000000, mc.decimals=2, mc.conf = 0.99, feedback=FALSE, ...){

	# valid passthroughs:  meanepg=200, reduction = 80, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, 
	# Check passthrough and stop if any others given


	# All this crap can be handeled by the fecrt.precision function so remove from here

	if(mc.conf >= 1) stop("mc.conf must be < 1")
	conf <- (1-mc.conf)/2
	lci <- 0+conf
	uci <- 1-conf

	if(feedback){
	if(any(.Platform$GUI == c("AQUA", "Rgui"))){
	#feedback <- FALSE
	warning("Printing the progress of the function using the GUI versions of R may massively increase the time taken, I suggest setting feedback=FALSE or using a command line version of R instead")
	}
	##### CHECK ON WINDOWS / UNIX
	}

	if(replicates < 1 | animals < 1) stop("Specified values for animals and replicates must be greater than 0")

	if(animals==1 & true.sample==FALSE) warning("NOTE:  The confidence calculated is for the population from which the animal is derived, to calculate the confidence for the true mean of this individual use true.sample=TRUE")

	lowerl <- (100-upper.interval)/100
	upperl <- (100-lower.interval)/100
	delta <- (100-reduction)/100

	if(lowerl < 0 | upperl < 0 | delta < 0 | lowerl > 1 | upperl > 1 | delta > 1) stop("reduction, lower.interval and upper.interval must all be between 0 and 100(%)")


	# Use fecrt.precision(type='mc') to get the mc output we need, then fecrt.precision(type='precision') to get power

	#	meanepg, reduction = 80, g.faeces=3, sensitivity=1/25, replicates=1, animals=10, pre.coeffvarrep=0.4, pre.coeffvarind=0.3, pre.coeffvargroup=0.7, post.coeffvarrep=0.4, post.coeffvarind=0.3, post.coeffvargroup=0.7, true.sample=FALSE, lower.interval = 0, upper.interval=90, maxiterations=1000000, mc.decimals=2, mc.conf = 0.99, feedback=FALSE
		p <<- fecrt.confidence(mean, thresh, g.faeces, sens, replicates, ani, coeffr, coeffi, cvg, coeffr, coeffi, cvg, FALSE, mc.decimals=NA, mc.output=TRUE)
	
		qm <- quantile((1-p)*100, 0.05)

	#	q <- fecrt.confidence(mean, qm, g.faeces, sens, replicates, ani, coeffr, coeffi, cvg, coeffr, coeffi, cvg, FALSE, mc.decimals=NA)
	
	#	pr <- quantile((1-q)*100, 0.95)
	
		confidence <<- fecrt.confidence(mean, red, g.faeces, sens, replicates, ani, coeffr, coeffi, cvg, coeffr, coeffi, cvg, TRUE, upper.interval=qm)$ci[2]
	
#		if(confidence>0.8) break
	#}

	results <- rbind(results, c(thresh, sens, ani, red, replicates))

	print(nrow(results))

	#}}}}


	#save(results, file=new_unique('fecrt confidence results', '.Rsave'))
	
	
	#}



	return(output)
}


count.power <- count_power
fec.power <- count_power
FEC.power <- count_power
