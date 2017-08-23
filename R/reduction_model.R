#' @name reduction_model
#' @aliases reduction.model fecrt.model FECRT.model
#' @title Generate an un-run MCMC model for FECRT data

#' @description
#' Does not require Just Another Gibbs Sampler (JAGS) (see http://mcmc-jags.sourceforge.net) but not much use otherwise
#'
#' @seealso \code{\link{reduction_analysis}}, \code{\link{count_model}}, \code{\link{bayescount}}, \code{\link[runjags]{autoextend.jags}}, \code{\link[runjags]{runjags}}

#' This function applies a Bayesian [zero-inflated] gamma Poisson (negative binomial) model to faecal egg count reduction test (FECRT) data to return possible values for the mean drug efficacy.  Pre-treatment data are assumed to arise from either a gamma-Poisson or zero-inflated gamma-Poisson distribution, with post-treatment data described by a separate gamma-Poisson or zero-inflated gamma-Poisson distribution.  The change in mean between these distributions is therefore the mean egg count reduction.  If paired.model=TRUE, a slightly different formulation is used whereby the observed count pre-treatment is assumed to follow a compound gamma-gamma-Poisson distribution with the variability within and between animals separated.  The post treatment mean for each animal is derived from the pre-treatment animal mean and FEC reduction.  This formulation allows data with non-random missing post-treatment counts to be analysed correctly, and also allows data with repeat counts from an individual to be analysed - providing a method of increasing the power of the method substantially.  Results are also obtained using non-parametric bootstrapping and the method advocated in the 1992 World Association for the Advancement of Veterinary Parasitology (W.A.A.V.P.) methods for the detection of anthelmintic resistance in nematodes of veterinary importance (unless the data contains repeat values or missing values).  Confidence intervals for the relevant statistics are printed to file if write.file = TRUE, and returned using a custom print method.
#'
#' Lower level running of JAGS and assessing the simulation for convergence and required run length is done using  \code{\link[runjags]{autorun.jags}}.  Note: The GUI interface for R in Windows may not continually refresh the output window, making it difficult to track the progress of the simulation (if silent.jags is FALSE).  To avoid this, you can run the function from the terminal version of R (located in the Program Files/R/bin/ folder).  
#'
#' *THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*

#' @details
#' Should probably discuss different methods here

#' @keywords models

#' @return
#' Returns a list of the results obtained using each method, printed using a custom print method.

#' @references M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent. Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009

#' @examples
#' reduction_model()

#' @param name a name for the analysis (character).  Missing by default (function will require it to be input).

#' @param pre.data the pre-treatment data.  Either a numeric vector, amatrix with repeated McMasters counts from the same sample (each row) indifferent columns, or an array with different animals in dimension 1,repeat samples from the same animal in dimension 2, repeated McMasterscounts from the same sample in dimension 3 (can be of length 1 if only 1count recorded).  If an array, then the use.paired must be TRUE. Nodefault.  Ignored if a value is specified for data.

#' @param post.data the post-treatment data, for acceptable formats seepre.data.

#' @param data either a path to a comma delimited csv file, or an existingR object containing the data.  Data can be specified in one of thefollowing ways.  If a matrix, the first column is pre-treatment data,the second column is the post-treatment data, and the third column (ifsupplied) indicates control (1) or treatment (0) animal.  Alternatively,the first column is animal names id animal.names = TRUE, with columns 2and 3 making up the pre- and post-treatment data and optionally column 4the treatment of control status.  If the data is specified as a list,then the first element is taken as the pre-treatment data, and thesecond element is taken as the post-treatment data.  These elements ofthe list can be provided as for pre.data and post.data.  If the data isin a comma delimited file, it must be in the format specified for asingle matrix.  Missing data, or unused elements of non-ragged arrays,may be represented using NA, which will be removed from the data beforeanalysis.  This argument is taken from the specified values of pre.dataand post.data by default.  If a value is specified for data thenarguments specified for pre.data and post.data are ignored.

#' @param animal.names either a character vector of names to be used foreach animal specified in the data, TRUE or FALSE. If TRUE, then animalnames are taken from the first column of the matrix supplied for data (amatrix must be supplied for data if animal.names is TRUE).  If FALSE,then animal names are generated automatically.  Default FALSE.

#' @param efficacy the target \% efficacy of the drug used.  Used tocalculate the probability of resistance with the MCMC and bootstrapmethods.  Default 95\%.

#' @param confidence the degree of confidence required with which to reportthe confidence limits for the true FEC reduction.

#' @param restrict.efficacy option to restrict the estimate of efficacyfrom the MCMC method to values of greater than 0 (i.e. the mean faecalegg count cannot increase after treatment).  If TRUE, then the prior forchange in mean egg count is uniform between 0 and 1, and if FALSE thenthe prior for change in mean egg count is uniform between 0 and 10.  Ifusing control animals, then restrict.efficacy applies to the change inmean of the treated animals only, allowing the efficacy to be below 0. Default TRUE.

#' @param control.animals indication of which animals are to be used ascontrols.  Should be either a vector of TRUE/FALSE (or 1/0)corresponding to whether each animal is a control (TRUE) or treatment(FALSE), or simply FALSE (the default) in which case all animals areassumed to be treated.  Ignored if data is specified as a matrix (or csvfile) with 3 columns, in which case the third column should reflect thetreatment status.  Default FALSE.

#' @param paired.model the FECRT can be run using two compound gammadistributions to describe the variability within and between animals inplace of the single gamma distribution combining the two sources ofvariability.  When using two compound gamma distributions forpre-treatment data, the post-treatment data are paired to thepre-treatment data by animal.  The advantage of using a singledistribution is that the model is more identifiable, and therefore islikely to converge more quickly and return errors less frequently.  Theadvantage of the paired model is that it allows repeat measurements tobe incorporated in the model, and non-randomly missing data (ie.protocols that involve post-treatment sampling of only animals with ahigh pre-treatment count) to be modelled appropriately.  The simplemodel cannot be used when using repeat samples within animal, and willprovide inaccurate results if animals are targeted for post-treatmentsampling based on their pre-treatment count.  Default uses the simplemodel (FALSE).

#' @param zero.inflation option to use a zero-inflated gamma Poisson inplace of the gamma Poisson distribution.  If TRUE, zero inflateddistributions are used pre- and post-treatment (with the prevalencefixed between pre- and post-treatment distributions).  Default FALSE.

#' @param divide.data count division factor to allow egg count data in eggsper gram to be used raw (numeric).  Default 1 (no transformation todata).

#' @param record.chains option to allow the MCMC chains to be recorded forfuture use.  If TRUE, the function returns the MCMC object as part ofthe return value (the MCMC object is not printed using the printmethod).  If write.file==TRUE, the results are also saved using afilename containing the name of the analysis.  Default FALSE.

#' @param write.file option to write the results of the analysis to a textfile using a filename containing the name of the analysis.  The contentsof the text file are identical to the return value of the function. Default FALSE.

#' @param bootstrap.iters the number of bootstrap iterations to use for thebootstrap method.  Default 10000.

#' @param plot.graph an option to plot the posterior true egg countreduction from the MCMC method graphically.  If write.file==TRUE thenthe graph is saved as a PDF, otherwise the graph is plotted on theactive device.  Default FALSE.

#' @param skip.mcmc option to omit the MCMC analysis, and return bootstrapand WAAVP method analysis results alone.  Default FALSE.

#' @param ... other options to be passed directly to \code{\link[runjags]{autorun.jags}}.

reduction_model <- function(name = NA, pre.data = NA, post.data = NA, data = list(pre=pre.data, post=post.data), animal.names = FALSE, efficacy=95, confidence=95, restrict.efficacy = TRUE, control.animals = FALSE, paired.model = FALSE, fix.variation = FALSE, fix.efficacy = TRUE, zero.inflation = FALSE, individual.analysis = FALSE, divide.data = 1, record.chains = FALSE, write.file = FALSE, bootstrap.iters=10000, plot.graph = TRUE, skip.mcmc = NULL, stat.method = c('MCMC','WAAVP','Bootstrap'), custom.model = NA, ...){

	  if(!is.null(skip.mcmc)){
	    warning('Use of the skip.mcmc argument is deprecated')
	    if(skip.mcmc) stat.method <- stat.method[stat.method!='MCMC']
	  }
	  skip.mcmc <- FALSE
  
	  print('Check methods match what is possible.  Separate paired and individual MCMC methods?  Probably not')
  
	  # create new reduction.model function to just return the model like count.model does
  
	  # fecrt should be reduction.analysis with an alias for fecrt.analysis and alias for fecrt
  
	if(!paired.model & !fix.efficacy) stop('You must use the paired model if fix.efficacy is FALSE')
	### Unfixed efficacy currently won't work with Txcont stuff:
	if(sum(control.animals)>0 & !fix.efficacy) stop('Control animals cannot be used if fix.efficacy is FALSE')
	if(individual.analysis & !(paired.model & !fix.efficacy)) stop('You must use the paired model with fix.efficacy = FALSE if individual.analysis = TRUE')

	runname <- name

	# Individual analysis is removed for now, unless I add a distribution of efficacy for the paired model later.
	# individual.analysis=paired.model, 
	#if(!paired.model & individual.analysis) stop("Individual analysis is only available using the paired model")

	lci <- 0+((1-(confidence/100))/2)
	uci <- 1-((1-(confidence/100))/2)

	div <- divide.data 

	passthrough <- list(...)
	if(is.null(passthrough$max.time)) passthrough$max.time <- "1hr"
	if(is.null(passthrough$interactive)) passthrough$interactive <- FALSE
	if(is.null(passthrough$plots)) passthrough$plots <- FALSE

	arguments <- formals(autorun.jags)
	newargs <- formals(add.summary)
	arguments <- c(arguments, newargs[!names(newargs) %in% names(arguments)])
	arguments <- arguments[! names(arguments) %in% c('...', 'runjags.object')]

	for(i in 1:length(passthrough)){
		if(is.null(arguments[[names(passthrough)[i]]])){
			warning(paste("'", names(passthrough)[i], "' is not an argument to autorun.jags and was ignored"))
		}else{
			arguments[[names(passthrough)[i]]] <- passthrough[[i]]
		}
	}

	jags <- eval(arguments$jags)
	silent.jags <- eval(arguments$silent.jags)

	test.jags <- testjags(jags, silent=TRUE)
	if(test.jags[[2]][1]==FALSE){
		cat("Unable to call JAGS using '", jags, "'\n", sep="")
		stop("Unable to call JAGS")
	}

	testwritable <- new_unique("test")
	if(testwritable=="Directory not writable"){
		cat("\nThe working directory is not writable.  Please change the working directory\n\n")
		stop("Directory not writable")
	}

	if(class(data)=="function") stop("The class of the object supplied for data is a function")

	if(class(data)!="character" & class(data)!="matrix"){
		
		if(class(data)!="list") stop("Data supplied in an incorrect format")
		if(length(data)!=2) stop("Data supplied in an incorrect format - the list must be of length two")
	
		pre.data <- data[[1]]
		post.data <- data[[2]]
	
		if(class(pre.data)=="integer") pre.data <- as.numeric(pre.data)
		if(class(post.data)=="integer") post.data <- as.numeric(post.data)
	
		if(class(pre.data)!=class(post.data)) stop("The pre and post treatment data must be provided in the same format")
	
		if(class(pre.data)=="array" & paired.model==FALSE) stop("The paired model is required when using repeat samples within animal.  Either specify the data as a matrix, or set paired.model=TRUE")

		if(class(pre.data)=="matrix"){
			dims <- dim(pre.data)
			pre.data <- array(t(pre.data), dim=c(dims[1],dims[2],1))
			dims <- dim(post.data)
			post.data <- array(t(post.data), dim=c(dims[1],dims[2],1))
		}
	
		if(class(pre.data)=="numeric" | class(pre.data)=="integer"){
			pre.data <- array(pre.data, dim=c(length(pre.data),1,1))
			post.data <- array(post.data, dim=c(length(post.data),1,1))
		}
	
		if(dim(pre.data)[1] != dim(post.data)[1]) stop("Unequal numbers of animals pre- and post-treatment")
	
		data <- list(pre.counts=pre.data, post.counts=post.data)
	
	}



	cat(" ", strwrap("--- FECRT: Analyse faecal egg cout reduction test data using a Bayesian distributional simulation model ---"), "", sep="\n")
	cat(strwrap("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*"), sep="\n")
	cat(strwrap("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*"), "", sep="\n")

	if(is.na(runname)){
		runname <- ask(prompt = "Please enter a name for this analysis:  ", type="character")
	}

	datain <- data
	datana <- data
	datana <- as.vector(na.omit(datana))
	length(datana) <- 1

	dataok=FALSE

	if(class(data)=="matrix" | class(data)=="list"){
		dataok <- TRUE
	}else{
		if(is.na(datana)==FALSE){
		exists <- try(file.exists(datana), silent=TRUE)
		if((class(exists)=="try-error")==FALSE){
			if(exists==TRUE){
				suppressWarnings(data <- try(read.csv(datain, header=FALSE), silent=TRUE))
			}
		}
		suppressWarnings(valid.data <- try((length(as.matrix(data[,1])) > 1), silent=TRUE))
		if((class(valid.data)=="try-error")==TRUE){
			cat("ERROR:  The path you have entered does not appear to be valid\n")
		}else{
			if(valid.data==FALSE){
				cat("ERROR:  Invalid path / data\n") 
			}else{
				dataok=TRUE
			}
		}
		cat("\n")
		}
	}
	while(dataok==FALSE){
		datain <- ask(prompt = "Please enter the path to a (comma delimited) CSV file containing the data (type 'exit' to quit):  ", type="character")
		if((datain=="exit")==TRUE){
			stop("User exited the program")
		}
		exists <- try(file.exists(datain), silent=TRUE)
		if((class(exists)=="try-error")==FALSE){
			if(exists==TRUE){
			data <- try(as.matrix(read.csv(datain, header=FALSE), silent=TRUE))
			if((class(data)=="try-error")==FALSE){
				valid.data <- try(length(data[,1]) > 1, silent=TRUE)
				if((class(valid.data)=="try-error")==TRUE){
					cat("ERROR:  The path you have entered does not appear to be valid\n")
				}else{
					if(valid.data==FALSE){
						cat("ERROR:  The data you have entered is of length less than 2\n")
					}else{
						dataok=TRUE
					}
				}
			}else{
				cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
			}
			}else{
				cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
			}
		}else{
			cat("ERROR:  The path you entered does not appear to be valid; please retry (path may be absolute or relative)\n")
		}
		cat("\n")

	}

	if(class(data)!="matrix" & identical(animal.names,TRUE)){
		warning("'animal.names' cannot be TRUE if the data is not provided as a matrix.  'animal.names' will be set to FALSE")
		animal.names <- FALSE
	}

	if(class(data)=="matrix"){
	
		if(ncol(data) < (2+identical(animal.names,TRUE)) | ncol(data) > (3+identical(animal.names,TRUE))) stop("If provided as a matrix, the data must have either 2, 3 or 4 columns (if animal.names==TRUE)")
	
		if(identical(animal.names,TRUE)){
			animal.names <- data[,1]
			data <- data[,2:ncol(data)]
		}
		if(ncol(data)==3) control.animals <- data[,3]

		pre.data <- array(data[,1], dim=c(length(data[,1]),1,1))
		post.data <- array(data[,2], dim=c(length(data[,2]),1,1))

		data <- list(pre.counts=pre.data, post.counts=post.data)
	}

	if(class(data$pre.counts) == "NULL" | class(data$post.counts) == "NULL") stop("An error occured while transforming the data to an appropriate format")

	if(!identical(animal.names, FALSE)){
		if(length(animal.names) != dim(data$pre.counts)[1]) stop("The length of the character vector animal.names does not match the length of the data provided")
		names <- animal.names
	}else{
		names <- paste("Animal ", 1:dim(data$pre.counts)[1], sep="")
	}

	N <- dim(data$pre.counts)[1]

	if(N > 100 & !skip.mcmc){
		if(arguments$interactive){
			if(individual.analysis){
				cat("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  ")
				returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may cause a crash)\n2 Set individual.analysis to FALSE\n3 Skip the MCMC analysis\n:", "numeric", c(1,3))
				cat("\n")
				if(returned==2){
					individual.analysis <- FALSE
				}
				if(returned==3) skip.mcmc <- TRUE
			}#else{
			#	cat("There are a large number of animals in the data which will take a long time to analyse using MCMC.  ")
			#	returned <- ask("Please choose from the following options:\n1 Continue with analysis (this may take a while)\n2 Skip the MCMC analysis\n:", "numeric", c(1,2))
			#	cat("\n")
			#	if(returned==2) skip.mcmc <- TRUE
			#}
		}else{
			if(individual.analysis){
				individual.analysis <- FALSE
				warning("There are a large number of animals in the data which will take a long time to analyse using MCMC, and may cause memory issues with individual analysis.  individual analysis was therefore not performed.")
			}
		}
	}

	pre <- data$pre.counts / div
	post <- data$post.counts / div

	if(identical(control.animals, FALSE)){
		control.animals <- replicate(N, FALSE)
	}

	control.animals <- as.integer(control.animals)

	if(length(control.animals)!=N) stop("The length of the control/treatment vector does not match the number of animals")

	if(length(pre)!=length(post)) warning("Pre and post treatment data are of different lengths")


	if(sum(post,na.rm=TRUE)<5 & !fix.variation) warning('The difficulty in estimating the variability of the small number of counts in the post treatment data may make the model unstable.  You should try running the model again with fix.variation = TRUE.')

	#### Not using controls as possibility to boost inference about pre.mean but required if all treatment animals:
	usecontrol <- sum(control.animals)>0

	print('totally changed ordering of data arrays - multiple formats of data entry needs serious testing and help files updating')

	pairedmodel <- paste("model{

		for(animal in 1:N){
		for(sample in 1:Pre.Samples[animal]){
			for(chamber in 1:Pre.Chambers[animal]){
				Pre[animal,sample,chamber] ~ dpois(xpre.lambda[animal,sample])
			}
			xpre.lambda[animal,sample] <- ", if(zero.inflation) "probpos[animal] * ", "ind.pre.mean[animal] * pre.gamma[animal,sample]
			pre.gamma[animal,sample] ~ dgamma(pre.disp, pre.disp)T(10^-200,)
		}
		for(sample in 1:Post.Samples[animal]){
			for(chamber in 1:Post.Chambers[animal]){
				Post[animal,sample,chamber] ~ dpois(xpost.lambda[animal,sample])
			}
			xpost.lambda[animal,sample] <- ", if(zero.inflation) "probpos[animal] * ", "ind.post.mean[animal] * post.gamma[animal,sample]
			post.gamma[animal,sample] ~ dgamma(post.disp[Txcont[animal]], post.disp[Txcont[animal]])T(10^-200,)
		}	

			ind.pre.mean[animal] <- pre.mean * animal.gamma[animal]

			animal.gamma[animal] ~ dgamma(animal.disp, animal.disp)T(10^-200,)

			ind.post.mean[animal] <- ind.pre.mean[animal] * ind.delta.mean[animal]

			", if(zero.inflation) "probpos[animal] ~ dbern(prob)", "

			", if(fix.efficacy) "ind.delta.mean[animal] <- delta.m[Txcont[animal]] " else "ind.delta.mean[animal] ~ dgamma(alpha, beta)T(10^-200,)#dbeta(alpha, beta)", "

		}

		", if(!fix.efficacy) "
		#alpha ~ dunif(0,1000)
		#beta <- (alpha - (delta.mean*alpha)) / delta.mean
		
		alpha ~ dgamma(0.01, 0.01)
		beta <- alpha/delta.mean
	
		sample.delta.m <- mean(ind.delta.mean[])
		", "

		# these represent within animal variability so can change with tx:
		pre.disp <- 1 / ia
		post.disp[1] <- pre.disp * delta.disp
		post.disp[2] <- pre.disp
		animal.disp <- 1 / iaa

		# Priors
		pre.mean ~ dgamma(0.001, 0.001) #dunif(0.0001,10000) #
		", if(zero.inflation) "prob ~ dunif(0,1)", "

		ia <- exp(logia)
		logia ~ dunif(-9.21,4.6)
		iaa <- exp(logiaa)
		logiaa ~ dunif(-9.21,4.6)

		ddispl <- ", {if(TRUE) exp(-9.21) else 0.02}, " / pre.disp;
		ddispu <- ", {if(TRUE) exp(4.6) else 100}, " / pre.disp;

		", if(usecontrol) "delta.mean <- delta.m[1] / delta.m[2]" else "delta.mean <- delta.m[1]", "

		delta.m[1] ~ ", {if(restrict.efficacy) "dbeta(1,1)" else "dunif(0, 10)"}, "
		delta.m[2] <- exp(ldelta.m)
		ldelta.m ~ dunif(-4,4)

		delta.disp ", {if(fix.variation) "<- 1" else "~ dlnorm(0, 0.01)T(0.001, 1000)"}, "

	}", sep="")

		print('swapped sample and repeat')


	singlemodel <- paste("model{

	for(animal in 1:N){
	for(sample in 1:Pre.Samples[animal]){
		for(chamber in 1:Pre.Chambers[animal]){
			Pre[animal,sample,chamber] ~ dpois(xpre.lambda[animal,sample])
		}
		xpre.lambda[animal,sample] <- ", if(zero.inflation) "probpos[animal] * ", "pre.mean * pre.gamma[animal,sample]
		pre.gamma[animal,sample] ~ dgamma(pre.disp, pre.disp)T(10^-200,)
	}
	for(sample in 1:Post.Samples[animal]){
		for(chamber in 1:Post.Chambers[animal]){
			Post[animal,sample,chamber] ~ dpois(xpost.lambda[animal,sample])
		}
		xpost.lambda[animal,sample] <- ", if(zero.inflation) "probpos[animal] * ", "post.mean[Txcont[animal]] * post.gamma[animal,sample]
		post.gamma[animal,sample] ~ dgamma(post.disp[Txcont[animal]], post.disp[Txcont[animal]])T(10^-200,)
	}	

	", if(zero.inflation) "probpos[animal] ~ dbern(prob)","
		}

		pre.disp <- 1 / ia

		post.mean[1] <- pre.mean * delta.m[1]
		post.mean[2] <- ", if(usecontrol) "pre.mean * delta.m[2]" else "pre.mean", "
		post.disp[1] <- pre.disp * delta.disp
		post.disp[2] <- pre.disp

		", if(usecontrol) "delta.mean <- delta.m[1] / delta.m[2]" else "delta.mean <- delta.m[1]", "

		# Priors
		pre.mean ~ dgamma(0.001, 0.001) #dunif(0.0001,10000) #
	", if(zero.inflation) "prob ~ dunif(0,1)","
		ia <- exp(logia)
		logia ~ dunif(-9.21,4.6)#4.6

		delta.m[1] ~ ", {if(restrict.efficacy) "dbeta(1,1)" else "dunif(0, 10)"}, "
		delta.m[2] <- exp(ldelta.m)
		ldelta.m ~ dunif(-4,4)
		delta.disp ", {if(fix.variation) "<- 1" else "~ dlnorm(0, 0.01)T(0.001, 1000)"}, "

		}", sep="")


	print('changed')
	probposinit <- as.integer(apply(pre, 1, sum) == 0)
	probposinit[is.na(probposinit)] <- 1
	probinit <- max(sum(probposinit) / length(probposinit), 0.1)

	preg <- matrix(1, nrow=N, ncol=dim(pre)[2])
	postg <- matrix(1, nrow=N, ncol=dim(post)[2])

	scpre <- pre
	scpost <- post
	for(i in 1:N){
		if(all(is.na(scpre[i,,]))) scpre[i,1,1] <- 0
		if(all(is.na(scpost[i,,]))) scpost[i,1,1] <- 0
	}

	pre.samples <- apply(scpre, 1, function(z) return(max(apply(z, 2, function(x) return(max(which(!is.na(x))))))))
	Pre.Chambers <- apply(scpre, 1, function(z) return(max(apply(z, 1, function(x) return(suppressWarnings(max(which(!is.na(x)))))))))
	post.samples <- apply(scpost, 1, function(z) return(max(apply(z, 2, function(x) return(max(which(!is.na(x))))))))
	Post.Chambers <- apply(scpost, 1, function(z) return(max(apply(z, 1, function(x) return(suppressWarnings(max(which(!is.na(x)))))))))

	for(i in 1:length(pre.samples)){
		if(pre.samples[i] < ncol(preg)){
			preg[i,(pre.samples[i]+1):ncol(preg)] <- NA
		}
	}
	for(i in 1:length(post.samples)){	
		if(post.samples[i] < ncol(postg)){
			postg[i,(post.samples[i]+1):ncol(postg)] <- NA
		}
	}

	preg <- preg[,1:max(pre.samples),drop=FALSE]
	postg <- postg[,1:max(post.samples),drop=FALSE]

	msg <- paste("Assessing the faecal egg count reduction for ", N, " animals", sep="")
	if(sum(control.animals)>0) msg <- paste(msg, " (including ", sum(control.animals==1), " control animals)", sep="")
	msg <- paste(msg, ".  There are ", {if(all(pre.samples==pre.samples[1])) pre.samples[1] else paste(min(pre.samples), "-", max(pre.samples), sep="")}, " pre samples per animal (with ", {if(all(Pre.Chambers==Pre.Chambers[1])) Pre.Chambers[1] else paste(min(Pre.Chambers), "-", max(Pre.Chambers), sep="")}, " McMasters chambers), and ", {if(all(post.samples==post.samples[1])) post.samples[1] else paste(min(post.samples), "-", max(post.samples), sep="")}, " post samples per animal (with ", {if(all(Post.Chambers==Post.Chambers[1])) Post.Chambers[1] else paste(min(Post.Chambers), "-", max(Post.Chambers), sep="")}, " McMasters chambers).  This will take some time...\n", sep="") 
	#browser()
	cat(strwrap(msg), sep='\n')

	#cat('We have ', N, ' animals, ', control.animals, ' controls, {', pre.samples, '} pre samples, {', post.samples, '} post samples, {', Pre.Chambers, '} pre chambers, and {', Post.Chambers, '} post chambers\n', sep='')


	datastring <- dump.format(list(Pre=pre, Post=post, N=N, Pre.Samples = pre.samples, Pre.Chambers = Pre.Chambers, Post.Samples = post.samples, Post.Chambers = Post.Chambers, Txcont=control.animals+1))

	if(!paired.model){

		model <- singlemodel

		inits1 <- list(probpos=probposinit, pre.mean = max(mean(pre, na.rm = TRUE) * 2, 1), delta.m=c(0.01,NA), ldelta.m=0, logia=log(0.1), pre.gamma=preg, post.gamma=postg, prob=probinit, alpha=1)

		inits2 <- list(probpos=replicate(length(pre), 1), pre.mean = max(mean(pre, na.rm = TRUE) * 2, 1), delta.m=c(0.1,NA), ldelta.m=0, logia=log(10), pre.gamma=preg, post.gamma=postg, prob=1, alpha=1)
	
		if(!fix.variation){
			inits1$delta.disp <- 2
			inits2$delta.disp <- 0.5
		}	
		
	}else{

		model <- pairedmodel
	
		idm <- replicate(N, 0.05)
	
		inits1 <- list(probpos=probposinit, ind.delta.mean=idm,  pre.mean = max(mean(pre, na.rm = TRUE) * 2, 1), animal.gamma=replicate(N,1), delta.m=c(0.01,NA), ldelta.m=0, logia=log(0.1), logiaa=log(0.1), pre.gamma=preg, post.gamma=postg, prob=probinit)

		inits2 <- list(probpos=replicate(length(pre), 1), ind.delta.mean=idm, pre.mean = max(mean(pre, na.rm = TRUE) * 2, 1), animal.gamma=replicate(N,1), delta.m=c(0.99,NA), ldelta.m=0, logia=log(10), logiaa=log(10), pre.gamma=preg, post.gamma=postg, prob=1)
	
		if(!fix.variation){
			inits1$delta.disp <- 2
			inits2$delta.disp <- 0.5
		}	
		
	}

	monitor=c("pre.mean", "delta.mean", "pre.disp")
	if(sum(control.animals)>0) monitor <- c(monitor, "delta.m")
	if(!fix.variation) monitor <- c(monitor, "delta.disp")
	if(paired.model) monitor <- c(monitor, "pre.disp", "animal.disp")
	#monitor <- c(monitor, "pre.disp", "animal.disp", "ind.pre.mean", "ind.post.mean")
	if(individual.analysis) monitor <- c(monitor, "ind.delta.mean")
	if(!fix.efficacy) monitor <- c(monitor, "alpha", "sample.delta.m")
	if(zero.inflation) monitor <- c(monitor, "prob")
	if(zero.inflation & individual.analysis) monitor <- c(monitor, "probpos")

	if(!skip.mcmc){
	
		arguments$model <- model
		arguments$inits <- c(dump.format(inits1), dump.format(inits2))
		arguments$n.chains <- length(arguments$inits)
		arguments$data <- datastring
		arguments$monitor <- monitor
		arguments$plots <- FALSE

		class(arguments) <- "list"
	
		results <- do.call(autorun.jags, arguments, quote=FALSE)
	
		#results <- autorun.jags(data=datastring, model=model, monitor=monitor, n.chains=2, inits=c(inits1, inits2), silent.jags = list(silent.jags=silent.jags, killautocorr=TRUE), plots = FALSE, thin.sample = TRUE, interactive=interactive, max.time=max.time, ...)

	if(results[1]=="Error"){
		#print(results)
		cat("An unexpected error occured during the simulation.  Ensure that the values of divide.data, pre.data and post.data provided are correct and re-run the functon.\n\n--- Simulation aborted ---\n\n")
		stop("An error occured during the simulation")
	}
	if(any(names(results)=="pilot.mcmc")){
		warning("The chains did not achieve convergence during the simulation, you should interpret the Bayesian results with extreme caution")
		converged <- FALSE
		results$mcmc <- results$pilot.mcmc
		results$summary <- results$pilot.summary
	}else{
		converged <- TRUE
	}

	}else{
		results <- list(mcmc="MCMC analysis not performed")
	}


	blankquant <- quantile(0, probs=c(lci, 0.5, uci))

	# No bootstrapping or WAAVP for paired model type data

	if(any(c(dim(post)[2:3], dim(pre)[2:3])>1)){
		cat("Bootstrap and WAAVP calculations are not available for repeated pre and/or post treatment egg counts\n")
		boot.reductions = bootquant = method.boot = waavpquant = method.waavp = bootprob <- NA
	}else{
	
		pre <- pre[,1,1]#apply(pre.data, 3, mean) # Should only be 1 datapoint
		post <- post[,1,1]#apply(post.data, 3, mean) # Should only be 1 datapoint
		# Just to check:
		#pre2 <- pre.data[1,1,]
		#post2 <- post.data[1,1,]
		#if(any(pre!=pre2) | any(post!=post2)) stop("An unexpected error occured while manipulating the data for the bootstrap/WAAVP methods")
	
		if(sum(control.animals)>0){
			warning("Control animals are ignored using the bootstrap and WAAVP calculations")
			pre <- pre[!control.animals]
			post <- post[!control.animals]
		}
	
		cat("Calculating the Bootstrap and WAAVP method analysis...\n")
		# Bootstrap:
	
		boot.pre <- matrix(data=sample(pre, length(pre)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(pre))
		boot.post <- matrix(data=sample(post, length(post)*bootstrap.iters, replace=TRUE), ncol=bootstrap.iters, nrow=length(post))
	
		boot.reductions <- (1 - apply(boot.post, 2, mean) / apply(boot.pre, 2, mean)) *100

		boot.reductions[boot.reductions==-Inf] <- NA
		bootprob <- sum(boot.reductions < efficacy, na.rm=TRUE) / sum(!is.na(boot.reductions)) * 100
	
		bootquant <- quantile(boot.reductions, probs=c(lci, 0.5, uci), na.rm=TRUE)
		if(is.na(bootprob)){
			method.boot <- "Returned error"
		}else{
			method.boot <- if(bootprob >= (uci*100)) "Confirmed resistant" else if(bootprob >= 50) "Probable resistant" else if(bootprob >= (lci*100)) "Possible resistant" else "Confirmed susceptible"
		}

		# WAAVP:
	
		waavpquant <- blankquant
	
		pre.data <- pre
		post.data <- post

		N <- length(pre.data)
		arith.pre <- mean(pre.data)
		arith.post <- mean(post.data)
		var.pre <- var(pre.data)
		var.post <- var(post.data)

		red <- 100*(1-arith.post/arith.pre)
		var.red <- var.post/(N*arith.post^2) + var.pre/(N*arith.pre^2)
		try(if(var.post==0 & arith.post==0) var.red <- var.pre/(N*arith.pre^2))

		upper.ci <- 100 * (1-(arith.post/arith.pre * exp(-2.048*sqrt(var.red))))
		lower.ci <- 100 * (1-(arith.post/arith.pre * exp(2.048*sqrt(var.red))))

		waavpquant[1] <- lower.ci
		waavpquant[2] <- red
		waavpquant[3] <- upper.ci
	
		if(any(is.na(c(red,efficacy,lower.ci)))){
			method.waavp <- "Returned error"
		}else{
		method.waavp <- "susceptible"
		if(red < 95 | lower.ci < 90){
			method.waavp <- "Suspected resistant"
			if(red < 95 & lower.ci < 90){
				method.waavp <- "Confirmed resistant"
			}
		}
	
		if(confidence!=95){
			waavpquant[] <- NA
			method.waavp <- "Not available"
		}
	}
	}

	if(!skip.mcmc){

	cat("Calculating the Bayesian method analysis...\n")

	reduction <- (1-unlist(results$mcmc[,"delta.mean"]))*100
	pre.mean <- (unlist(results$mcmc[,"pre.mean"]))
	if(!fix.variation) deltashape <- unlist(results$mcmc[,"delta.disp"]) else deltashape <- NA
	if(zero.inflation) zi <- (1-unlist(results$mcmc[,"prob"]))*100 else zi <- NA

	vars <- nvar(results$mcmc)

	n <- (vars-(length(monitor)-(1+zero.inflation)))/(1+zero.inflation)

	medred <- median(reduction)
	ci <- HPDinterval(as.mcmc(reduction), prob=(0.01*confidence))
	l95red <- ci[1]
	u95red <- ci[2]
	mcmcquant <- blankquant
	mcmcquant[1] <- l95red
	mcmcquant[3] <- u95red
	mcmcquant[2] <- medred

	meanquant <- blankquant
	meanquant[2] <- median(pre.mean*divide.data)
	ci <- HPDinterval(as.mcmc(pre.mean*divide.data), prob=(0.01*confidence))
	meanquant[1] <- ci[1]
	meanquant[3] <- ci[2]

	ziquant <- blankquant
	ziquant[2] <- if(zero.inflation) median(zi) else NA
	ci <- if(zero.inflation) HPDinterval(as.mcmc(zi), prob=(0.01*confidence)) else c(NA,NA)
	ziquant[1] <- ci[1]
	ziquant[3] <- ci[2]

	# post.disp = pre.disp * delta.disp
	# 1/pod^2 = 1/prd^2 * dd
	# 1/pod^2 = 1/prd^2 * dd
	# prd^2 = pod^2 * dd

	dshapequant <- blankquant
	if(fix.variation){
		dshapequant[] <- NA
	}else{
		dshapequant[2] <- median(deltashape)
		ci <- HPDinterval(as.mcmc(deltashape), prob=(0.01*confidence))
		dshapequant[1] <- ci[1]
		dshapequant[3] <- ci[2]
	}

	mcmcprob <- sum(reduction < efficacy) / length(reduction) * 100

	method.mcmc <- if(mcmcprob >= (uci*100)) "Confirmed resistant" else if(mcmcprob >= 50) "Probable resistant" else if(mcmcprob >= (lci*100)) "Possible resistant" else "Confirmed susceptible"

	#class <- (prob.res > 0.025) + (prob.res > 0.50) + (prob.res > 0.975) + 1
	#method3 <- switch(class, "1"="susceptible", "2"="possible", "3"="probable", "4"="resistant")

	}else{
		method.mcmc = mcmcquant = mcmcprob = ziquant = dshapequant = meanquant = indredquant = ind.prob.inf = ind.zi.probs = converged = efficacy = restrict.efficacy <- NA
	}

	cat("Finished calculations\n")

	s <- try(output <- list(results.mcmc=method.mcmc, quant.mcmc=mcmcquant, results.boot=method.boot, quant.boot=bootquant, results.waavp=method.waavp, quant.waavp = waavpquant, prob.mcmc=mcmcprob, prob.boot=bootprob, ziquant=ziquant, dshapequant=dshapequant, meanquant=meanquant, confidence=confidence, converged=converged, efficacy=efficacy, restrict.efficacy=restrict.efficacy, animal.names=names, name=name, model=model))
	if(class(s)=='try-error') browser()#print('remove browser()')
	if(record.chains) output <- c(output, list(runjags=results))

	#browser()

	if(write.file){
		cat("Writing results to file\n")
		filename <- new_unique(paste(name, ".results",sep=""), ".txt", ask=arguments$interactive)
		print.fecrt.results(output, filename=filename)
		filename <- new_unique(paste(name, ".graph",sep=""), ".pdf", ask=arguments$interactive)
		if(!skip.mcmc){
			pdf(file=filename, width=10, height=6)
			par(mfrow=c(1,2))
			hist(reduction, col="red", breaks="fd", main="Posterior distribution", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100), border="red")
			abline(v=efficacy, lty="dashed")
			plot.ecdf(reduction, ylab="Cumulative probability", main="Cumulative distribution", xlab="True FEC reduction (%)", col="red")
			abline(v=efficacy, lty="dashed")
			abline(h=output$prob.mcmc/100, lty="dashed")		
			dev.off()
		}
		if(record.chains) save(results, file=new_unique(paste("fecrt.", name, sep=""), ".Rsave", ask=FALSE))
		cat("Results file written successfully\n")
	}else{
		if(plot.graph){
			if(!skip.mcmc){
				hist(reduction, col="red", breaks="fd", main="Posterior distribution for the true reduction", freq=FALSE, ylab="Probability", xlab="True FEC reduction (%)", xlim=c(0, 100))
				abline(v=efficacy)
			}
		}
	}



	cat("\nAnalysis complete\n\n", sep="")
	cat(strwrap("*PLEASE NOTE:  THIS SOFTWARE IS INTENDED FOR EDUCATIONAL PURPOSES ONLY AND SHOULD NOT BE RELIED UPON FOR REAL WORLD APPLICATIONS*"), sep="\n")
	cat(strwrap("*ANALYSING DATA USING MCMC SAMPLING CAN PRODUCE MISLEADING RESULTS IF USED INAPPROPRIATELY*"), "", sep="\n")
	cat("--- End ---\n\n")

	class(output) <- "fecrt_results"
	return(output)

}

reduction.model <- reduction_model
fecrt.model <- reduction_model
FECRT.model <- reduction_model
