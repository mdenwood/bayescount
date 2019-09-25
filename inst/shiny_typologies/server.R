# Should not be called in this code, but bayescount needs to have been
# installed this way before deployment:
# devtools::install_github('mdenwood/bayescount')

library('shiny')
library('bayescount')
library('dplyr')
library('tibble')
library('tidyr')
library('ggplot2')

options(stringsAsFactors=FALSE)

# Permanent settings:
citers <- 10^3

pubdate <- "2019-09-24"



frameworks <- list(
	Denwood = paste0('<p>Under the Denwood et al classification system, typologies are classified as:</p><p>Typologies 1a/1b/1c:  Reduced Efficacy</p><p>Typologies 2a/2b/2c:  Inconclusive</p><p>Typology 3:  Borderline Efficacy</p><p>Typologies 4a/4b/4c:  Adequate Efficacy</p>'),
	Kaplan = paste0('<p>Under the Kaplan et al classification system, typologies are classified as:</p><p>Typologies 1a/1b:  Resistant</p><p>Typologies 1c/3:  Low Resistant</p><p>Typologies 2a/2b/2c:  Inconclusive</p><p>Typologies 4a/4b/4c:  Susceptible</p>'),
	Coles = paste0('<p>Under the (modified) Coles et al classification system, typologies are classified as:</p><p>Typologies 1a/1b/2a/2b:  Resistant</p><p>Typologies 1c/2c/3/4a:  Suspected Resistant (here re-classified as inconclusive)</p><p>Typologies 4b/4c:  Susceptible</p>')
)

classifications <- list(
	Denwood = list(reduced = c('1a','1b','1ab','1c'), inconclusive = c('2a','2b','2c','SumPost0'), borderline = '3', adequate = c('4a','4b','4c','4bc')),
	Kaplan = list(resistant = c('1a','1b','1ab'), inconclusive = c('2a','2b','2c','SumPost0'), lowresistant = c('1c','3'), susceptible = c('4a','4b','4c','4bc')),
	Coles = list(resistant = c('1a','1b','1ab','2a','2b'), inconclusive = c('1c','2c','3','4a','SumPost0'), susceptible = c('4b','4c','4bc'))
	)
	

# In-built k settings:
kparameters <- tribble(~Species, ~k1, ~k2, ~cor,
			"Sheep", 1.2, 0.7, 0.5,
			"Cattle", 1.1, 0.5, 0.1,
			"Calves", 1.5, 0.4, 0.1,
			"Equine", 0.6, 0.6, 0.5
			)
## TODO: build this from kparameters directly
kstring <- list(
	Sheep = "<p>Pre-set over-dispersion parameters for Sheep are:</p><p>k1: 1.2</p><p>k2: 0.7</p><p>correlation: 0.5</p>",
	Cattle = "<p>Pre-set over-dispersion parameters for Cattle are:</p><p>k1: 1.1</p><p>k2: 0.5</p><p>correlation: 0.1</p>",
	Calves = "<p>Pre-set over-dispersion parameters for Calves are:</p><p>k1: 1.5</p><p>k2: 0.4</p><p>correlation: 0.1</p>",
	Equine = "<p>Pre-set over-dispersion parameters for Horses are:</p><p>k1: 0.6</p><p>k2: 0.6</p><p>correlation: 0.5</p>",
	Custom = "<p>Enter custom over-dispersion parameters using the sliders on the left</p>"
	)

getsimres <- function(simpars){

	if(simpars$kpreset!="Custom"){
		simpars$k1 <- NULL
		simpars$k2 <- NULL
		simpars$cor <- NULL
		simpars <- c(simpars, as.list(kparameters %>% filter(Species==simpars$kpreset) %>% select(-Species)))
	}
	
	tsp <- as.data.frame(simpars) %>%
		mutate(ti=target, ta=(target-delta), ko1=k1, ko2=k2, kd=k2/k1) %>%
		mutate(kc = sqrt(ko1*ko2) / cor, ku1 = ((kc+1)*ko1)/(kc-ko1), ku2 = ((kc+1)*ko2)/(kc-ko2))

	# Get lower and upper as function of ti and ta:
	upper <- 1 - max(0, tsp$ta-0.05)
	lower <- 1 - min(1, tsp$ti+0.05)

	# Make by a sensible number to give around 25 values depending on upper-lower:
	by <- seq(0,upper-lower,length.out=25)[2]
	choices <- c(0.0001, 0.00025, 0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05)
	by <- choices[which.min(abs(by-choices))]

	rvs <- seq(lower,upper,by=by)
	if(0 %in% rvs)
		rvs <- unique(c(0, choices[1:6], rvs))
	rvs <- unique(c(rvs, (1-tsp$ta)+c(-by/2,0,by/2), (1-tsp$ti)+c(-by/2,0,by/2)))

	args <- as.list((tsp %>% select(N, mu, k1, k2, kc, ta, ti)))
	args <- c(args, list(rvs=rvs, iters=citers, approx=1, priors=c(0,0), useml=TRUE, plot=FALSE))

	res <- do.call('powersim_paired', args)
	
	res$Method[res$Method=='Levecke'] <- 'Gamma'
	res$Method[res$Method=='MLE'] <- 'Asymptotic'
	
	return(list(simpars=simpars, simres=res))
}

inits <- getsimres(list(StudyType='Simple Paired', N=10, mu=20, target=99/100, delta=4/100, pthresh=0.025, kpreset='Sheep'))


function(input, output, session) {
	
	rv <- reactiveValues(simpars=inits$simpars, storedres=inits$simres)
	
	observeEvent(input$simulate, {
		
		simpars <- list(StudyType=input$type, N=input$N, mu=input$mu, target=input$target/100, delta=input$delta/100, pthresh=input$pthresh, kpreset=input$kpreset, k1=10^input$k1, k2=10^input$k2, cor=input$cor)
		
		newres <- getsimres(simpars)
		
		rv$storedres <- newres$simres
		rv$simpars <- newres$simpars
		
	})
	
	fluidPage(
		output$kvalues <- renderText(kstring[[input$kpreset]]),
		output$classifications <- renderText(paste0(frameworks[[input$statframe]], if(input$statmethod %in% c("Gamma","WAAVP","Asymptotic")) paste0("<p>NOTE: any datasets where sum(post)==0 will also be classified as inconclusive using the ", input$statmethod, " method as 95% CI are incalculable</p>")))
	)
 	
	output$typologies <- renderImage({
		list(src="typologies.png", contentType = 'image/png', width = 400, height = 600, alt = "Typologies image")
	}, deleteFile=FALSE)
	
	output$noninc_plot <- renderPlot({
		
		res <- rv$storedres
		
		res$Classification <- factor(res$Classification)
		levels(res$Classification) <- classifications[[input$statframe]]
		
		methodsusing <- input$statmethod
		
		typres <- res %>%
			filter(Method %in% methodsusing) %>%
			mutate(Efficacy = 1-Reduction) %>%
			group_by(Method, Efficacy, Classification) %>%
			summarise(Number = n()) %>%
			# Need to join with full dataset:
			full_join( crossing(Efficacy=unique(1-res$Reduction), Method= methodsusing, Classification=factor(levels(res$Classification), levels=levels(res$Classification))) , by = c("Method", "Efficacy", "Classification")) %>%
			# Any blank numbers are zeros:
			mutate(Number = ifelse(is.na(Number), 0, Number)) %>%
			group_by(Method, Efficacy) %>%
			mutate(Proportion = Number / sum(Number)) %>%
			ungroup()

		theme_set(theme_light())
		
		ti <- rv$simpars$target
		ta <- ti - rv$simpars$delta
		incres <- typres %>%
			filter(Classification=='inconclusive') %>%
			mutate(region = ifelse(Efficacy <= ta, 'one', ifelse(Efficacy < ti, 'two', 'three')))
	
		ggplot(incres, aes(x=Efficacy, y=1-Proportion, group=region)) +
			geom_rect(xmin=ta, xmax=ti, ymin=0, ymax=1, fill='grey80', col='grey80', alpha=0.5) +
		#	geom_line(stat='smooth', method='loess', formula=y~x, se=FALSE, alpha=0.5, lwd=2) +
		#	geom_line(stat='smooth', method='lm', formula=y~poly(x,2), se=FALSE, alpha=0.5, lwd=2) +
			geom_line(alpha=0.5, lwd=2) +
			geom_point() +
			geom_vline(xintercept=ti, lty='dashed') +
			geom_vline(xintercept=ta, lty='dotted') +
			geom_hline(yintercept=0.8, lty='solid') +
			facet_grid(rows="Method") +
			theme(legend.pos='none') +
			scale_y_continuous(breaks=seq(0,1,by=0.2), limits=c(0,1))
	})
	
	output$typologies_plot <- renderPlot({

		res <- rv$storedres

		res$Classification <- factor(res$Classification)
		levels(res$Classification) <- classifications[[input$statframe]]

		methodsusing <- input$statmethod
		
		typres <- res %>%
			filter(Method %in% methodsusing) %>%
			mutate(Efficacy = 1-Reduction) %>%
			group_by(Method, Efficacy, Classification) %>%
			summarise(Number = n()) %>%
			# Need to join with full dataset:
			full_join( crossing(Efficacy=unique(1-res$Reduction), Method= methodsusing, Classification=factor(levels(res$Classification), levels=levels(res$Classification))) , by = c("Method", "Efficacy", "Classification")) %>%
			# Any blank numbers are zeros:
			mutate(Number = ifelse(is.na(Number), 0, Number)) %>%
			group_by(Method, Efficacy) %>%
			mutate(Proportion = Number / sum(Number)) %>%
			ungroup()

		theme_set(theme_light())
		
		ti <- rv$simpars$target
		ta <- ti - rv$simpars$delta
		ggplot(typres, aes(x=Efficacy, y=Proportion, col=Classification)) +
			geom_rect(xmin=ta, xmax=ti, ymin=0, ymax=1, fill='grey80', col='grey80', alpha=0.5) +
		#	geom_line(stat='smooth', method='lm', formula=y~poly(x,10), se=FALSE, alpha=0.5, lwd=2) +
		#	geom_line(stat='smooth', method='loess', formula=y~x, se=FALSE, alpha=0.5, lwd=2) +
			geom_line(alpha=0.5, lwd=2) +
			geom_point() +
			geom_vline(xintercept=ti, lty='dashed') +
			geom_vline(xintercept=ta, lty='dotted') +
			geom_hline(yintercept=0.025, lty='solid') +
			facet_grid(rows="Method") +
			theme(legend.pos='bottom')
		
	})
	
	output$parameters <- renderText(paste0('<p>Showing results for the ', input$statmethod, ' method and ', input$statframe, ' et al framework, with parameters: N = ', rv$simpars$N, ', mu = ', rv$simpars$mu, ', k1 = ', rv$simpars$k1, ', k2 = ', rv$simpars$k2, ', correlation = ', rv$simpars$cor, ', target = ', rv$simpars$target, ', delta = ', rv$simpars$delta, ', alpha = ', rv$simpars$pthresh))
	
	output$footer <- renderText(paste0('<p align="center">Power calculation tool for FECRT data by Matthew Denwood (last updated ', pubdate, ')<br><a href="http://www.fecrt.com/", target="_blank">Learn more (opens in a new window)</a></p>'))
}
