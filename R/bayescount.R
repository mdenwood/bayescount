#' @name bayescount
#' @aliases bayescountpackage bayescount-package
#' @docType package
#' @title Analysis and power calculations for faecal egg count (FEC) and faecal egg count reduction test (FECRT) data using computationally intensive statistical methods
#'
#' @description
#' The bayescount package analyses the over-dispersed count datasets(such as those commonly encountered in parasitology) usingcomputationally intensive statistical techniques such asnon-parametric bootstrapping and Bayesian MCMC.  The FEC data analysisreturns information on the mean count, coefficient of variation andzero inflation (true prevalence) present in the data.  The faecal eggcount reduction test (FECRT) analysis calculates efficacy as well aspost-treatment changes in count variation, and allows the analysis ofrepeat observations from the same animals.  Power and sample sizecalculations for FECRT and FEC studies (given any expected parameter values) are also provided.  Requires Just Another Gibbs Sampler (JAGS)for MCMC analysis functions (see http://mcmc-jags.sourceforge.net)
#'
#' \emph{FEC Analysis}
#'
#' You can do that
#' 
#' \emph{FECRT analysis}
#' 
#' That too
#' 
#' \emph{Power}
#' 
#' And that
#' 
#' @references
#' M. J. Denwood, S. W. J. Reid, S. Love, M. K. Nielsen, L. Matthews, I. J. McKendrick, and G. T. Innocent. Comparison of three alternative methods for analysis of equine Faecal Egg Count Reduction Test data. Prev. Vet. Med. (2009), doi:10.1016/j.prevetmed.2009.11.009
#' 
#' Denwood, M. J. (2010). A quantitative approach to improving the analysis of faecal worm egg count data. University of Glasgow. Retrieved from http://www.gla.ac.uk/media/media_149338_en.pdf
#' 
#' @seealso \code{\link{fecrt}}, \code{\link{bayescount}}, \code{\link{bayescount.single}}, \code{\link{fec.power}}, \code{\link{fecrt.power}}, \code{\link[runjags]{runjags}}

NULL
