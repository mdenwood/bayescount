#' @name lnormal_params
#' @aliases normal_params lnormal.params normal.params
#' @title Calculate the Log-Normal Mean and Standard Deviation Using the Normal Mean and Standard Deviation and Vice Versa
#'
#' @description
#' Function to calculate the equivalent values for the mean and standard deviation of a normal distribution from the mean and standard deviation of the log-normal distribution and vice versa.  Outputs from this function can be used with the dnorm() or dlnorm functions, and with the equivalent distributions in JAGS.
#'
#' @seealso \code{\link{bayescount}} \code{\link[stats]{Lognormal} \code{\link[stats]{Normal}
#'
#' @keywords methods
#'
#' @examples
#' lmean <- 2.5
#' lsd <- 0.2
#' mean <- normal.params(lmean,lsd)[[1]]
#' sd <- normal.params(lmean,lsd)[[2]]
#' 
#' curve(dlnorm(x, lmean, lsd), from=0, to=25)
#' dev.new()
#' curve(dnorm(x, mean, sd), from=0, to=25)
#' 
#' mean <- 10
#' sd <- 2
#' lmean <- lnormal.params(mean,sd)[[1]]
#' lsd <- lnormal.params(mean,sd)[[2]]
#' 
#' curve(dnorm(x, mean, sd), from=0, to=20)
#' dev.new()
#' curve(dlnorm(x, lmean, lsd), from=0, to=20)



normal_params <- function(log.mean, log.sd, coeff.variation=sqrt(exp(log.sd^2)-1)){

dispersion <- coeff.variation
log.sd <- sqrt(log(dispersion^2 +1))

lsd <- log.sd
lmu <- log.mean

mu <- exp(lmu + ((lsd^2) / 2))
sd <- sqrt(exp((2*lmu)+(lsd^2)) * (exp(lsd^2) - 1))
dispersion <- sd / mu

return(list(mean=mu, sd=sd, coeff.variation=dispersion))
}

lnormal_params <- function(mean, sd, coeff.variation=sd/mean){

tau <- coeff.variation
lsd <- sqrt(log(tau^2 +1))
lmu <- log(mean) - ((lsd^2) / 2)
dispersion <- tau

return(list(log.mean=lmu, log.sd=lsd, coeff.variation=dispersion))
}

normal.params <- normal_params
lnormal.params <- lnormal_params