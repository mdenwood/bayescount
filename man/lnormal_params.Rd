% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normal_lnormal_params.R
\name{lnormal_params}
\alias{lnormal_params}
\alias{normal_params}
\alias{lnormal.params}
\alias{normal.params}
\title{Calculate the Log-Normal Mean and Standard Deviation Using the Normal Mean and Standard Deviation and Vice Versa}
\usage{
normal_params(log.mean, log.sd, coeff.variation = sqrt(exp(log.sd^2) - 1))
}
\description{
Function to calculate the equivalent values for the mean and standard deviation of a normal distribution from the mean and standard deviation of the log-normal distribution and vice versa.  Outputs from this function can be used with the dnorm() or dlnorm functions, and with the equivalent distributions in JAGS.
}
\examples{
lmean <- 2.5
lsd <- 0.2
mean <- normal.params(lmean,lsd)[[1]]
sd <- normal.params(lmean,lsd)[[2]]

curve(dlnorm(x, lmean, lsd), from=0, to=25)
dev.new()
curve(dnorm(x, mean, sd), from=0, to=25)

mean <- 10
sd <- 2
lmean <- lnormal.params(mean,sd)[[1]]
lsd <- lnormal.params(mean,sd)[[2]]

curve(dnorm(x, mean, sd), from=0, to=20)
dev.new()
curve(dlnorm(x, lmean, lsd), from=0, to=20)
}
\keyword{methods}
