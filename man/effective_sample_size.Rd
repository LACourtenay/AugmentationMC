% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC_Trace.R
\name{effective_sample_size}
\alias{effective_sample_size}
\title{Calculate Effective Sample Size of MCMC Trace}
\usage{
effective_sample_size(traceObject)
}
\arguments{
\item{traceObject}{an MCMC_Trace object}
}
\value{
Values of the Effective Sample Size (ESS) for each dimension of the trace
}
\description{
Calculation of the effective sample size for mean estimation,
adjusted for autocorrelation
}
\details{
This function calculates the effective sample size for mean estimation
from an MCMC trace by estimating the spectral density of the trace at frequency 0.
Said estimation is performed by fitting an autoregressive model to each dimension
of the trace, framed as a time-series. ESS thus quantifies the number of samples
with the same standard error as the original sample, so as to compute the reliability
of the posterior distribution drawn from the chain.
}
\section{Note}{


The smaller the ESS, the poorer the posterior distribution of the MCMC is considered
to be at estimating the original parameters of the distribution.

A definition of a small ESS is hard to define, with not all authors coming to a
consensus. Likewise, ESS calculations for univariate chains tend to present
much higher values than individual ESS values for each dimension of multiple chains.
For a univariate chain, values above 1000 are recommendable,
while for multivariate chains, values above 100 (or preferably 200) are recommendable
per dimension.

For an alternative ESS calculation, we strongly recommend the functions available
in the mcmcse R library
}

\section{Bibliography}{


Vats, D.; Flegal, J.M.; Jones, G.L. (2015) Multivariate output analysis for Markov
Chain Monte Carlo. arXiv: 1512.07713

Elvira, V.; Martino, L.; Robert, C. (2018) Rethinking the Effective Sample Size.
arXiv: 1809.04129v1
}

\examples{

data("single_moon")

trial_1 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 1)
trial_2 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.5)
trial_3 <- MCMC(as.matrix(single_moon[,1]), 10000, step_size = 0.1)


# example of convergence
effective_sample_size(trial_1)

# example of poor convergence
effective_sample_size(trial_2)

# example of no convergence
effective_sample_size(trial_3)

#

}
\seealso{
\code{\link{MCMC}}, \code{\link{MCMC_Trace}}, \code{\link{burn_in}},
\code{\link{sample_from_trace}},
\code{\link{TOST}}, \code{\link{autocorrelation}},
\code{\link{trace_plot}}
}
\author{
Lloyd A. Courtenay
}
