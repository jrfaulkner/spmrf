% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_zeta_phylo.r
\name{set_zeta_phylo}
\alias{set_zeta_phylo}
\title{Set zeta hyperparameter for coalescent data}
\usage{
set_zeta_phylo(phylo, ncell, upBound = NULL, alpha = 0.05, order = 1)
}
\arguments{
\item{phylo}{A list containing a numeric vector of coalescent times (\code{coal_times}), a numeric vector of sampling times (\code{samp_times}, and a integer vector of number of samples taken at each sampling time (\code{n_sampled})}

\item{ncell}{The number of grid cells from the uniformly spaced grid over which effective population size is to be estimated.}

\item{upBound}{Upper bound on the expected value of the marginal standard deviations of the log latent trend (field) parameters. This value is rarely known \emph{a priori} and here is assumed to equal the standard deviation of the observed data (skyline estimates in this case) unless otherwise specified.}

\item{alpha}{The probability of exceeding \code{upBound}.}

\item{order}{The order of the SPMRF model (1 or 2).}
}
\value{
A numeric scalar value for the hyperparmeter \code{zeta}, where \code{zeta} > 0.
}
\description{
Calculate a value for the hyperparameter on the global smoothing prior for coalescent data
}
\details{
This function can be used to calculate reasonable values for the hyperparameter \code{zeta}, which controls the scale (and median) of the half-Cauchy prior on the global smoothing parameter for the latent field of trend parameters (on the log scale) of an \code{spmrf} model applied to coalescent data.

  Making \code{alpha} smaller will decrease the size of \code{zeta}, which will result in smoother latent trends if the information in the data does not overcome the prior information. 

The methods for calculation of the hyperparameter \code{zeta} are outlined in Faulkner and Minin (2018) and are based on methods introduced by Sorbye and Rue (2014) for setting hyperparameters for the precision of Gaussian Markov random field priors.
}
\references{
Faulkner, J. R., and V. N. Minin. 2018. Locally adaptive smoothing with Markov random fields and shrinkage priors. \emph{Bayesian Analysis} 13(1):225-252.

Sorbye, S. and H. Rue.  2014.  Scaling intrinsic Gaussian Markov random field priors in spatial modelling. \emph{Spatial Statistics} 8:39-51.
}
\seealso{
\code{\link{spmrf}}, \code{\link{set_zeta}}
}
