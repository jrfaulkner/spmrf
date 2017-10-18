\name{munich}
\alias{munich}
\docType{data}
\title{Munich Rental Data}
\description{Cost per square meter for rental units in Munich, Germany in 2003.}
\usage{data("munich")}
\format{
  A data frame with 2,035 observations on the following 4 variables.
  \describe{
    \item{\code{rent}}{a numeric vector indicating the montly cost per square meter of rent in Euros}
    \item{\code{fsize}}{a numeric vector indicating the floor size of a rental unit in square meters}
    \item{\code{year}}{an integer vector indicating the year of construction of a rental unit}
    \item{\code{location}}{an integer vector with an identification number for the district within Munich where a rental unit is located}
  }
}
\details{These data were originally presented by Fahrmeir et al. (2004) and were later analyzed by others including Rue and Held (2005).}

\source{The full data set with more predictor variables is directly available from the R-INLA website: \url{https://folk.ntnu.no/hrue/r-inla.org/examples/munich/Munich.txt}  }

\references{Fahrmeir, L., T. Kneib, and S. Lang. 2004. Penalized structured additive regression for space-time data: a Bayesian perspective. \emph{Statistica Sinica} \bold{14}, 731--761.

Rue, H., and L. Held. 2005. \emph{Gaussian Markov random fields: theory and applications}. Chapman and Hall/CRC.
}
\examples{ \dontrun{
data(munich)

## convert to list for passage to spmrf
rent_dat <- list(y=munich$rent, xvar1=munich$fsize)

## scatterplot of rent vs. floor size
plot(munich$fsize, munich$rent)

}}
\keyword{datasets}
