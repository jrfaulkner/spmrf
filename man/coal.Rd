\name{coal}
\alias{coal}
\docType{data}
\title{Coal Mine Disaster Data}
\description{Number of coal mine disasters per year for the years 1851-1962.}
\usage{data("coal")}
\format{
  A data frame with 112 observations on the following 2 variables.
  \describe{
    \item{\code{year}}{a numeric vector indicating the years of events}
    \item{\code{events}}{a numeric vector indicating the number of events per year}
  }
}
\details{These data were originally presented as the time intervals between successive coal-mining disasters by Maguire et al. (1952). Later corrections were given by Jarrett (1979) and Raftery and Akman (1986).  We use the version presented by Raftery and Akman (1986) and convert the time intervals to number of disasters per year. A disaster event is defined as an accident involving 10 or more deaths. The first disaster was recorded in March of 1851 and the last in March of 1962, with 191 total event times during the period 1 January, 1851 through 31 December, 1962. Visual inspection of the data suggests a decrease in rate of disasters over time, but it is unclear by eye alone whether this change is abrupt or gradual. The decrease in disasters is associated with a few changes in the coal industry at the time. A sharp decline in labor productivity at the end of the 1880`s is thought to have decreased the opportunity for disasters, and the formation of the Miner`s Federation, a labor union, in late 1889 brought added safety and protection to the workers (Raftery and Akman, 1986).
}
\source{The raw data for time between events is provided in Jarrett (1979). Further corrections were commented on by Raftery and Akman (1986).}
\references{Jarrett, R. G. 1979. A note on the intervals between coal-mining disasters. \emph{Biometrika} \bold{66}, 191--193.

Maguire, B. A., Pearson, E. S, and Wynn, A. H. A. 1952. The time intervals between industrial accidents. \emph{Biometrika} \bold{39}, 168--180.

Raftery, A.E and Akman, V.E. 1986. Bayesian analysis of a Poisson process with a change-point. \emph{Biometrika} \bold{73}, 85--89.
}
\examples{ \dontrun{
data(coal)

## convert to list for passage to spmrf
coal_dat <- list(J=nrow(coal), y=coal$events)

## plot the series with lines
plot(coal$year, coal$events, type="o")

}}
\keyword{datasets}
