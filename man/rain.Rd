\name{rain}
\alias{rain}
\docType{data}
\title{Tokyo Rainfall Data}
\description{Binomial data for the number of years with rain exceeding 1 mm in Tokyo, Japan, where data are binned by day of year, and the years are 1951-1989.}
\usage{data("rain")}
\format{
  A data frame with 366 observations on the following 3 variables.
  \describe{
    \item{\code{day}}{a numeric vector} indicating the day of year. Includes 366 days to accomodate leap years.
    \item{\code{n.years}}{a numeric vector} indicating the number of years a particular day occured.  There are 39 years and 10 of those are leap years, so day 366 only occurs 10 times.
    \item{\code{n.rain}}{a numeric vector} indicating the number of years on which at least 1 mm of rain was recorded for a particular day of year.
  }
}
\details{The raw data were binary indicators of when daily rainfall exceeded 1 mm in Tokyo, Japan, over the course of 39 consecutive years (1951-1989). These indicators were combined by day of year across years, resulting in a sample size of m = 39 for each of 365 out of 366 possible days, and a size of m = 10 for the additional day that occurred in each of the 10 leap years.  These are therefore binomial data where the number of years on each day represent the number of trials and the number of years with rain represent the number of successes.
}
\source{The raw precipitation data from which the counts were derived were obtained from NOAA`s National Center for Climate Information (\url{https://www.ncdc.noaa.gov}). A smaller subset of these data (1983-1984) were initially analyzed by Kitagawa (1987).
}
\references{NOAA`s National Center for Climate Information \url{https://www.ncdc.noaa.gov}.

Kitagawa, G. 1987. Non-Gaussian state-space modeling of nonstationary time series. \emph{Journal of the American Statistical Association} \bold{82}, 1032--1041.
}
\examples{\dontrun{
data(rain)

## create data list for passage to bsmfit:
tokyo_dat <- list(J=nrow(rain), N=rain$n.years, y=rain$n.rain)

}}
\keyword{datasets}
