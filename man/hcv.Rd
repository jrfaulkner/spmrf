\name{hcv}
\alias{hcv}
\docType{data}
\title{Hepatitis C Coalescent Data}
\description{Data describing a fixed phylogenetic tree estimated from genetic sequence data from a Hepatitis C outbreak in Egypt.}
\usage{data("hcv")}
\format{
  A data frame with 62 observations on the following 3 variables.
  \describe{
    \item{\code{samp_times}}{a numeric vector with the unique set of sample times (all zero)}
    \item{\code{n_sampled}}{a numeric vector with the number sampled at each sample time}
    \item{\code{coal_times}}{a numeric vector with the estimated coalescent times}
  }
}
\details{
The Hepatitis C virus (HCV) is a blood-borne RNA viral infection that occurs exclusively in humans. HCV is often asymptomatic, but can lead to liver disease and liver failure. The original data are 63 DNA sequences of type 4 with 411 base pairs from the E1 region of the
HCV genome that were collected in 1993 in Egypt (Ray et al., 2000).  These data have been analyzed by many authors, including Pybus et al. (2003) and Faulkner et al. (2018). The data presented here are coalescent times derived from a maximum clade credibility genealogical tree that was estimated with phylogenetic analysis.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Faulkner, J.R., A.F. Magee, B. Shapiro, and V.N. Minin. 2018. Locally-adaptive Bayesian nonparametric inference for phylodynamics. arXiv preprint.

Pybus, O.G., A.J. Drummond, T. Nakano, B.H. Robertson, and A. Rambaut. 2003. The epidemiology and iatrongenic transmission of hepatitis C virus in Egypt: a Bayesian coalescent approach. Molecular Biology and Evolution 20:381-387.

Ray, S.C., R.R. Arthur, A. Carella, J. Bukh, and D.L. Thomas. 2000. Genetic epidemiology of hepatitis C virus throughout Egypt. The Journal of Infectious Diseases 182:698–707.
}
\examples{
data(hcv)
str(hcv)
make_coalescent_data(hcv)
}
\keyword{datasets}
