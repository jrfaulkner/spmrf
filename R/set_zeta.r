#' Calculate a value for the hyperparameter on the global smoothing prior
#' 
#' This function can be used to calculate reasonable values for the hyperparameter \code{zeta}, which controls the scale (and median) of the half-Cauchy prior on the global smoothing parameter for the latent field of trend parameters of an \code{spmrf} model. 
#' 
#' @param yvec A vector of observations on the original scale of measurement.
#' @param mvec For binomial reponse variables only. Is a vector of 'trials' associated with each observed number of 'successes' represented in \code{yvec}.
#' @param linkfun The link function associated with the transformation of the expected value of the response (in a generalized linear models sense). Current options are "identity", "log", "logit", and "probit". 
#' @param ncell The number of grid cells. If there is only one observation per grid location (e.g., observation time or covariate value), then this is equal to the total number of observations. Otherwise is equal to the number of unique location values.
#' @param upBound Upper bound on the expected value of the marginal standard deviations of the latent trend (field) parameters. This value is rarely known \emph{a priori} and here is assumed to equal the standard deviation of the observed data unless otherwise specified.
#' @param alpha The probability of exceeding \code{upBound}. 
#' @param order The order of the SPMRF model (1, 2, or 3).
#' 
#' @details   The   Making \code{alpha} smaller will decrease the size of \code{zeta}, which will result in smoother latent trends if the information in the data does not overcome the prior information. 
#' 
#' The methods for calculation of the hyperparameter \code{zeta} are outlined in Faulkner and Minin (2017) and are based on methods introduced by Sorbye and Rue (2014) for setting hyperparameters for the precision of Gaussian Markov random field priors.

#' 
#' @return A numeric scalar value for the hyperparmeter \code{zeta}, where \code{zeta} > 0.
#' @references Faulkner, J. R., and V. N. Minin. 2017. Locally adaptive smoothing with Markov random fields and shrinkage priors. \emph{Bayesian Analysis} advance publication online.
#' 
#' Sorbye, S. and H. Rue.  2014.  Scaling intrinsic Gaussian Markov random field priors in spatial modelling. \emph{Spatial Statistics} 8:39-51. 
#' @seealso \code{\link{spmrf}}
#' @export


set_zeta <- function(yvec, mvec=NULL, linkfun="identity", ncell=NULL, upBound=NULL, alpha=0.05, order=1){
	if (!(linkfun %in% c("identity", "log", "logit")) ) stop("link function must be either identity, log, or logit")
	if (is.null(ncell)) ncell <- length(yvec)
  if (linkfun=="identity") vld <- var(yvec)
	if (linkfun=="log") vld <- var(log(yvec + 0.5)) 
	if (linkfun=="logit") {
		if (is.null(mvec)) mvec <- rep(1, length(yvec))
		pp <- yvec/mvec ; pp[pp==1] <- 0.995; pp[pp==0] <- 0.005
		vld <- var(qlogis(pp))
	}	
	if (order==1)	cmr <- varRef(ncell, 1,  vld, order=1)
	if (order==2)	cmr <- varRef(ncell, 1,  vld, order=2)
	if (order>=3) stop("Orders >= 3 are currently not supported")
	sref <- exp(mean(0.5*log(cmr)))
  if (is.null(upBound)) upBound <- sqrt(vld)
  zz <- upBound/(sref*(tan((pi/2)*(1-alpha))))
  return(zz)
}


