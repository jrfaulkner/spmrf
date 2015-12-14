#' Extract latent process parameters and associated Bayesian credible intervals
#'
#' Extract theta (latent process) parameters from an object with posterior draws and associated Bayesian credible intervals.
#' @param mfit An object containing posterior draws from a \code{bnps} model fit or \code{rstan::stan} model fit.  The object can be of class 'stanfit', 'array', 'matrix', or 'data.frame'.
#' @param obstype Distribution name of observations.  This controls the back-transformation of the process parameters.  Possible obstypes are 'normal', 'poisson', or 'binomial'.
#' @param alpha Controls level for (1-alpha) Bayesian credible intervals.
#' @return Returns the posterior median and posterior (1-alpha) quantiles of theta.
#' @seealso rstan, rstan::as.array, rstan::as.matrix, rstan::as.data.frame
#' @export

extract_theta <- function(mfit, obstype="normal", alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a bnps or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class stanfit, array, matrix, or data.frame.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit")  tmp.th <- extract(mfit, "theta")[[1]]
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	tmp.th1 <- mfit[ , 1, ]
	if (nca > 1){
	  for (jj in 2:nca){
	    tmp.th1 <- rbind(tmp.th, mfit[ ,jj,])
	  }
	}
	tmp.th <- tmp.th1[ ,grep(x=dimnames(tmp.th1)[[2]], pattern="theta")]
	dimnames(tmp.th) <- dimnames(tmp.th)[[2]]
  }
  if (class(mfit)[1]=="matrix") tmp.th <- mfit[ , grep(x=dimnames(mfit), pattern="theta")]
  if (class(mfit)[1]=="data.frame") tmp.th <- as.matrix(mfit[ , grep(x=names(mfit), pattern="theta")])

  plow <- alpha/2
  phigh <- 1 - alpha/2

  if (obstype=="normal"){
   tmp.md <- apply(tmp.th, 2, median)
   tmp.l <- apply(tmp.th, 2, quantile, probs=plow)
   tmp.u <- apply(tmp.th, 2, quantile, probs=phigh)
  }
  if (obstype=="poisson"){
   tmp.md <- exp(apply(tmp.th, 2, median))
   tmp.l <- exp(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- exp(apply(tmp.th, 2, quantile, probs=phigh))
  }
  if (obstype=="binomial"){
   tmp.md <- plogis(apply(tmp.th, 2, median))
   tmp.l <- plogis(apply(tmp.th, 2, quantile, probs=plow))
   tmp.u <- plogis(apply(tmp.th, 2, quantile, probs=phigh))
  }
   out <- list(postmed=tmp.md, bci.lower=tmp.l, bci.upper=tmp.u)
   out

}

