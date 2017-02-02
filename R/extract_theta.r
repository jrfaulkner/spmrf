#' Extract latent process parameters and associated Bayesian credible intervals
#'
#' Extract theta (latent process) parameters from an object with posterior draws and associated Bayesian credible intervals.
#' @param mfit An object containing posterior draws from a \code{spmrf} model fit or \code{stan} model fit.  The object can be of class \code{stanfit}, \code{array}, \code{matrix}, or \code{data.frame}.
#' @param obstype Character string with the name of the probability distribution of the observations.  This controls the back-transformation of the process parameters.  Possible values for \code{obstype} are 'normal', 'poisson', or 'binomial'.
#' @param alpha Controls level for 100*(1-\code{alpha})\% Bayesian credible intervals. Values must be 0 < \code{alpha} < 1.
#' @return Returns a list with the posterior median and posterior (1-\code{alpha}) quantiles of the theta parameter vector.
#' @seealso \code{\link[rstan]{stan}}, \code{\link[rstan]{as.array.stanfit}}, \code{\link[rstan]{as.matrix.stanfit}}, \code{\link[rstan]{as.data.frame.stanfit}}, \code{\link{spmrf}}
#' @export

extract_theta <- function(mfit, obstype="normal",  alpha=0.05){

  if (missing(mfit)) stop("Must specify object with posterior draws from a spmrf or stan model fit.")
  if ( !(class(mfit)[1] %in% c("array", "matrix", "data.frame", "stanfit") ) ) stop("Object must be of class 'stanfit', 'array', 'matrix', or 'data.frame'.  Object must be or be generated from a stan model fit object.")
  if ( !(obstype %in% c("normal", "poisson", "binomial") ) ) stop("Argument 'obstype' must be 'normal', 'poisson', or 'binomial'.")
  if (!(0 < alpha & alpha < 1)) stop("Must specify 'alpha' between 0 and 1.")

  if (class(mfit)[1]=="stanfit") {
  	 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
  	tmp.th <- rstan::extract(mfit, "theta")[[1]]
  }
  if (class(mfit)[1]=="array")  {
    nca <- dim(mfit)[2]
	  tmp.th1 <- mfit[ , 1, ]
	  if (nca > 1){
	    for (jj in 2:nca){
	      tmp.th1 <- rbind(tmp.th1, mfit[ ,jj,])
	    }
	  }
	  ath <-  grep(x=dimnames(tmp.th1)[[2]], pattern="theta")
	  zth <-  grep(x=dimnames(tmp.th1)[[2]], pattern="ztheta")
	  thind <- setdiff(ath, zth)
	  tmp.th <- tmp.th1[ , thind]
  }
  if (class(mfit)[1]=="matrix"){
  	ath <-  grep(x=colnames(mfit), pattern="theta")
  	zth <-  grep(x=colnames(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- mfit[ , thind]
  }
  if (class(mfit)[1]=="data.frame") {
  	ath <-  grep(x=names(mfit), pattern="theta")
  	zth <-  grep(x=names(mfit), pattern="ztheta")
  	thind <- setdiff(ath, zth)
  	tmp.th <- as.matrix(mfit[ , thind])
  }
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

