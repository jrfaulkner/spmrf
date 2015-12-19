#' Plot traces of posterior draws of a parameter
#'
#' Plots the posterior draws of a specified parameter by MCMC iteration with options for visual differentiation of separate chains.
#' @param postob An array containing posterior draws from a \code{bnps} model fit or \code{stan} model fit.  Model fit objects of type \code{stanfit} can be coerced into arrays using either the \code{extract} or \code{as.array.stanfit} functions from the \pkg{rstan} package.
#' @param vname Name of variable to plot.  Name must correspond to a dimname in the array.
#' @param pscale Plotting scale. Options are 'original', 'log', or 'inv' (inverse).
#' @param stack Indicator to plot chains separately and stack the plots.
#' @param colset Color scheme for plotting traces by chain.  Schemes 'color' and 'gray' will produce separate colors or gray shades for each chain, and scheme 'black' will make all chain traces black.
#' @return Returns a \code{plot} object
#' @seealso \code{\link[rstan]{stanfit}}, \code{\link[rstan]{extract}}, \code{\link[rstan]{as.array.stanfit}}, \code{\link{bnps}}
#' @export


plot_trace <- function(postob, vname, pscale="original", stack=FALSE, colset="color"){

   # Error checks
   if (missing(postob)) stop("Need to specify object containing posterior draws.")
   if (class(postob)!="array") stop("Posterior draws must be stored in an array.  See 'extract' or 'as.array' functions in rstan.")
   ## add check for vname specification - must be a character and be in dimnames of postob
   if (missing(vname)) stop("Need to specify variable name of parameter to plot.")
   if (!(vname %in% dimnames(postob)$par)) stop("Variable name must be a parameter name in postob.")
   if (!(pscale %in% c("original", "log", "inv"))) stop("Check specification of 'pscale'. Must be either 'original', 'log', or 'inv'.")
   if (!(colset %in% c("color", "gray", "black"))) stop("Check specification of 'colset'. Must be either 'color', 'gray', or 'black'.")

    # set up plots
	nc <- dim(postob)[2] #number of chains
	nit <- dim(postob)[1] #number of iterations per chain
	tnit <- nc*nit  #total iterations
	exv <- postob[, , vname]
	vrng <- range(exv)
	chtck <-
	if (colset=="color") vcv <- rainbow(nc)
	if (colset=="black") vcv <- rep("black", nc)
	if (colset=="gray") {
		glev <- seq(0, 1 - 1/nc, by=1/nc)
		vcv <- gray(glev)
	}
  if (stack==FALSE){
  	if (pscale=="original") {
  		plot(1:tnit, 1:tnit, type="n", ylim=vrng,
  			main=paste("trace of",vname, "by chain") , xlab="iteration", ylab=vname)
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, exv[ ,ii], lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  	if (pscale=="log") {
  		plot(1:tnit, 1:tnit, type="n", ylim=log(vrng),
  			main=paste("trace of log",vname, "by chain") , xlab="iteration", ylab=paste("log",vname) )
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, log(exv[ ,ii]), lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  	if (pscale=="inv") {
  		plot(1:tnit, 1:tnit, type="n", ylim=rev(1/vrng),
  			main=paste("trace of inverse",vname, "by chain") , xlab="iteration", ylab=paste("1/",vname,sep=""))
  		for (ii in 1:nc){
  			tx <- (nit*(ii-1)+1):(ii*nit)
  			lines(tx, 1/exv[ ,ii], lwd=2, col=vcv[ii])
  			if (ii < nc) abline(v=ii*nit, lty=3, col="gray30", lwd=2)
  		}
  	}
  } # end stack false
	if (stack==TRUE){
		par(mfrow=c(4,1), oma=c(2, 2, 2, 0), mar=c(3,2,2,1) )
		if (pscale=="original") {
			for (ii in 1:nc){
			  plot(1:nit, exv[ ,ii], ylim=vrng, main="", xlab="", ylab="",
			  	type="l", lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=vname)
			mtext(side=3, outer=T, line=1, text=paste("trace of",vname, "by chain") )

		}
		if (pscale=="log") {
			for (ii in 1:nc){
			  plot(1:nit, log(exv[ ,ii]), ylim=log(vrng), main="", xlab="", ylab="",
			  	type="l",lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=paste("log",vname))
			mtext(side=3, outer=T, line=1, text=paste("trace of log",vname, "by chain") )
		}
		if (pscale=="inv") {
			for (ii in 1:nc){
			  plot(1:nit, 1/exv[ ,ii], ylim=rev(1/vrng), main="", xlab="", ylab="",
			  	type="l", lwd=2, col=vcv[ii])
			}
			mtext(side=1, outer=T, line=1, text="iteration")
			mtext(side=2, outer=T, line=1, text=paste("1/",vname,sep=""))
			mtext(side=3, outer=T, line=1, text=paste("trace of inverse",vname, "by chain")  )
		}
  }  #end stack true
  par(mfrow=c(1,1), oma=c(0, 0, 0, 0), mar=c(5,4,4,2) )

}


