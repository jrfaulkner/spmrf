#' Fit Bayesian nonparametric adaptive smoothing model
#'
#' Fit Bayesian nonparametric adaptive temporal smoothing models using Hamiltonian Monte Carlo with the \code{stan} engine from the \pkg{rstan} package.
#' @param prior A character string specifying which prior to use on order-\emph{k} differences. Choices are "horseshoe", "laplace", and "normal".
#' @param likelihood A character string specifying the probability distribution of the observation variable. Current choices are "normal", "poisson", and "binomial".
#' @param order Numeric value specifying order of differencing (1, 2, or 3).
#' @param zeta The hyperparameter for the global smoothing parameter gamma.  This is the scale parameter of a half-Cauchy distribution.  Smaller values will result in more smoothing depending on the prior specification and the strength of the data. Values must be > 0.
#' @param fit An instance of S4 class \code{stanfit} derived from a previous fit; defaults to NA. If \code{fit} is not NA, the compiled model associated with the fitted result is re-used; thus the time that would otherwise be spent recompiling the C++ code for the model can be saved.
#' @param data A named list providing the data for the model. See details below.
#' @param pars A vector of character strings specifying parameters of interest; defaults to NA indicating all parameters in the model. If \code{include} = TRUE, only samples for parameters given in \code{pars} are stored in the fitted results. Conversely, if \code{include} = FALSE, samples for all parameters except those given in \code{pars} are stored in the fitted results.
#' @param chains A positive integer specifying number of chains; defaults to 4.
#' @param iter A positive integer specifying how many iterations for each chain (including warmup). The default is 2000.
#' @param warmup A positive integer specifying number of warmup (aka burnin) iterations. This also specifies the number of iterations used for stepsize adaptation, so warmup samples should not be used for inference. The number of warmup should not be larger than \code{iter} and the default is \code{iter}/2.
#' @param thin A positive integer specifying the period for saving sample; defaults to 1.
#' @param ... Additional arguments passed to \code{rstan::stan}.  Do not include \code{file}, \code{model_code}, or \code{init} in this list.
#'
#' @details This function first internally creates \code{stan} model code and a function to generate initial parameter seeds using the information specified in the \code{prior}, \code{likelihood}, \code{order}, and \code{zeta} arguments.  It then passes these to the \code{stan} function.  The \code{bnps} function will take additional arguments passed to \code{stan} \emph{except for} the arguments \code{file}, \code{model_code}, and \code{init}.  These arguments are not accepted because they may conflict with model forms expected by \code{bnps}.  See the \code{stan} function documentation for more details.
#'
#' \code{stan} does all of the work of fitting a Stan model and returning the results as an instance of \code{stanfit}. First, it translates the Stan model to C++ code. Second, the C++ code is compiled into a binary shared object, which is loaded into the current \code{R} session (an object of S4 class \code{stanmodel} is created). Finally, samples are drawn and wrapped in an object of S4 class \code{stanfit}, which provides functions such as \code{print}, \code{summary}, and \code{plot} to inspect and retrieve the results of the fitted model.
#'
#' \code{bnps} can also be used to sample again from a fitted model under different settings (e.g., different \code{iter}) by providing argument \code{fit}. In this case, the compiled C++ code for the model is reused.
#' @return An object of class \code{stanfit}.  See \code{stanfit} and \code{stan} for more details.
#' @references Faulkner, J. R., and V. N. Minin. 2015. Bayesian trend filtering: adaptive temporal smoothing with shrinkage priors. \emph{arXiv preprint}.
#' @seealso \code{\link[rstan]{stan}}, \code{\link[rstan]{stanfit}}, \code{\link{get_model}}, and \code{\link{get_init}}
#' @export


bnps <- function(prior="horseshoe", likelihood="normal", order=1, zeta=0.01, fit=NA, data, pars=NA, chains=4, iter=2000, warmup=floor(iter/2), thin=1, ...)  {

		# check for rstan
		 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
		if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  	if (!(likelihood %in% c("normal", "poisson", "binomial"))) stop("Must specify likelihood of type 'normal', 'poisson' or 'binomial'.")
		if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")
		if (zeta <= 0) stop("zeta must be > 0.")
	  if (is.null(data) | class(data)!="list") stop("Must specify input data set as a list.")

		stlist <- list(...)
		tmp.dat <<- data
	  mcode <- get_model(prior=prior, likelihood=likelihood, order=order, zeta=zeta)
	  finits <- get_init(prior=prior, likelihood=likelihood, order=order)

	  mname <- paste(prior, likelihood, order, sep="_")

	  if ("init" %in% names(stlist) ) {
	  	if (!is.na(stlist$init)){
	  		stop("Do not specify an 'init' file - one will be automatically generated.")
	  	}
	  }

	  if ("file" %in% names(stlist) ) {
	  	if (!is.na(stlist$file)){
	  		stop("External file specified for model code. Use bnps arguments to create model or use previously compiled bnps fit object.")
	  	}
	  }
	  if ("model_code" %in% names(stlist)  ) {
	  	if (!is.na(stlist$file) & stlist$model_code!=""){
	  		stop("Argument 'model_code' specified. Use bnps arguments to create model or use previously compiled bnps fit object.")
	  	}
	  }

	  if (class(fit)=="stanfit") {
	  		#warning("Warning: Make sure 'prior','likelihood', and 'order' match those of 'fit' object!")
	  	  sfit <- stan( model_name=mname, fit=fit, data=data, pars=pars, chains=chains, iter=iter, warmup=warmup, thin=thin, init=finits, ...)
	  }

	  if (class(fit)!="stanfit"){
	  	 if (class(fit)!="logical") stop("Must specify fit as 'stanfit' object or as NA.")
	     if (class(fit)=="logical") {
	     	   if (!is.na(fit)) stop("Must specify fit as 'stanfit' object or as NA.")
	     		 if (is.na(fit)) {
	     		 	    sfit <- stan(model_name=mname, model_code = mcode, fit=fit, data=data, chains=chains, iter=iter, warmup=warmup, thin=thin, init=finits, ...)
	     		 }
	     	}
	  }

	  rm(tmp.dat, envir=.GlobalEnv)

	  return(sfit)
}


