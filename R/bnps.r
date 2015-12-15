#' Fit Bayesian nonparametric adaptive smoothing model with stan
#'
#' @param prior A character string specifying which prior to use on order-k differences. Choices are "horseshoe", "laplace", and "normal".
#' @param likelihood A character string specifying the probability distribution of the observation variable. Current choices are "normal", "poisson", and "binomial".
#' @param order Numeric value specifying order of differencing (1, 2, or 3).
#' @param zeta The hyperparameter for the global smoothing parameter gamma.  This is the scale parameter of a half-Cauchy distribution.  Smaller values will result in more smoothing depending on the prior specification and the strength of the data. Values must be > 0.
#' @return An object of class 'stanfit'.
#' @seealso get_model, get_init, rstan::stan
#' @export


bnps <- function(prior="horseshoe", likelihood="normal", order=1, zeta=0.01, ...)  {

		# check for rstan
		 if (!requireNamespace("rstan", quietly = TRUE)) {
    			stop("Package 'rstan' needed for this function to work. Please install it.", call. = FALSE)
  	 }
		if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  	if (!(likelihood %in% c("normal", "poisson", "binomial"))) stop("Must specify likelihood of type 'normal', 'poisson' or 'binomial'.")
		if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")
		if (zeta <= 0) stop("zeta must be > 0.")

		stlist <- list(...)
		tmp.dat <<- stlist$data
	  mcode <- get_model(prior=prior, likelihood=likelihood, order=order, zeta=zeta)
	  finits <- get_init(prior=prior, likelihood=likelihood, order=order)

		# remove any model or init spec from ...
		# check if is an initialization

		sfit <- stan(model_code = mcode, init=finits, ...)
	  rm(tmp.dat, envir=.GlobalEnv)

	  return(sfit)
}


