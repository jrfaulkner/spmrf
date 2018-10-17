#' Generate function for initial parameter values
#'
#' Creates a function which generates initial parameter seeds for passing to \code{stan} internally via \code{spmrf}.
#' @param prior A character string specifying which prior to use on order-\emph{k} differences. Choices are "horseshoe", "laplace", and "normal".  Note that "laplace" priors are currently not available for coalescent likelihoods.
#' @param likelihood A character string specifying the distribution of the observations. Choices are "normal", "poisson", "binomial", or "coalescent".
#' @param order Numeric value specifying order of differencing (1, 2, or 3). Note that order 3 is currently not available for coalescent likelihoods.
#' @return An object of class \code{function} for generating initial parameter seeds for passage to \code{stan} via \code{spmrf}. Parameters are appropriate to a model with specified prior, likelihood, and order.
#' @seealso \code{\link{spmrf}}, \code{\link[rstan]{stan}}, \code{\link{get_model}}
#' @export

get_init <- function(prior="horseshoe", likelihood="normal", order=1) {
  # prior: ("horseshoe", "laplace", "normal")
  # likelihood: ("normal", "poisson", "binomial", "coalescent")
  # order: (1, 2, 3)

	if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  if (!(likelihood %in% c("normal", "poisson", "binomial", "coalescent"))) stop("Must specify likelihood of type 'normal', 'poisson', 'binomial', or 'coalescent'.")
	if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")
  if (likelihood=="coalescent" & prior=="laplace") stop("Laplace priors are currently not supported with coalescent likelihoods")
  if (likelihood=="coalescent" & order==3) stop("Order 3 models are currently not supported with coalescent likelihoods")
  
  ## --- ORDER 1 --- ##
	# Horseshoe
	initf_H_1 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-1, 0.25, .75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zgam = zgm,
			ztau = ztau
			)
	}'

	##  Laplace (Double Exponential)
	initf_L_1 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-1, 0.25, .75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zgam = zgm,
			ztau2 = zta2
			)
	}'


	## Normal  (GMRF with constant precision)
	initf_N_1 <-  'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zgm <- runif(1, 0.25, 0.75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zgam = zgm
			)
	}'


	## --- ORDER 2 --- ##
	# Horseshoe
	initf_H_2 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zpt2 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-2, 0.25, .75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zptau2 = zpt2,
			zgam = zgm,
			ztau = ztau
			)
	}'



	##  Laplace (Double Exponential)
	initf_L_2 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zpt2 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-2, 0.25, 0.75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zptau2 = zpt2,
			zgam = zgm,
			ztau2 = zta2
			)
	}'

	##  Normal (GMRF)
	initf_N_2 <-  'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zgm <- runif(1, 0.25, 0.75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zgam = zgm
			)
	}'


	## --- ORDER 3 --- ##

	# Horseshoe
	initf_H_3 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zpt2 <- runif(1, 0.25, 0.75)
		zpt3 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-3, 0.25, .75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zptau2 = zpt2,
			zptau3 = zpt3,
			zgam = zgm,
			ztau = ztau
			)
	}'


	##  Laplace (Double Exponential)
	initf_L_3 <- 'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zpt2 <- runif(1, 0.25, 0.75)
		zpt3 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-3, 0.25, 0.75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zptau2 = zpt2,
			zptau3 = zpt3,
			zgam = zgm,
			ztau2 = zta2
			)
	}'

	##  Normal (GMRF)
	initf_N_3 <-  'function(chain_id=1, dat=tmp.dat) {
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zgm <- runif(1, 0.25, 0.75)
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zgam = zgm
			)
	}'

	## find model template
	options(scipen=25)
	ms <- data.frame(matrix(NA, 9, 3))
	ms[,1] <- c("initf_H_1", "initf_L_1", "initf_N_1",
						"initf_H_2", "initf_L_2", "initf_N_2",
						"initf_H_3", "initf_L_3", "initf_N_3" )
	ms[,2] <- rep(1:3,rep(3,3))
	ms[,3] <- rep(c("horseshoe", "laplace", "normal"), 3)
	names(ms) <- c("modname", "order", "prior")
	mL <- list(initf_H_1, initf_L_1, initf_N_1,
						initf_H_2, initf_L_2, initf_N_2,
						initf_H_3, initf_L_3, initf_N_3 )
	nm <- nrow(ms)
	mm <- (1:nm)[ms$prior==prior & ms$order==order]
	tmp.a <- mL[[mm]]


	# replace likelihood-related statements
	if (likelihood=="normal"){
		tmp.b <- sub("ZSIGSET", "zsig <- runif(1, .25, .75)", x=tmp.a)
		tmp.b <- sub("ZSIG", "zsigma = zsig,", x=tmp.b)
	}

	if (likelihood=="poisson"){
		tmp.b <- sub("ZSIGSET", "", x=tmp.a)
		tmp.b <- sub("ZSIG", "", x=tmp.b)
	}


	if (likelihood=="binomial"){
		tmp.b <- sub("ZSIGSET", "", x=tmp.a)
		tmp.b <- sub("ZSIG", "", x=tmp.b)
	}


	if (likelihood=="coalescent" & prior=="normal" & order==1) {
    tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
     zdel <- rnorm(dat$J-1, 0, 2)
	   th1 <- rnorm(1, dat$log_mu, sd=1.5)
	   zgm <- runif(1, 0.25, 0.75)
	   list(zdelta =  zdel,
	   theta1 = th1,
	   zgam = zgm
	   )
     }'
	}

	if (likelihood=="coalescent" & prior=="normal" & order==2) {
	  tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
        zdel <- rnorm(dat$J-1, 0, 4)
	  th1 <- rnorm(1, dat$log_mu, sd=.5)
	  zpt2 <- runif(1, 0.25, 0.75)
	  zgm <- runif(1, 0.4, 0.6)
	  
	  list(zdelta =  zdel,
	  theta1 = th1, 
	  zptau2 = zpt2,
	  zgam = zgm
	  )
	}' 
	}
	  
	if (likelihood=="coalescent" & prior=="horseshoe" & order==1) {
	 tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
    zdel <- rnorm(dat$J-1, 0, 2)
	 th1 <- rnorm(1, dat$log_mu, sd=1.5)
	 zgm <- runif(1, 0.25, 0.75)
	 ztau <- runif(dat$J-1, 0.25, .75)
	 
	 list(zdelta =  zdel,
	 theta1 = th1,
	 zgam = zgm,
	 ztau = ztau
	 )
	}'
	}

	if (likelihood=="coalescent" & prior=="horseshoe" & order==2) {
	  tmp.b <- 'function(chain_id=1, dat=tmp.dat) {
        zdel <- rnorm(dat$J-1, 0, 4)
	  th1 <- rnorm(1, dat$log_mu, sd=.5)
	  zpt2 <- runif(1, 0.25, 0.75)
	  zgm <- runif(1, 0.4, 0.6)
	  ztau <- runif(dat$J-2, 0.4, .6)
	  
	  list(zdelta =  zdel,
	  theta1 = th1, 
	  zptau2 = zpt2,
	  zgam = zgm,
	  ztau = ztau
	  )
	}'

	}	  
		
	ptmp <- parse(text=tmp.b)
	outf <- eval(ptmp)
	return(outf)

}  #end function



