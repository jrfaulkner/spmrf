#' Generate function for initial parameter values
#'
#' Creates a function for passing to stan which generates initial parameter seeds.
#' @param prior A character string specifying which prior to use on order-k differences. Choices are "horseshoe", "laplace", and "normal".
#' @param likelihood A character string specifying the distribution of the observations. Choices are "normal", "poisson", and "binomial".
#' @param order Numeric value specifying order of differencing (1, 2, or 3).
#' @return A function for generating initial parameter seeds for passage to stan. Parameters are appropriate to a model with specified prior, likelihood, and order.
#' @seealso get_model, rstan::stan
#' @export

get_init <- function(prior="horseshoe", likelihood="normal", order=1) {
  # prior: ("horseshoe", "laplace", "normal")
  # likelihood: ("normal", "poisson", "binomial")
  # order: (1, 2, 3)

	if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  if (!(likelihood %in% c("normal", "poisson", "binomial"))) stop("Must specify likelihood of type 'normal', 'poisson' or 'binomial'.")
	if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")

  ## --- ORDER 1 --- ##
	# Horseshoe
	initf_H_1 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-1, 0.25, .75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-1)
		tth[1] <- om1*zth1 + mt1
		for (zz in 1:(dat$J-1)){
			tu[zz] <- gm*tan(ztau[zz]*pi/2)
			tth[zz+1] <- zdel[zz]*tu[zz] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zgam = zgm,
			ztau = ztau,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			tau = tu
			)
	}'

	##  Laplace (Double Exponential)
	initf_L_1 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-1, 0.25, .75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-1)
		tth[1] <- om1*zth1 + mt1
		for (zz in 1:(dat$J-1)){
			tu[zz] <- gm*sqrt(-2*log(1-zta2[zz]))
			tth[zz+1] <- zdel[zz]*tu[zz] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zgam = zgm,
			ztau2 = zta2,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			tau = tu
			)
	}'


	## Normal  (GMRF with constant precision)
	initf_N_1 <-  'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		tth <- numeric(dat$J)
		tth[1] <- om1*zth1 + mt1
		for (zz in 1:(dat$J-1)){
			tth[zz+1] <- gm*zdel[zz] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zgam = zgm,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1
			)
	}'


	## --- ORDER 2 --- ##
	# Horseshoe
	initf_H_2 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zpt2 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-2, 0.25, .75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- gm*sqrt(1/3)*tan(zpt2*pi/2)
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-2)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		for (zz in 1:(dat$J-2)){
			tu[zz] <- gm*tan(ztau[zz]*pi/2)
			tth[zz+2] <- zdel[zz+1]*tu[zz] + 2*tth[zz+1]-tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zptau2 = zpt2,
			zgam = zgm,
			ztau = ztau,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2,
			tau = tu
			)
	}'



	##  Laplace (Double Exponential)
	initf_L_2 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zpt2 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-2, 0.25, 0.75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- gm*sqrt(-(2/3.0)*log(1-zpt2))
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-1)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		for (zz in 1:(dat$J-2)) {
			tu[zz] <- gm*sqrt(-2*log(1-zta2[zz]))
			tth[zz+2] <- zdel[zz+1]*tu[zz] + 2*tth[zz+1]-tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zptau2 = zpt2,
			zgam = zgm,
			ztau2 = zta2,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2,
			tau = tu
			)
	}'

	##  Normal (GMRF)
	initf_N_2 <-  'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- (gm/sqrt(3))
		tth <- numeric(dat$J)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		for (zz in 1:(dat$J-2)){
			tth[zz+2] <- gm*zdel[zz+1] + 2*tth[zz+1] - tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zgam = zgm,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2
			)
	}'


	## --- ORDER 3 --- ##

	# Horseshoe
	initf_H_3 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zpt2 <- runif(1, 0.25, 0.75)
		zpt3 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		ztau <- runif(dat$J-3, 0.25, .75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- gm*sqrt(1/10)*tan(zpt2*pi/2)
		pt3 <- gm*sqrt(3/10)*tan(zpt3*pi/2)
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-3)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		tth[3] <- pt3*zdel[2] + 2*tth[2] - tth[1]
		for (zz in 1:(dat$J-3)){
			tu[zz] <- gm*tan(ztau[zz]*pi/2)
			tth[zz+3] <- zdel[zz+2]*tu[zz] + 3*tth[zz+2] - 3*tth[zz+1] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zptau2 = zpt2,
			zptau3 = zpt3,
			zgam = zgm,
			ztau = ztau,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2,
			ptau3 = pt3,
			tau = tu
			)
	}'


	##  Laplace (Double Exponential)
	initf_L_3 <- 'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zpt2 <- runif(1, 0.25, 0.75)
		zpt3 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		zta2 <- runif(dat$J-3, 0.25, 0.75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- gm*sqrt(-(1/5.0)*log(1-zpt2))
		pt3 <- gm*sqrt(-(3/5.0)*log(1-zpt3))
		tth <- numeric(dat$J)
		tu <- numeric(dat$J-3)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		tth[3] <- pt3*zdel[2] + 2*tth[2] - tth[1]
		for (zz in 1:(dat$J-3)) {
			tu[zz] <- gm*sqrt(-2*log(1-zta2[zz]))
			tth[zz+3] <- zdel[zz+2]*tu[zz] + 3*tth[zz+2] - 3*tth[zz+1] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zptau2 = zpt2,
			zptau3 = zpt3,
			zgam = zgm,
			ztau2 = zta2,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2,
			ptau3 = pt3,
			tau = tu
			)
	}'

	##  Normal (GMRF)
	initf_N_3 <-  'function(chain_id=1, dat=tmp.dat) {
		LGYSET
		MUYSET
		SDYSET
		ZSIGSET
		zdel <- runif(dat$J-1, -1, 1)
		zth1 <- runif(1, -1, 1)
		zm1 <- runif(1, -1, 1)
		zom1 <- runif(1, 0.25, 0.75)
		zgm <- runif(1, 0.25, 0.75)
		SIGSET
		gm <- 0.1*tan(zgm*pi/2)
		mt1 <- sdy*zm1 + muy
		om1 <- tan(zom1*pi/2)
		pt2 <- gm*sqrt(1/10)
		pt3 <- gm*sqrt(3/10)
		tth <- numeric(dat$J)
		tth[1] <- om1*zth1 + mt1
		tth[2] <- pt2*zdel[1] + tth[1]
		tth[3] <- pt3*zdel[2] + 2*tth[2] - tth[1]
		for (zz in 1:(dat$J-3)){
		   tth[zz+3] <- gm*zdel[zz+2] + 3*tth[zz+2] - 3*tth[zz+1] + tth[zz]
		}
		list(ZSIG
			zdelta =  zdel,
			ztheta1 = zth1,
			zmuth1 = zm1,
			zomega1 = zom1,
			zgam = zgm,
			SIG
			gam = gm,
			theta = tth,
			muth1 = mt1,
			omega1 = om1,
			ptau2 = pt2,
			ptau3 = pt3,
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
		tmp.b <- sub("LGYSET", "", x=tmp.a)
		tmp.b <- sub("MUYSET", "muy <- mean(dat$y)", x=tmp.b)
		tmp.b <- sub("SDYSET", "sdy <- sd(dat$y)", x=tmp.b)
		tmp.b <- sub("ZSIGSET", "zsgm <- runif(1, .25, .75)", x=tmp.b)
		tmp.b <- sub("SIGSET", "sgm <- 5*tan(zsgm*pi/2)", x=tmp.b)
		tmp.b <- sub("ZSIG", "zsigma = zsig,", x=tmp.b)
		tmp.b <- sub("SIG", "sigma = sig,", x=tmp.b)
	}

	if (likelihood=="poisson"){
		tmp.b <- sub("LGYSET", "", x=tmp.a)
		tmp.b <- sub("MUYSET", "muy <- log(mean(dat$y))", x=tmp.b)
		tmp.b <- sub("SDYSET", "sdy <- sd(log(dat$y))", x=tmp.b)
		tmp.b <- sub("ZSIGSET", "", x=tmp.b)
		tmp.b <- sub("SIGSET", "", x=tmp.b)
		tmp.b <- sub("ZSIG", "", x=tmp.b)
		tmp.b <- sub("SIG", "", x=tmp.b)
	}

	bintrans <- 'pp <- dat$y/dat$N
	pp[pp==0.0] <- 0.005
	pp[pp==1.0] <- 0.995
	logitp <- qlogis(pp) \n'


	if (likelihood=="binomial"){
		tmp.b <- sub("LGYSET", bintrans, x=tmp.a)
		tmp.b <- sub("MUYSET", "muy <- qlogis(mean(pp))", x=tmp.b)
		tmp.b <- sub("SDYSET", "sdy <- sd(logitp)", x=tmp.b)
		tmp.b <- sub("ZSIGSET", "", x=tmp.b)
		tmp.b <- sub("SIGSET", "", x=tmp.b)
		tmp.b <- sub("ZSIG", "", x=tmp.b)
		tmp.b <- sub("SIG", "", x=tmp.b)
	}


	ptmp <- parse(text=tmp.b)
	outf <- eval(ptmp)
	return(outf)

}  #end function



