#' Generate model code for passage to stan
#'
#' Generates a text string of code describing the model of interest for passage to the \code{stan} function in the \pkg{rstan} package.  This function is called internally in the \code{bnps} function.
#' @param prior A character string specifying which prior to use on order-k differences. Choices are "horseshoe", "laplace", and "normal".
#' @param likelihood A character string specifying the probability distribution of the observation variable. Current choices are "normal", "poisson", and "binomial".
#' @param order Numeric value specifying order of differencing (1, 2, or 3).
#' @param zeta The hyperparameter for the global smoothing parameter gamma.  This is the scale parameter of a half-Cauchy distribution.  Smaller values will result in more smoothing depending on the prior specification and the strength of the data. Values must be > 0.
#' @return A character string of code describing the model of interest for passage to the function \code{stan} in the \code{rstan} package.
#' @details This function can be used to generate a text string containing code in the \code{stan} model syntax.  This function is called by the \code{bnps} function internally, so it is not necessary to use \code{get_model} external to \code{bnps}.
#' @seealso \code{\link{bnps}},  \code{\link[rstan]{stan}}, \code{\link{get_init}}
#' @export


get_model <- function(prior="horseshoe",  likelihood="normal", order=1,  zeta=0.01){
	# prior: ("horseshoe", "laplace", "normal")
	# likelihood: ("normal", "poisson", "binomial")
	# order: (1, 2, 3)

	if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  if (!(likelihood %in% c("normal", "poisson", "binomial"))) stop("Must specify likelihood of type 'normal', 'poisson' or 'binomial'.")
	if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")
	if (zeta <= 0) stop("zeta must be > 0.")

	###  MODEL TEMPLATES  #########
	# These are modified by the getModel() function before passing to stan

	###-----  ORDER 1 -----#####

	###  HORSESHOE PRIOR

	H_1_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau[J-1];
		real <lower=0, upper=1> zgam;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		vector[J-1] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		theta[1] <- 2*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
		   tau[j] <- gam*tan(ztau[j]*pi()/2);
	 	   theta[j+1] <- zdelta[j]*tau[j] + theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		ztau ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'


	###  LAPLACE PRIOR
	L_1_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau2[J-1];
		real <lower=0, upper=1> zgam;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		vector[J-1] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		theta[1] <- 2*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
	  	   tau[j] <- gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+1] <- zdelta[j]*tau[j] + theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		ztau2 ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'



	###   NORMAL PRIOR (GMRF)

	N_1_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> zgam;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		theta[1] <- 2*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
	 	   theta[j+1] <- gam*zdelta[j] + theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'


	###-----  ORDER 2 -----#####

	###  HORSESHOE

	H_2_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau[J-2];
		real <lower=0, upper=1> zgam;
		real <lower=0, upper=1> zptau2;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		vector[J-2] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*(1/sqrt(3.0))*tan(zptau2*pi()/2);
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
	   	   tau[j] <- gam*tan(ztau[j]*pi()/2);
	 	   theta[j+2] <- zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		zptau2 ~ uniform(0, 1);
		ztau ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'


	###  LAPLACE

	L_2_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau2[J-2];
		real <lower=0, upper=1> zgam;
		real <lower=0, upper=1> zptau2;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		vector[J-2] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*sqrt(-(2.0/3.0)*log(1-zptau2));
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
		   tau[j] <- gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+2] <- zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		zptau2 ~ uniform(0, 1);
		ztau2 ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'


	###   NORMAL
	N_2_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> zgam;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*sqrt(1.0/3.0);
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
	 	   theta[j+2] <- gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'



	###-----  ORDER 3 -----#####

	###  HORSESHOE
	H_3_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau[J-3];
		real <lower=0, upper=1> zgam;
		real <lower=0, upper=1> zptau2;
		real <lower=0, upper=1> zptau3;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		real <lower=0> ptau3;
		vector[J-3] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*sqrt(1.0/10.0)*tan(zptau2*pi()/2);
		ptau3 <- gam*sqrt(3.0/10.0)*tan(zptau3*pi()/2);
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		theta[3] <- ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
	   	   tau[j] <- gam*tan(ztau[j]*pi()/2);
	 	   theta[j+3] <- zdelta[j+2]*tau[j] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
		}
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		zptau2 ~ uniform(0, 1);
		zptau3 ~ uniform(0, 1);
		ztau ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'


	###  LAPLACE
	L_3_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> ztau2[J-3];
		real <lower=0, upper=1> zgam;
		real <lower=0, upper=1> zptau2;
		real <lower=0, upper=1> zptau3;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		real <lower=0> ptau3;
		vector[J-3] tau;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*sqrt(-(1.0/5.0)*log(1-zptau2));
		ptau3 <- gam*sqrt(-(3.0/5.0)*log(1-zptau3));
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		theta[3] <- ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
		   tau[j] <- gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+3] <- zdelta[j+2]*tau[j] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
	    }
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		zptau2 ~ uniform(0, 1);
		zptau3 ~ uniform(0, 1);
		ztau2 ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'




	###   NORMAL
	N_3_temp <- '
	  data {
	    int<lower=0> J; // number of observations
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
	  }
	  parameters {
	    real zdelta[J-1];
	    real ztheta1;
		real <lower=0, upper=1> zgam;
		//SIGPARM
	  }
	  transformed parameters{
	    vector[J] theta;
		real <lower=0> gam;
		real <lower=0> ptau2;
		real <lower=0> ptau3;
		//SIGTPARM

		//SIGSET
		gam <- ZETAVAL*tan(zgam*pi()/2);
		ptau2 <- gam*sqrt(1.0/10.0);
		ptau3 <- gam*sqrt(3.0/10.0);
		theta[1] <- 2*sdy*ztheta1 + muy;
		theta[2] <- ptau2*zdelta[1] + theta[1];
		theta[3] <- ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
	 	   theta[j+3] <- gam*zdelta[j+2] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
	    }
	  }
	  model {
		//ZSIGSTATE
		zgam ~ uniform(0, 1);
		ztheta1 ~ normal(0, 1);
		zdelta ~ normal(0, 1);
		LIKESTATE
	  }
	'

	## find model template
	options(scipen=25)
	ms <- data.frame(matrix(NA, 9, 3))
	ms[,1] <- c("H_1_temp", "L_1_temp", "N_1_temp",
						"H_2_temp", "L_2_temp", "N_2_temp",
						"H_3_temp", "L_3_temp", "N_3_temp" )
	ms[,2] <- rep(1:3,rep(3,3))
	ms[,3] <- rep(c("horseshoe", "laplace", "normal"), 3)
	names(ms) <- c("modname", "order", "prior")
	mL <- list(H_1_temp, L_1_temp, N_1_temp,
						H_2_temp, L_2_temp, N_2_temp,
						H_3_temp, L_3_temp, N_3_temp)
	nm <- nrow(ms)
	mm <- (1:nm)[ms$prior==prior & ms$order==order]
	tmp.a <- mL[[mm]]

	# replace likelihood-related statements
	if (likelihood=="normal"){
		tmp.b <- sub("YSTATE", "real y[J];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", "", x=tmp.b)
		tmp.b <- sub("//NSTATE", "", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy <- mean(y);", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy <- sd(y);", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "real <lower=0, upper=1> zsigma;", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "real <lower=0> sigma;", x=tmp.b)
		tmp.b <- sub("//SIGSET", "sigma <- 5.0*tan(zsigma*pi()/2);", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "zsigma ~ uniform(0,1); ", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "y ~ normal(theta, sigma); ", x=tmp.b)
	}
	if (likelihood=="poisson"){
		tmp.b <- sub("YSTATE", "int <lower=0> y[J];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", "real logy[J];\n real ry[J]; \n for (j in 1:J) {\n logy[j] <- log(y[j]+0.5);\n ry[j] <- 1.0*y[j]; \n}\n", x=tmp.b)
		tmp.b <- sub("//NSTATE", "", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy <- log(mean(ry));", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy <- sd(logy);", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "y ~ poisson_log(theta);", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGSET", "", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "", x=tmp.b)

	}
	bintrans <- 'real <lower=0,upper=1> pp[J];
	real logp[J];

	for (j in 1:J){
		pp[j] <- (y[j]+0.0)/(N[j]+0.0);
		if (pp[j]==0.0)
			pp[j] <- 0.005;
		if (pp[j]==1.0)
			pp[j] <- 0.995;
		logp[j] <- logit(pp[j]);
	}'
	if (likelihood=="binomial"){
		tmp.b <- sub("YSTATE", "int <lower=0> y[J];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", bintrans, x=tmp.b)
		tmp.b <- sub("//NSTATE", "int <lower=1> N[J];", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy <- logit(mean(pp));", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy <- sd(logp);", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "y ~ binomial_logit(N, theta);", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGSET", "", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "", x=tmp.b)
	}

	## replace gamma
	tmp.c <- sub(pattern="ZETAVAL", replacement=zeta, x=tmp.b)
	return(tmp.c)
}

















