#' Generate model code for passage to stan
#'
#' Generates a text string of code describing the model of interest for passage to the \code{stan} function in the \pkg{rstan} package.  This function is called internally in the \code{spmrf} function.
#' @param prior A character string specifying which prior to use on order-k differences. Choices are "horseshoe", "laplace", and "normal". Note that "laplace" priors are currently not available for coalescent likelihoods.
#' @param likelihood A character string specifying the probability distribution of the observation variable. Current choices are "normal", "poisson", "binomial", and "coalescent".
#' @param order Numeric value specifying order of differencing (1, 2, or 3). Note that order 3 is currently not available for coalescent likelihoods.
#' @param zeta The hyperparameter for the global smoothing parameter gamma.  This is the scale parameter of a half-Cauchy distribution.  Smaller values will result in more smoothing depending on the prior specification and the strength of the data. Values must be > 0.
#' @return A character string of code describing the model of interest for passage to the function \code{stan} in the \code{rstan} package.
#' @details This function can be used to generate a text string containing code in the \code{stan} model syntax.  This function is called by the \code{spmrf} function internally, so it is not necessary to use \code{get_model} external to \code{spmrf}.
#' @seealso \code{\link{spmrf}},  \code{\link[rstan]{stan}}, \code{\link{get_init}}
#' @export


get_model <- function(prior="horseshoe",  likelihood="normal", order=1,  zeta=0.01){
	# prior: ("horseshoe", "laplace", "normal")
	# likelihood: ("normal", "poisson", "binomial", "coalescent")
	# order: (1, 2, 3)

	if (!(prior %in% c("horseshoe", "laplace", "normal"))) stop("Must specify prior of type 'normal', 'laplace' or 'horseshoe'.")
  if (!(likelihood %in% c("normal", "poisson", "binomial", "coalescent"))) stop("Must specify likelihood of type 'normal', 'poisson', 'binomial', or 'coalescent'.")
	if (!(order %in% c(1,2,3))) stop("Model must be of order 1, 2, or 3.")
  if (likelihood=="coalescent" & prior=="laplace") stop("Laplace priors are currently not supported with coalescent likelihoods")
  if (likelihood=="coalescent" & order==3) stop("Order 3 models are currently not supported with coalescent likelihoods")
	if (zeta <= 0) stop("zeta must be > 0.")

	###  MODEL TEMPLATES  #########
	# These are modified by the getModel() function before passing to stan

	###-----  ORDER 1 -----#####

	###  HORSESHOE PRIOR

	H_1_temp <- '
	  data {
      int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
      vector [N] xvar1;  //locations for observations
      vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
      YSTATE  // response for obs i
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		theta[1] = 5*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
		   tau[j] = gam*tan(ztau[j]*pi()/2);
	 	   theta[j+1] = zdelta[j]*tau[j]*sqrt(duxvar1[j]) + theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		theta[1] = 5*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
	  	   tau[j] = gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+1] = zdelta[j]*tau[j]*sqrt(duxvar1[j]) + theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		theta[1] = 5*sdy*ztheta1 + muy;
		for (j in 1:(J-1)){
	 	   theta[j+1] = gam*zdelta[j]*sqrt(duxvar1[j]) + theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
  	vector <lower=0> [J-2] drat;
	  vector <lower=0> [J-2] sdrat;
		real muy;
		real sdy;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
 	  for (k in 1:(J-2)){
	    drat[k] = duxvar1[k+1]/duxvar1[k]; 
	    sdrat[k] = sqrt(0.5*square(duxvar1[k+1])*(duxvar1[k] + duxvar1[k+1]));
	  }
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*(1/sqrt(3.0))*tan(zptau2*pi()/2);
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
	   	   tau[j] = gam*tan(ztau[j]*pi()/2);
	 	   theta[j+2] = zdelta[j+1]*tau[j]*sdrat[j] + (1+drat[j])*theta[j+1]-drat[j]*theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
		real muy;
		real sdy;
  	vector <lower=0> [J-2] drat;
	  vector <lower=0> [J-2] sdrat;
		//LGYSTATE
		MUYSTATE
		SDYSTATE
 	  for (k in 1:(J-2)){
	    drat[k] = duxvar1[k+1]/duxvar1[k]; 
	    sdrat[k] = sqrt(0.5*square(duxvar1[k+1])*(duxvar1[k] + duxvar1[k+1]));
	  }
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*sqrt(-(2.0/3.0)*log(1-zptau2));
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
		   tau[j] = gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+2] = zdelta[j+1]*tau[j]*sdrat[j] + (1+drat[j])*theta[j+1]-drat[j]*theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
		YSTATE  // response for obs j
		//NSTATE
	  }
	  transformed data {
  	vector <lower=0> [J-2] drat;
	  vector <lower=0> [J-2] sdrat;
		real muy;
		real sdy;
		//LGYSTATE
 	  for (k in 1:(J-2)){
	    drat[k] = duxvar1[k+1]/duxvar1[k]; 
	    sdrat[k] = sqrt(0.5*square(duxvar1[k+1])*(duxvar1[k] + duxvar1[k+1]));
 	  }
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*sqrt(1.0/3.0);
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		for (j in 1:(J-2)){
 	   theta[j+2] = gam*sdrat[j]*zdelta[j+1] + (1+drat[j])*theta[j+1] - drat[j]*theta[j]; 
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*sqrt(1.0/10.0)*tan(zptau2*pi()/2);
		ptau3 = gam*sqrt(3.0/10.0)*tan(zptau3*pi()/2);
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		theta[3] = ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
	   	   tau[j] = gam*tan(ztau[j]*pi()/2);
	 	   theta[j+3] = zdelta[j+2]*tau[j] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*sqrt(-(1.0/5.0)*log(1-zptau2));
		ptau3 = gam*sqrt(-(3.0/5.0)*log(1-zptau3));
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		theta[3] = ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
		   tau[j] = gam*sqrt(-2*log(1-ztau2[j]));
	 	   theta[j+3] = zdelta[j+2]*tau[j] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
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
	    int<lower=0> N; // number of observations
	    int<lower=0> J; // number of grid cells
			vector [N] xvar1;  //locations for observations
	    vector [J-1] duxvar1;  //distances between unique locations
      int<lower=0> xrank1[N]; //rank order of location for each obs
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
		gam = ZETAVAL*tan(zgam*pi()/2);
		ptau2 = gam*sqrt(1.0/10.0);
		ptau3 = gam*sqrt(3.0/10.0);
		theta[1] = 5*sdy*ztheta1 + muy;
		theta[2] = ptau2*zdelta[1] + theta[1];
		theta[3] = ptau3*zdelta[2] + 2*theta[2] - theta[1];
		for (j in 1:(J-3)){
	 	   theta[j+3] = gam*zdelta[j+2] + 3*theta[j+2] - 3*theta[j+1] + theta[j];
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

	
	
	## Coalescent -- Normal Order 1
	N_1_temp_coal <- '
	functions {
	  real coal_loglik_lp(vector ft, vector yy, int nc, vector aik, vector da) {
	    vector [nc] ll;
	    real sll;
	    ll = -1.0*yy .* ft - da .* aik .* exp(-ft); 
	    sll = sum(ll);
	    return sll ;
	  }
	
	}
	
	data {
	  int <lower=1> J; //number of grid points (theta params)
	  int <lower=1> N; //number of grid subsections
	  vector <lower=0>[N] y; //auxillary coal indicator
	  int <lower=1>  gridrep [J]; //number of reps per theta
	  vector <lower=0> [N] Aik; // active lineage combinatoric coefficient
	  vector [N] dalpha; //delta alpha - subgrid widths 
	  real log_mu; //mle for const Ne on log scale
	 }
	
	parameters {
	  vector [J-1] zdelta;
	  real ztheta1;	
	  real <lower=0, upper=1> zgam;
	}
	
	transformed parameters {
	  vector [J] theta;
	  vector [N] ftheta;
	  real <lower=0> gam;
	
	  gam = ZETAVAL*tan(zgam*pi()/2);
	  theta[1] = 10*ztheta1 + log_mu; 
	  for (j in 1:(J-1)){
	    theta[j+1] = gam*zdelta[j] + theta[j]; 
	  }
	  { int cnt;
	    cnt = 0;
	    for (j in 1:J){
	      for (k in 1:gridrep[j]){
	        cnt = cnt + 1;
	        ftheta[cnt] = theta[j];
	      }
	    }
	  }
	}
	
	model {
	  zgam ~ uniform(0, 1);
	  ztheta1 ~ normal(0, 1);
	  zdelta ~ normal(0, 1);
	  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha); 
	} 
	'
	
## Coalescent -- Horseshoe Order 1
H_1_temp_coal <- '
 functions {
    real coal_loglik_lp(vector ft, vector yy, int nc, vector aik, vector da) {
      vector [nc] ll;
      real sll;
      ll = -1.0*yy .* ft - da .* aik .* exp(-ft); 
      sll = sum(ll);
      return sll ;
    }
  }

 data {
  int <lower=1> J; //number of grid points (theta params)
  int <lower=1> N; //number of grid subsections
  vector <lower=0>[N] y; //auxillary coal indicator
  int <lower=1>  gridrep [J]; //number of reps per theta
  vector <lower=0> [N] Aik; // active lineage combinatoric coefficient
  vector [N] dalpha; //delta alpha - subgrid widths 
  real log_mu; //mle for const Ne on log scale
 }

 parameters {
  vector [J-1] zdelta;
  real ztheta1;	
  vector <lower=0, upper=1> [J-1] ztau;
  real <lower=0, upper=1> zgam;
 }

 transformed parameters {
  vector [J] theta;
  vector [N] ftheta;
  real <lower=0> gam;
  vector[J-1] tau;

  gam = ZETAVAL*tan(zgam*pi()/2);
  theta[1] = 10*ztheta1 + log_mu;
  for (j in 1:(J-1)){
    tau[j] = gam*tan(ztau[j]*pi()/2);	
    theta[j+1] = zdelta[j]*tau[j] + theta[j]; 
  }
  { int cnt;
    cnt = 0;
    for (j in 1:J){
      for (k in 1:gridrep[j]){
        cnt = cnt + 1;
        ftheta[cnt] = theta[j];
      }
    }
  }
}

 model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1); 
  ztheta1 ~ normal(log_mu, 10);
  zdelta ~ normal(0, 1);
  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
} 
'

## Coalescent -- Normal Order 2
N_2_temp_coal <- ' 
  functions {
  	real coal_loglik_lp(vector ft, vector yy, int nc, vector aik, vector da) {
    vector [nc] ll;
    real sll;
    ll = -1.0*yy .* ft - da .* aik .* exp(-ft); 
    sll = sum(ll);
    return sll ;
  }
 }

 data {
  int <lower=1> J; //number of grid points (theta params)
  int <lower=1> N; //number of grid subsections
  vector <lower=0>[N] y; //auxillary coal indicator
  int <lower=1>  gridrep [J]; //number of reps per theta
  vector <lower=0> [N] Aik; // active lineage combinatoric coefficient
  vector [N] dalpha; //delta alpha - subgrid widths 
  real log_mu; //mle for const Ne on log scale
 }

 parameters {
  vector [J-1] zdelta;
  real ztheta1;	
  real <lower=0, upper=1> zgam;
 }

 transformed parameters {
  vector[J] theta;
  vector[N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;

  gam = ZETAVAL*tan(zgam*pi()/2);
	theta[1] = 10*ztheta1 + log_mu; 
	theta[2] = sqrt(0.5)*gam*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
    theta[j+2] = gam*zdelta[j+1] + 2*theta[j+1]-theta[j];
  }

  { int cnt;
    cnt = 0;
    for (j in 1:J){
      for (k in 1:gridrep[j]){
        cnt = cnt + 1;
        ftheta[cnt] = theta[j];
      }
    }
  }
 }

 model {
  zgam ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);
  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
 }  
'

## Coalescent -- Normal Order 2
H_2_temp_coal <- ' 
 functions {
	real coal_loglik_lp(vector ft, vector yy, int nc, vector aik, vector da) {
    vector [nc] ll;
    real sll;
    ll = -1.0*yy .* ft - da .* aik .* exp(-ft); 
    sll = sum(ll);
    return sll ;
  }
 }

 data {
  int <lower=1> J; //number of grid points (theta params)
  int <lower=1> N; //number of grid subsections
  vector <lower=0>[N] y; //auxillary coal indicator
  int <lower=1>  gridrep [J]; //number of reps per theta
  vector <lower=0> [N] Aik; // active lineage combinatoric coefficient
  vector [N] dalpha; //delta alpha - subgrid widths 
  real log_mu; //mle for const Ne on log scale
 }

 parameters {
  vector [J-1] zdelta;
  real ztheta1;	
  vector <lower=0, upper=1>[J-2]  ztau;
  real <lower=0, upper=1> zgam;
  real <lower=0, upper=1> zptau2;
 }

 transformed parameters {
  vector[J] theta;
  vector[N] ftheta;
  real <lower=0> gam;
  real <lower=0> ptau2;
  vector[J-2] tau;

  gam = ZETAVAL*tan(zgam*pi()/2);
  ptau2 = (gam/sqrt(2.0))*tan(zptau2*pi()/2);
  theta[1] = 10*ztheta1 + log_mu ;
  theta[2] = ptau2*zdelta[1] + theta[1];
  for (j in 1:(J-2)){
    tau[j] = gam*tan(ztau[j]*pi()/2);
    theta[j+2] = zdelta[j+1]*tau[j] + 2*theta[j+1]-theta[j];
  }

  { int cnt;
    cnt = 0;
    for (j in 1:J){
      for (k in 1:gridrep[j]){
        cnt = cnt + 1;
        ftheta[cnt] = theta[j];
      }
    }
  }
 }

 model {
  zgam ~ uniform(0, 1);
  ztau ~ uniform(0, 1); 
  zptau2 ~ uniform(0, 1);
  ztheta1 ~ normal(0, 1);
  zdelta ~ normal(0, 1);
  target += coal_loglik_lp(ftheta, y, N, Aik, dalpha);
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
		tmp.b <- sub("YSTATE", "real y[N];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", "", x=tmp.b)
		tmp.b <- sub("//NSTATE", "", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy = mean(y);", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy = sd(y);", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "real <lower=0, upper=1> zsigma;", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "real <lower=0> sigma;", x=tmp.b)
		tmp.b <- sub("//SIGSET", "sigma = 5.0*tan(zsigma*pi()/2);", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "zsigma ~ uniform(0,1); ", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "for (i in 1:N){\n y[i] ~ normal(theta[xrank1[i]], sigma); \n}\n", x=tmp.b)
	}
	if (likelihood=="poisson"){
		tmp.b <- sub("YSTATE", "int <lower=0> y[N];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", "real logy[N];\n real ry[N]; \n for (j in 1:N) {\n logy[j] = log(y[j]+0.5);\n ry[j] = 1.0*y[j]; \n}\n", x=tmp.b)
		tmp.b <- sub("//NSTATE", "", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy = log(mean(ry));", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy = sd(logy);", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "for (i in 1:N){\n y[i] ~ poisson_log(theta[xrank1[i]]);\n}\n", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGSET", "", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "", x=tmp.b)

	}
	bintrans <- 'real <lower=0,upper=1> pp[N];
	real logp[N];

	for (i in 1:N){
		pp[i] = (y[i]+0.0)/(sizeN[i]+0.0);
		if (pp[i]==0.0)
			pp[i] = 0.005;
		if (pp[i]==1.0)
			pp[i] = 0.995;
		logp[i] = logit(pp[i]);
	}\n'
	if (likelihood=="binomial"){
		tmp.b <- sub("YSTATE", "int <lower=0> y[N];", x=tmp.a)
		tmp.b <- sub("//LGYSTATE", bintrans, x=tmp.b)
		tmp.b <- sub("//NSTATE", "int <lower=1> sizeN[N];", x=tmp.b)
		tmp.b <- sub("MUYSTATE", "muy = logit(mean(pp));", x=tmp.b)
		tmp.b <- sub("SDYSTATE", "sdy = sd(logp);", x=tmp.b)
		tmp.b <- sub("LIKESTATE" , "for (i in 1:N){\n y[i] ~ binomial_logit(sizeN[i], theta[xrank1[i]]);\n}\n", x=tmp.b)
		tmp.b <- sub("//SIGPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGTPARM", "", x=tmp.b)
		tmp.b <- sub("//SIGSET", "", x=tmp.b)
		tmp.b <- sub("//ZSIGSTATE" , "", x=tmp.b)
	}
	
	if (likelihood=="coalescent" & prior=="normal" & order==1) tmp.b <- N_1_temp_coal
	if (likelihood=="coalescent" & prior=="normal" & order==2) tmp.b <- N_2_temp_coal
	if (likelihood=="coalescent" & prior=="horseshoe" & order==1) tmp.b <- H_1_temp_coal
	if (likelihood=="coalescent" & prior=="horseshoe" & order==2) tmp.b <- H_2_temp_coal
	
	## replace gamma
	tmp.c <- sub(pattern="ZETAVAL", replacement=zeta, x=tmp.b)
	return(tmp.c)
}

















