## ----eval=FALSE----------------------------------------------------------
#  # First make sure spmrf is installed and loaded
#  library(spmrf)
#  
#  # number of grid locations
#  nloc <- 50
#  
#  # number of observations per grid location
#  nobs_per_loc <- sample(1:5, size=nloc, replace=T)
#  
#  # location ID variable
#  locID <- rep(1:nloc, nobs_per_loc)
#  
#  # value of underlying trend function for each obs
#  mu.vec <- 10 + 20*sin(2*pi*locID/nloc )
#  
#  # generate normal obs
#  y.vec <- rnorm(n=length(mu.vec), mean=mu.vec, sd=2)
#  

## ----eval=FALSE----------------------------------------------------------
#  # Set up the data list
#  multi_dat <- list(y = y.vec, xvar1 = locID)
#  

## ----eval=FALSE----------------------------------------------------------
#  # Load munich data
#  data(munich)
#  
#  # Set up data list for spmrf
#  # for using floor size as the covariate
#  mun_dat <- list(y = munich$rent, xvar1 = munich$fsize)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Calculate number of grid locations for covariate
#  nfs <- length(unique(munich$fsize))  # 134 floor space values
#  
#  # Calculate zeta values associated with the covariate
#  zeta_fs <- set_zeta(yvec = munich$rent, linkfun = "identity", ncell = nfs, alpha=0.05, order=2)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  library(rstan)
#  
#  # Parameters to keep
#  pars.N <- c("theta", "gam")
#  pars.L <- c("theta", "tau", "gam")
#  pars.H <- c("theta", "tau", "gam")
#  
#  # MCMC settings
#  nchain <- 4
#  ntotsamp <- 2500
#  nthin <- 5
#  nburn <- 1500
#  niter <- (ntotsamp/nchain)*nthin + nburn
#  
#  # Run models
#  mfit.N <- spmrf(prior="normal", likelihood="normal", order=2, data=mun_dat, par=pars.N,
#  					zeta = zeta_fs, chains=nchain, warmup=nburn, thin=nthin, iter=niter)
#  mfit.L <- spmrf(prior="laplace", likelihood="normal", order=2, data=mun_dat, par=pars.L,
#  					zeta = zeta_fs, chains=nchain, warmup=nburn, thin=nthin, iter=niter)
#  mfit.H <- spmrf(prior="horseshoe", likelihood="normal", order=2, data=mun_dat, par=pars.H,zeta = zeta_fs, chains=nchain, warmup=nburn, thin=nthin, iter=niter, control=list(adapt_delta=0.995, max_treedepth=12))
#  
#  # Extract posterior draws
#  mout.N <- as.array(mfit.N)
#  mout.L <- as.array(mfit.L)
#  mout.H <- as.array(mfit.H)
#  
#  # Get posterior summary for theta
#  th.N <- extract_theta(mfit.N)
#  th.L <- extract_theta(mfit.L)
#  th.H <- extract_theta(mfit.H)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Print parameter summaries
#  print(mfit.N, pars=pars.N)
#  print(mfit.L, pars=pars.L)
#  print(mfit.H, pars=pars.H)
#  
#  # Some example trace plots for the horseshoe model
#  plot_trace(mout.H, "theta[21]", pscale="original", stack=TRUE, colset="black")
#  plot_trace(mout.H, "tau[20]", pscale="log", stack=TRUE, colset="black")
#  plot_trace(mout.H, "gam", pscale="log", stack=TRUE, colset="black")
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # make variable with unique floor sizes sorted
#  ufsize <- sort(unique(munich$fsize))
#  
#  # set range of y-axis
#  yrng <- c(3, 15)
#  
#  png(filename='rent_plots.png', width=1500, height=500, res=200)
#    par(mfrow=c(1,3), mar=c(2,1.5,1.5,1), oma=c(2,2,0,0))
#    plot_trend(theta=th.N, obstype="normal", uxvar=ufsize, main="Normal",	xlab="", ylab="",                     ylim=yrng, trend.lwd=2)
#    points(ufsize, rep(3, length(ufsize)), pch="|", col="gray50")
#    plot_trend(theta=th.L, obstype="normal", uxvar=ufsize, main="Laplace", xlab="", ylab="",                    ylim=yrng, trend.lwd=2)
#    points(ufsize, rep(3, length(ufsize)), pch="|", col="gray50")
#    plot_trend(theta=th.H, obstype="normal", uxvar=ufsize, main="Horseshoe",	xlab="",                          ylab="", ylim=yrng, trend.lwd=2)
#    points(ufsize, rep(3, length(ufsize)), pch="|", col="gray50")
#   legend(x="topright", legend=c("Posterior Median", "95% BCI"), col=c("blue","lightblue"),          lwd=3, bty="n", cex=1)
#   mtext(side=1, outer=T, line=1, text="Floor size", font=2, cex=0.8)
#   mtext(side=2, outer=T, line=1, text="Rent", font=2, cex=0.8)
#  dev.off()
#  
#  

