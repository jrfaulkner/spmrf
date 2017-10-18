## ----eval=FALSE----------------------------------------------------------
#  # load coal data
#  data(coal)
#  
#  # set up data list for spmrf
#  coal_dat <- list(J = nrow(coal), y = coal$events)
#  

## ----eval=FALSE----------------------------------------------------------
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
#  cfit.N <- spmrf(prior="normal", likelihood="poisson", order=1, data=coal_dat, par=pars.N,
#  					chains=nchain, warmup=nburn, thin=nthin, iter=niter)
#  cfit.L <- spmrf(prior="laplace", likelihood="poisson", order=1, data=coal_dat, par=pars.L,
#  					chains=nchain, warmup=nburn, thin=nthin, iter=niter)
#  cfit.H <- spmrf(prior="horseshoe", likelihood="poisson", order=1, data=coal_dat, par=pars.H,
#  					chains=nchain, warmup=nburn, thin=nthin, iter=niter, control=list(adapt_delta=0.995, max_treedepth=12))
#  
#  # Extract posterior draws
#  cout.N <- as.array(cfit.N)
#  cout.L <- as.array(cfit.L)
#  cout.H <- as.array(cfit.H)
#  
#  # Get posterior summary for theta
#  th.N <- extract_theta(cfit.N, obstype="poisson")
#  th.L <- extract_theta(cfit.L, obstype="poisson")
#  th.H <- extract_theta(cfit.H, obstype="poisson")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Print parameter summaries
#  print(cfit.N, pars=pars.N)
#  print(cfit.L, pars=pars.L)
#  print(cfit.H, pars=pars.H)
#  
#  # Some example trace plots for the horseshoe model
#  plot_trace(cout.H, "theta[10]", pscale="original", stack=TRUE, colset="black")
#  plot_trace(cout.H, "tau[10]", pscale="log", stack=TRUE, colset="black")
#  plot_trace(cout.H, "gam", pscale="log", stack=TRUE, colset="black")
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  yrng <- c(0,5.5)
#  png(filename='coal_plots.png', width=1500, height=400, res=200)
#    par(mfrow=c(1,3), mar=c(2,1.5,1.5,1), oma=c(2,2,0,0))
#    plot_trend(theta=th.N, obstype="poisson", obsvar=coal$events, xvar=coal$year, main="Normal",
#  		xlab="", ylab="", ylim=yrng)
#    plot_trend(theta=th.L, obstype="poisson", obsvar=coal$events, xvar=coal$year, main="Laplace",
#  		xlab="", ylab="", ylim=yrng)
#    plot_trend(theta=th.H, obstype="poisson", obsvar=coal$events, xvar=coal$year, main="Horseshoe",
#  		xlab="", ylab="", ylim=yrng)
#   legend(x="topright", legend=c("Median", "95% BCI"), col=c("blue","lightblue"), lwd=3, bty="n", cex=1)
#   mtext(side=1, outer=T, line=1, text="Year", font=2, cex=0.8)
#   mtext(side=2, outer=T, line=1, text="Accidents per year", font=2, cex=0.8)
#  dev.off()
#  
#  

