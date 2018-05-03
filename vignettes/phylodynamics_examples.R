## ----eval=FALSE----------------------------------------------------------
#  library(spmrf)
#  library(rstan)

## ----eval=FALSE----------------------------------------------------------
#  
#  # set random number seed for reproducibility
#  set.seed(5)
#  
#  # Generate simulated data
#  nsamp <- 500   #number of samples total
#  nstart <- 50   #number of samples at time zero
#  nbndry <- 101  #number of grid cell boundaries (number of cells plus 1)
#  samp.end <- 8  #last potential sample time
#  samptv <- c(0, sort(runif(nsamp-nstart, 0, samp.end)) ) #vector of sample times
#  nsampv <- c(nstart, rep(1, nsamp-nstart))  #vector of number sampled
#  # simulate coalescent times
#  coaldat <- coaltimeSim(samp_times = samptv, n_sampled = nsampv,
#  			traj = bottleNeck_traj, lower_bound = 0.1, ne.max=1, ne.min=0.1, bstart=6, bend=4  )
#  
#  # Calculate a grid for estimation
#  sgrid <- makeGrid(coal_times = coaldat$coal_times, samp_times = samptv, Ngrid=nbndry)
#  
#  # Make data set for input to spmrf
#  cdat <-  make_coalescent_data(samp_times = samptv, n_sampled = nsampv, coal_times = coaldat$coal_times, grid = sgrid$grid)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Calculate trajectory over the grid
#  truetraj <- bottleNeck_traj(t = sgrid$midpts, lower_bound = 0.1, ne.max=1, ne.min=0.1, bstart=6, bend=4)
#  
#  # Plot
#  plot(sgrid$midpts, truetraj, xlim=range(rev(sgrid$midpts)), ylim=c(0, 1.2), type="l", col="red",
#       xlab="Time before present", ylab="Effective population size")
#  
#  png(file="vignettes/figure/phylo_bottleneck_traj.png", width=1000, height=600)
#  par(font.axis=2, font.lab=2)
#  plot(sgrid$midpts, truetraj, xlim=range(rev(sgrid$midpts)), ylim=c(0, 1.2), type="l", col="red",
#       xlab="Time before present", ylab="Effective population size", lwd=3)
#  dev.off()
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Set hyperparameter for global scale
#  zeta <- set_zeta_phylo(phylo = coaldat, ncell = 100, alpha = 0.05, order = 1)
#  
#  # Parameters to keep
#  pars.G <- c("theta", "gam")
#  pars.H <- c("theta", "tau", "gam")
#  
#  # MCMC settings
#  nchain <- 4    #number of chains
#  ntotsamp <- 2000  #total number of samples to keep across all chains
#  nthin <- 2        #thinning level
#  nburn <- 1000     #warm-up / burn-in iterations per chain
#  niter <- (ntotsamp/nchain)*nthin + nburn  #total iterations to run
#  
#  # Run models
#  fit.G <- spmrf(prior="normal", likelihood="coalescent", order=1, data=cdat, par=pars.G,
#  					chains=nchain, warmup=nburn, thin=nthin, iter=niter, control=list(adapt_delta=0.98, max_treedepth=12), zeta=zeta)
#  fit.H <- spmrf(prior="horseshoe", likelihood="coalescent", order=1, data=cdat, par=pars.H,
#  					chains=nchain, warmup=nburn, thin=nthin, iter=niter, control=list(adapt_delta=0.995, max_treedepth=12), zeta=zeta)
#  
#  # Extract posterior draws
#  pout.G <- as.array(fit.G)
#  pout.H <- as.array(fit.H)
#  
#  # Get posterior summary for theta
#  th.G <- extract_theta(fit.G, obstype="coalescent")
#  th.H <- extract_theta(fit.H, obstype="coalescent")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Print parameter summaries
#  print(fit.G, pars=pars.G)
#  print(fit.H, pars=pars.H)
#  
#  # Some example trace plots for the horseshoe model
#  plot_trace(pout.H, "theta[10]", pscale="original", stack=TRUE, colset="black")
#  plot_trace(pout.H, "tau[10]", pscale="log", stack=TRUE, colset="black")
#  plot_trace(pout.H, "gam", pscale="log", stack=TRUE, colset="black")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  xrng <- rev(range(sgrid$midpts))
#  yrng <- c(0,1.8)
#  
#  png(filename='phylo_bottleneck_posterior_plots.png', width=1500, height=600, res=200)
#    par(mfrow=c(1,2), mar=c(2,1.5,1.5,1), oma=c(2,2,0,0))
#    plot_trend(theta=th.G, obstype="coalescent", xvar=sgrid$midpts, main="GMRF-1",
#  		xlab="", ylab="", xlim=xrng, ylim=yrng)
#    lines(sgrid$midpts, truetraj, lwd=2, lty=2, col="red")
#    plot_trend(theta=th.H, obstype="coalescent", xvar=sgrid$midpts, main="HSMRF-1",
#  		xlab="", ylab="", xlim=xrng, ylim=yrng)
#    lines(sgrid$midpts, truetraj, lwd=2, lty=2, col="red")
#   legend(x="topright", legend=c("Median", "95% BCI", "Truth"), col=c("blue","lightblue", "red"), lwd=3, lty=c(1,1,2), bty="n", cex=0.8)
#   mtext(side=1, outer=T, line=1, text="Time before present", font=2, cex=0.8)
#   mtext(side=2, outer=T, line=1, text="Effective population size", font=2, cex=0.8)
#  dev.off()
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # load the hcv data
#  
#  
#  
#  
#  

