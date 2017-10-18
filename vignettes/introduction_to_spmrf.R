## ----eval=FALSE----------------------------------------------------------
#  library(rstan)

## ---- eval=FALSE---------------------------------------------------------
#  intstall_github("jrfaulkner/spmrf", build_vignettes=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(spmrf)

## ----eval=FALSE,  fig.height=4, fig.width=5, fig.align='center'----------
#  # Set up function
#  nx <- 100
#  xv <- (1:nx)/nx
#  mpc <- rep(25, nx)
#  mpc[xv >= 0.2 & xv < 0.4] <- 10
#  mpc[xv >= 0.4 & xv < 0.6] <- 35
#  mpc[xv >= 0.6] <- 15
#  
#  # Generate data
#  set.seed(3)
#  pc.sd <- 4.5
#  pc.norm <- rnorm(n=nx, mean=mpc, sd=pc.sd)
#  
#  # Create data list for passage to spmrf
#  pcdat.norm <- list(J = nx, y = pc.norm)
#  
#  # Plot function and data
#  plot(1:nx, pc.norm, xlab="t", ylab="y", main="Piecewise Constant Function")
#  lines(1:nx, mpc, lwd=2, col="red")
#  legend('topright', legend="Truth", lty=1, lwd=2, col="red", bty="n")
#  

## ----eval=FALSE----------------------------------------------------------
#  # Parameters to keep
#  # for model with normal prior
#  pars.N <- c("theta", "gam", "sigma")
#  # for model with Laplace prior
#  pars.L <- c("theta", "tau", "gam", "sigma")
#  # for model with horseshoe prior
#  pars.H <- c("theta", "tau", "gam", "sigma")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Compile model objects
#  # -- GMRF
#  ipcfit.N <- spmrf(prior="normal", likelihood="normal", order=1, zeta=0.01,
#  							data=pcdat.norm, chains=0)
#  # -- Laplace
#  ipcfit.L <- spmrf(prior="laplace", likelihood="normal",  order=1, zeta=0.01,
#  							data=pcdat.norm, chains=0)
#  # -- Horseshoe
#  ipcfit.H <- spmrf( prior="horseshoe", likelihood="normal", order=1, zeta=0.01,
#  							data=pcdat.norm, chains=0)
#  		

## ----eval=FALSE----------------------------------------------------------
#  # MCMC run settings
#  nchain <- 4
#  ntotsamp <- 2000
#  nthin <- 5
#  nburn <- 1000
#  niter <- (ntotsamp/nchain)*nthin + nburn
#  
#  # Model with normal prior
#  pcfit.N <- spmrf(prior="normal", likelihood="normal", order=1, fit=ipcfit.N, data=pcdat.norm,
#  							par=pars.N,	chains=nchain, warmup=nburn, thin=nthin, iter=niter,
#  							control=list(adapt_delta=0.96, max_treedepth=12))
#  # Model with Laplace prior
#  pcfit.L <- spmrf(prior="laplace", likelihood="normal",  order=1,  fit=ipcfit.L, data=pcdat.norm,
#  							par=pars.L,	chains=nchain, warmup=nburn, thin=nthin, iter=niter,
#  							control=list(adapt_delta=0.96, max_treedepth=12))
#  # Model with horseshoe prior
#  pcfit.H <- spmrf( prior="horseshoe", likelihood="normal", order=1,  fit=ipcfit.H, data=pcdat.norm,
#  							par=pars.H, chains=nchain, warmup=nburn, thin=nthin, iter=niter,
#  							control=list(adapt_delta=0.995, max_treedepth=12))
#  

## ----eval=FALSE----------------------------------------------------------
#  # extract posterior samples
#  pcout.N <- as.array(pcfit.N)
#  pcout.L <- as.array(pcfit.L)
#  pcout.H <- as.array(pcfit.H)
#  
#  # extract posterior median and 95% BCIs for theta
#  theta.N <- extract_theta(pcfit.N, obstype="normal", alpha=0.05)
#  theta.L <- extract_theta(pcfit.L, obstype="normal", alpha=0.05)
#  theta.H <- extract_theta(pcfit.H, obstype="normal", alpha=0.05)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # example trace plots for horseshoe model
#  plot_trace(pcout.H, "theta[20]", pscale="original", stack=TRUE, colset="black")
#  plot_trace(pcout.H, "tau[20]", pscale="log", stack=TRUE, colset="black")
#  plot_trace(pcout.H, "gam", pscale="log", stack=TRUE, colset="black")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Create plots of posterior trend for all 3 models
#  yrng <- c(0,45)
#  png(filename='pc_plots.png', width=1500, height=420, res=200)
#    par(mfrow=c(1,3), mar=c(2,1.5,1.5,1), oma=c(2,2,0,0))
#    plot_trend(theta=theta.N, obstype="normal", obsvar=pc.norm, xvar=1:nx, pt.cex=0.9, main="Normal",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mpc, lwd=2, lty=2, col="red")
#    plot_trend(theta=theta.L, obstype="normal", obsvar=pc.norm, xvar=1:nx, pt.cex=0.9, main="Laplace",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mpc, lwd=2, lty=2, col="red")
#    plot_trend(theta=theta.H, obstype="normal", obsvar=pc.norm, xvar=1:nx, pt.cex=0.9, main="Horseshoe",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mpc, lwd=2, lty=2, col="red")
#    legend(x="topright", legend=c("Truth", "Median", "95% BCI"), col=c("red", "blue","lightblue"),
#  			lty=c(2,1,1), lwd=c(2,2,3), bty="n", cex=0.8)
#    mtext(side=1, outer=T, line=1, text="t", font=2, cex=0.8)
#    mtext(side=2, outer=T, line=1, text="y", font=2, cex=0.8)
#  dev.off()
#  

## ----eval=FALSE----------------------------------------------------------
#  # Set up generating function
#  nx <- 100
#  xv <- (1:nx)/nx
#  ygfun <- function(x){
#     sin(x) + 2*exp(-30*x^2)
#  }
#  gseq <- seq(-2,2,length=101)
#  mvs1 <- 20 + 10*ygfun(gseq)
#  mvs <- mvs1[-1]
#  
#  # Generate data
#  set.seed(3)
#  vs.sd <- 4.5
#  vs.norm <- rnorm(n=nx, mean=mvs, sd=vs.sd)
#  vsdat.norm <- list(J = nx, y = vs.norm)
#  
#  # Plot function and data
#  plot(1:nx, vs.norm, ylim=c(0,45), xlab="t", ylab="y", col="gray50", main="Varying Smooth Function")
#  lines(1:nx, mvs, lty=2, lwd=2, col="red")
#  legend('topright', legend="Truth", lty=2, lwd=2, col="red", bty="n")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # Compile and run spmrf models for varying smooth (vs) funtion
#  # -- GMRF
#  vsfit.N <- spmrf(prior="normal", likelihood="normal", order=2, data=vsdat.norm,
#  								par=pars.N,	chains=nchain, warmup=nburn, thin=nthin, iter=niter,
#  							  control=list(adapt_delta=0.96, max_treedepth=12))
#  # -- Laplace
#  vsfit.L <- spmrf(prior="laplace", likelihood="normal",  order=2, data=vsdat.norm,
#  								par=pars.L,	chains=nchain, warmup=nburn, thin=nthin, iter=niter, ,
#  							  control=list(adapt_delta=0.96, max_treedepth=12))
#  # -- Horseshoe
#  vsfit.H <- spmrf(prior="horseshoe", likelihood="normal", order=2, data=vsdat.norm,
#  								par=pars.H,	chains=nchain, warmup=nburn, thin=nthin, iter=niter, ,
#  							  control=list(adapt_delta=0.995, max_treedepth=12))
#  
#  # extract posteriors
#  vsout.N <- as.array(vsfit.N)
#  vsout.L <- as.array(vsfit.L)
#  vsout.H <- as.array(vsfit.H)
#  
#  # create posterior summaries for theta					
#  thvs.N <- extract_theta(vsfit.N, obstype="normal", alpha=0.05)
#  thvs.L <- extract_theta(vsfit.L, obstype="normal", alpha=0.05)
#  thvs.H <- extract_theta(vsfit.H, obstype="normal", alpha=0.05)
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  # print parameter summaries
#  print(vsout.N, pars=pars.N)
#  print(vsout.L, pars=pars.L)
#  print(vsout.H, pars=pars.H)
#  
#  # example trace plots for horseshoe model
#  plot_trace(vsout.H, "theta[20]", pscale="original", stack=TRUE, colset="color")
#  plot_trace(vsout.H, "tau[20]", pscale="log", stack=FALSE, colset="black")
#  plot_trace(vsout.H, "gam", pscale="log", stack=FALSE, colset="gray")
#  plot_trace(vsout.H, "sigma", pscale="log", stack=TRUE, colset="black")
#  

## ----eval=FALSE----------------------------------------------------------
#  
#  yrng <- c(0,45)
#  png(filename='vs_plots.png', width=1500, height=420, res=200)
#    par(mfrow=c(1,3), mar=c(2,1.5,1.5,1), oma=c(2,2,0,0))
#    plot_trend(theta=thvs.N, obstype="normal", obsvar=vs.norm, xvar=1:nx, pt.cex=0.9, main="Normal",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mvs, lwd=2, lty=2, col="red")
#    legend(x="topleft", legend=c("Truth", "Median", "95% BCI"), col=c("red", "blue","lightblue"),
#  			lty=c(2,1,1), lwd=c(2,2,3), bty="n", cex=0.8)
#    plot_trend(theta=thvs.L, obstype="normal", obsvar=vs.norm, xvar=1:nx, pt.cex=0.9, main="Laplace",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mvs, lwd=2, lty=2, col="red")
#    plot_trend(theta=thvs.H, obstype="normal", obsvar=vs.norm, xvar=1:nx, pt.cex=0.9, main="Horseshoe",
#  		xlab="", ylab="", ylim=yrng)
#    lines(1:nx, mvs, lwd=2, lty=2, col="red")
#    legend(x="topleft", legend=c("Truth", "Median", "95% BCI"), col=c("red", "blue","lightblue"),
#  			lty=c(2,1,1), lwd=c(2,2,3), bty="n", cex=0.8)
#    mtext(side=1, outer=T, line=1, text="t", font=2, cex=0.8)
#   mtext(side=2, outer=T, line=1, text="y", font=2, cex=0.8)
#  dev.off()
#  
#  

## ----eval=FALSE----------------------------------------------------------
#  # Sampling setting
#  nchain <- 4
#  ntotsamp <- 2000
#  nthin <- 5
#  nburn <- 1000
#  niter <- (ntotsamp/nchain)*nthin + nburn
#  chnlist <- 1:4
#  
#  # Compile fit object
#  ivsfit.H <- spmrf(prior="horseshoe", likelihood="normal", order=2, data=vsdat.norm,
#  						chains=0, par=pars.H)
#  
#  # Run model separately per chain -- one per core
#  tmp.sflist.H <- mclapply(1:nchain, mc.cores = nchain, function(xx) spmrf(prior="horseshoe",
#  							likelihood="normal", order=2, fit=ivsfit.H, data=vsdat.norm, par=pars.H,
#  							chains=1, warmup=nburn,	thin=nthin, iter=niter, chain_id=chnlist[xx], refresh=-1,
#  							control=list(adapt_delta=0.99, max_treedepth=12)) )
#  
#  # Convert list to stanfit object
#  vsfit.H <- sflist2stanfit(tmp.sflist.H)
#  
#  # Extract posteriors into array
#  vsout.H <- as.array(vsfit.H)
#  
#  

