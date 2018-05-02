

#' Classic Skyline Estimates
#' 
#' Classic Skyline estimates for isochronous or heterochronous data. 
#' 
#' @param phy a list containing a numeric vector of coalescent times (\code{coal_times}), a numeric vector of sampling times (\code{samp_times}, and a integer vector of number of samples taken at each sampling time (\code{n_sampled}).
#'   
#' @return A list containing a vector of indices \code{k}, a vector of midpoints between coalescent times \code{mid_ctime}, and a vector of skyline population estimates \code{theta}.
#' 
#' @export
#' 
skyLine <- function(phy){
  # classic skyline for isochronous or heterochronous data.
  # phy is a list from summarize_phylo() with 
  # coal_times, samp_times, and n_sampled
  # outputs vector of Ne estimates sorted from most recent.
  st <- phy$samp_times 
  ns <- phy$n_sampled 
  ct <- phy$coal_times 
  sm <- data.frame(type='s', time=st, nadd=ns)
  cm <- data.frame(type='c', time=ct, nadd=-1)
  zm <- merge(sm, cm, all=T)
  zm <- zm[order(zm$time), ]
  zm$ncount <- cumsum(zm$nadd)
  zn <- nrow(zm)
  zuid <- numeric(zn)
  zuid[zm$type=='c'] <- 1:nrow(cm)
  for (j in (zn-1):1) {
    if (zuid[j]==0) zuid[j] <- zuid[j+1]
  }
  zm$zuid <- zuid
  wk <- diff(zm$time)
  dm <- data.frame(uid=zm$zuid[-1], ncount=zm$ncount[-zn], wk=round(wk,8) )
  dm$thsub <- dm$ncount*(dm$ncount-1)*dm$wk/2
  adm <- aggregate(dm[,"thsub"], by=list(uid=dm$uid), sum)
  zc1 <- c(0,cm$time[-nrow(cm)])
  zc2 <- cm$time
  ctmid <- zc1 + (zc2-zc1)/2
  
  return(list(k=adm$uid, mid_ctime=ctmid, theta=adm$x))
}


#' Make grid for coalescent data
#' 
#' Makes grid for effective population size estimation from coalescent data.
#' @export
makeGrid <- function(coal_times, samp_times, Ngrid){
  gbds <- range(c(coal_times, samp_times))
  grd <- seq(gbds[1],gbds[2],length=Ngrid)
  intl <- grd[2] - grd[1]
  mids <- grd[-1] - intl/2
  return(list(grid=grd, midpts=mids))	
}


## log likelihood for a constant Ne
const_coal_loglike <- function(theta, y, C, D){
  -theta*sum(y) - exp(-theta)*sum(C*D)
} 


## squared exponential cov fun for ihPois_GP
#' @export
expCovih <- function(ssv, pvec){
  # plist is list of parms {sigmaf2, ll}
  sigmaf2 <- pvec[1]
  ll <- pvec[2]
  nn <- length(ssv)
  cmat <- matrix(0, nn, nn)
  for (i in 1:nn){
    for (j in 1:nn){
      cmat[i,j] <- sigmaf2*exp(-0.5*((i-j)/ll)^2)
    }
  }
  cmat
}

muGPconst <- function(ssv, cmu=0){
  rep(cmu, length(ssv))
}



#' Simulate from inhomogeneous, heterochronous coalescent.
#'
#' For simulating coalescent times using the method of thinning. 
#' This function is based on code written by Michael Karcher for the \code{phylodyn} package.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param lower_bound numeric lower limit on \code{traj} function on its support.
#' @param ... additional arguments to be passed to \code{traj} function.
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' 
#' @examples
#' coaltimeSim(0:2, 3:1, unif_traj, lower_bound=10, level=10)
#' @export
coaltimeSim <- function(samp_times, n_sampled, traj, lower_bound, ...)
{
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower_bound/traj(time, ...))
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}



#' Simulate from inhomogeneous, heterochronous coalescent driven by Gaussian process.
#' 
#' For simulating coalescent times using the method of thinning. 
#' This function is based on code written by Michael Karcher for the \code{phylodyn} package.
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param trajvec vector of population sizes (should be dense) 
#' @param lower_bound numeric lower limit on \code{traj} function on its support.
#' @param tvec grid boundaries for trajvec. Length is length(trajvec)+1
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
coaltimeSimGP <- function(samp_times, n_sampled, trajvec, tvec, lower_bound)
{
  # trajvec has gp trajectory (should be dense)
  # tvec is grid boundaries for trajvec (length(trajvec) + 1 )
  coal_times = NULL
  lineages = NULL
  
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  
  tup <- tvec[-1]
  tlo <- tvec[-length(tvec)]
  tint <- tvec[2]-tvec[1]
  tmids <- tvec[-1]-tint/2
  tadj <- diff(range(tvec))/(10*(length(tvec)-1))
  
  
  while (time <= max(samp_times) || active_lineages > 1)
  {
    if (active_lineages == 1)
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    
    time = time + rexp(1, 0.5*active_lineages*(active_lineages-1)/lower_bound)
    # convert time to index of trajvec
    dtime <- time
    if (sum(time==tlo)>0) dtime <- time + tadj  #adjujst for cell boundaries
    if (time >= max(tvec)) dtime <- max(tvec) - tadj	 # set long times equal to max
    tind <- which(dtime > tlo & dtime < tup)
    
    if (curr < length(samp_times) && time >= samp_times[curr + 1])
    {
      curr = curr + 1
      active_lineages = active_lineages + n_sampled[curr]
      time = samp_times[curr]
    }
    else if (runif(1) <= lower_bound/trajvec[tind] )
    {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
  }
  
  return(list(coal_times = coal_times, lineages = lineages,
              intercoal_times = c(coal_times[1], diff(coal_times)),
              samp_times = samp_times, n_sampled = n_sampled))
}


