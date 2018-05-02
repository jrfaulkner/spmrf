


piecewiseConst_traj <- function(tvec, levs, breaks,...){
  ## assumes num breaks = num levs - 1
  ## assumes first lev lies on interval tvec[1]:breaks[1]
  ## and last lev is on interval breaks[nb]:tvec[nt]
  nt <- length(tvec)
  nl <- length(levs)
  outv <- numeric(nt)
  for (i in 1:nl){
    if (i==1) bint <- tvec < breaks[1]
    if (i==nl) bint <- tvec >= breaks[i-1]
    if (i > 1 & i < nl) bint <- tvec >= breaks[i-1] & tvec < breaks[i]
    outv[bint] <- levs[i]
  }
  return(outv)
}

#' Bottleneck trajectory
#' 
#' Piecewise constant trajectory with two levels where the lower level represents a population bottleneck.
#' @export
bottleNeck_traj <- function(t, ne.max=100, ne.min=10, bstart=2, bend=1,...) {
  out <- rep(0,length(t))
  out[t <= bend] <- ne.max
  out[t > bend & t < bstart] <- ne.min
  out[t >= bstart] <- ne.max
  return(out)
}


#' Boom-bust trajectory
#' 
#' A boom-bust scenario where population experiences a sharp increase and sharp decline. This function has been used elsewhere and is also known as the 'Mexican hat' function.
#' @export
boomBust_traj <- function(x, mf, sf,soff, sscale, snht, bloc, bscale, bht, tmx=NULL){
  if (!is.null(tmx)) {
    if (length(x)==1) {
      if(x > tmx) x <- tmx
    }
    if (length(x)>1) x[x>=tmx] <- tmx
  }	 
  yg <- snht*sin((soff-x)/sscale) + bht*exp(-((x-bloc)^2)/bscale)
  out <- mf + sf*yg
  out
}

#' Piecewise-exponential trajectory
#' 
#' @export
piecewiseExp_traj <- function(t, tbreaks, lnvals, reverse=FALSE, tmx=NULL){
  # tbreaks must include first and last times of t
  # lnvals are log Ne values at sbreaks, tmx is max value of t, above which traj is constant
  ntb <- length(tbreaks)
  lnsv <- numeric(length(t))
  lnsv[1] <- lnvals[1]
  if (!is.null(tmx)) {
    if (length(t)==1) {
      if(t > tmx) t <- tmx
    }
    if (length(t)>1) t[t>=tmx] <- tmx
  }	 
  for (j in 1:(ntb-1)) {
    b1 <- (lnvals[j+1]-lnvals[j])/(tbreaks[j+1]-tbreaks[j])
    b0 <- lnvals[j]-b1*tbreaks[j]
    segid <- t > tbreaks[j] & t <= tbreaks[j+1]
    lnsv[segid] <- b0 + b1*t[segid]
  }
  out <- exp(lnsv)
  if (reverse==TRUE) out <- rev(out)
  return(out)
}

