
#' Make data list for coalescent data
#' 
#' Generates a list containing the data elements necessary to fit a SPMRF model to coalescent data, where effective population sizes are estimated on a uniform grid.
#' This function is based on code written by Michael Karcher for the \code{phylodyn} package and expanded for use with \code{spmrf}.
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param coal_times vector of coalescent times
#' @param grid vector of grid boundaries (length is number of grid cells plus 1)
#'   
#' @return A list containing 
#' \itemize{
#' \item{\code{J}}{ number of grid cells}
#' \item{\code{N}}{ number of subgrid cells}
#' \item{\code{y}}{ vector  of length \code{N} of binary indicators for whether a subgrid cell ends in a coalescent time or not}
#' \item{\code{gridrep}}{ vector  of length \code{J} containing number of subintervals in each grid cell}
#' \item{\code{Aik}} { vector of length \code{N} containing the binomial coefficient or 'coalescent factor' for each subgrid cell. This is equal to choose(x.k, 2), where x.k is the number of active lineages in subgrid cell k. }
#' \item{\code{coalind}} { vector of length \code{N} containing indicator for the coalescent event associated with each subgrid cell}
#' \item{\code{ncoal}} { number of coalescent times, which is equal to the number of samples minus 1}
#' \item{\code{ncoalv}} { vector of length \code{J} containing the number of coalescent events in each grid cell }
#' \item{\code{cstart}} { vector of length \code{ncoal} containing the index of the starting subgrid cell in the range of subgrid cells covering the time between coalescent events }
#' \item{\code{cend}} {vector of length \code{ncoal} containing the index of the ending subgrid cell in the range of subgrid cells covering the time between coalescent events   }
#' \item{\code{dalpha}}  { vector of length \code{N} containing the width of each subgrid cell }
#' \item{\code{rep.idx}} { matrix of \code{J} by 2 where the first column contains for the index for the subgrid cell at which each grid cell starts and the second column contains the index for the subgrid cell at which each grid cell ends.  }
#' \item{\code{log_mu}} { scalar value containing the natural log of the MLE for a constant population size estimated from the coalescent data}
#' }
#' 
#' @export


#J=ng, N=ny, y=y, gridrep=gridrep, Aik=C, coalind=coalind,ncoal=ncoal,ncoalv=ncoalv,cstart=cstart,
#cend=cend,gridstart=grstart,gridend=grend, dalpha=D, rep.idx=rep_idx, log_mu=mle

make_coalescent_data <- function(samp_times, n_sampled, coal_times, grid)
{
  ns <- length(samp_times)
  nc <- length(coal_times)
  ng <- length(grid)-1
  
  if (length(samp_times) != length(n_sampled))
    stop("samp_times vector of differing length than n_sampled vector.")
  
  if (length(coal_times) != sum(n_sampled) - 1)
    stop("Incorrect length of coal_times: should be sum(n_sampled) - 1.")
  
  if (max(samp_times, coal_times) > max(grid))
    stop("Grid does not envelop all sampling and/or coalescent times.")
  
  t <- sort(unique(c(samp_times, coal_times, grid)))  #combined times
  alin <- rep(0, length(t))   #number of active lineages
  
  for (i in 1:ns)
    alin[t >= samp_times[i]] <- alin[t >= samp_times[i]] + n_sampled[i]
  
  for (i in 1:nc)
    alin[t >= coal_times[i]] <- alin[t >= coal_times[i]] - 1
  
  #print(l)
  
  if (sum((alin < 1) & (t >= min(samp_times))) > 0)
    stop("Number of active lineages falls below 1 after the first sampling point.")
  
  mask <- alin > 0
  t <- t[mask]  # drop times with zero lineages
  alin <- head(alin[mask], -1) #drops first
  
  gridrep <- rep(0, ng)  #numbers of subintervals within a grid cell
  for (i in 1:ng)
    gridrep[i] <- sum(t > grid[i] & t <= grid[i+1])
  
  C <- 0.5 * alin * (alin-1)  #binomial coefficient for active lineages
  D <- diff(t)  #time step width
  
  y <- rep(0, length(D))  #indicator for coal event within sub interval
  y[t[-1] %in% coal_times] <- 1
  ny <- length(y)
  ncoal <- sum(y)
  
  rep_idx <- cumsum(gridrep)
  rep_idx <- cbind(start=rep_idx-gridrep+1,end=rep_idx)
  
  coalind <- integer(ny)
  cnt <- 1
  for (j in 1:ny){
    coalind[j] <- cnt
    if (y[j]==1) cnt <- cnt+1
  }
  cstart <- integer(ncoal)
  cend <- integer(ncoal)
  for (k in 1:ncoal){
    tmpi <- which(coalind==k)
    cstart[k] <- min(tmpi)
    cend[k] <- max(tmpi)
  }
  mle <- optimize(f=const_coal_loglike, interval=c(log(.00000001), log(1E6)), y=y, 
                  C=C, D=D, maximum=TRUE)$max
  ncoalv <- numeric(ng)
  for (jj in 1:ng){
    ncoalv[jj] <- sum(y[rep_idx[jj,1]:rep_idx[jj,2]])
  }
  
  return(list(J=ng, N=ny, y=y, gridrep=gridrep, Aik=C, coalind=coalind,ncoal=ncoal,ncoalv=ncoalv, cstart=cstart,cend=cend, dalpha=D, rep.idx=rep_idx, log_mu=mle))
}
