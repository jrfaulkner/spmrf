#' Calculate reference marginal variance of field parameters
#' 
#' This is an internal function used by the \code{set_zeta} function
#' 
#' @keywords internal




varRef <- function(nn, kap, omg2, order=1){
	#nn is number of thetas, kap = 1/gam^2 (precision), omg2=omega^2 (variance of theta1)
	#returns vector of marginal variances
	nnv <- 1:nn
	if (order==1) out <- omg2 + (nnv-1)*(1/kap)
	if (order==2) out <- omg2 + (1/kap)*nnv*(nnv-1)*(2*nnv-1)/6
	if (order >= 3) out <- NULL
	return(out)
}





