#' Fit Bayesian nonparametric adaptive smoothing model
#'

bnps <- function(initialize=FALSE, likelihood="normal", prior="normal", order=1, zeta=0.01,
									 file, model_name = "anon_model", model_code = "",
  									fit = NA, data = list(), pars = NA, chains = 4,
  									iter = 2000, warmup = floor(iter/2), thin = 1,
  									init = "random", ...)  {


		  mcode <- get_model(prior=prior, likelihood=likelihood, order=order, zeta=zeta)
		  finits <- get_init(prior=prior, likelihood=likelihood, order=order)





}
