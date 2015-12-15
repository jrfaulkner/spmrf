# bnps

This repository houses the `R` package `bnps`, which is used for fitting Bayesian nonparametric adaptive smoothing models as described in Faulkner and Minin (2015).  The `bnps` package interfaces with Stan, which is a C++ package for performing Bayesian inference using Hamiltonian Monte Carlo (see [http://mc-stan.org/](http://mc-stan.org/)).  The `bnps` package therefore depends on the `rstan` package to fit models.

## Installation
1. Install package dependency `rstan` and install package `devtools` using `install_packages` function.
2. Load `devtools` using `library(devtools)`.
3. Install `bnps` from GitHub using either
  1. `install_github("jrfaulkner/bnps/")` or
  2. `install_github("jrfaulkner/bnps/", build_vignettes=TRUE)` if you want the vignette documentation which provides examples of using `bnps`.  Note that building vignettes will make the load take a little longer.

## Vignettes
The following vignettes provide some examples using the `bnps` package with step-by-step instructions and R code. 

1. **Introduction_to_bnps**
2. **coal_mine_example**
3. **tokyo_rain_example**

## References
Faulkner, J. R., and V. N. Minin. 2015. Bayesian trend filtering: adaptive temporal smoothing with shrinkage priors. arXiv preprint.
