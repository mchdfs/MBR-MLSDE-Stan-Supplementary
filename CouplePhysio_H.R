
rm(list=ls())

library(rstan)
rstan_options(auto_write = TRUE)

# Number of chains
nchains <- 2
# Number of iterations of the MCMC sampling (warmup + sampling)
niter <- 1000

# Optional - initial values for the MCMC algorithm
initpara <- list()
j <- 0
repeat{  
  initpara <- append(initpara,
                     list(mu_a1=rnorm(1,0.08*20,0.5),
                          mu_a2=rnorm(1,0.04*20,0.5),
                          mu_b1=rnorm(1,0.03*20,0.5),
                          mu_b2=rnorm(1,0.08*20,0.5),
                          sigma_a1=rnorm(1,0.18*20,0.5),
                          sigma_a2=rnorm(1,0.16*20,0.5),
                          sigma_b1=rnorm(1,0.18*20,0.5),
                          sigma_b2=rnorm(1,0.12*20,0.5)))
  j <- j+1
  if(j == nchains){break}
}


# Load Data
load("StanDataFerrerH_reduced20_sync.RData")

stan_data <- c(stan_data,rel_tol=1e-3,abs_tol=1e-3,max_step=1e3)

# Optional - testing whether the stan file can be compiled successfully
# stanDso<-rstan::stan_model(file="CouplePhysio_H.stan")

combinedchains<-stan(file='CouplePhysio_H.stan',data=stan_data,iter=niter,
                init=list(initpara),
                chains=nchains,
                # cores=nchains, # for parallel computation
                # seed=2020, # for reproducible chains
                verbose=TRUE)
# Quick print of the results
print(combinedchains)

# Sumamry table
round(summary(combinedchains)[[1]][,1],3)

# Interactive shiny app with the MCMC samples
launch_shinystan(combinedchains)

q("no")





