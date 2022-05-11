
library(rstan)
rstan_options(auto_write = TRUE)

N=30
T=30
dt=0.1
nT=T/dt
trial=1
load(paste0("SimData",N,"T",T,"_",trial,".RData"))

stanDso<-stan_model(file="SampleStanScript.stan")

stan_data <- list(y = r2darray[-1,,],
                  N=N,
                  Tmax=nT,
                  T = rep(nT,N),
                  t0 = 0,
                  ts = array(rep(seq(dt,T,by=dt),N),dim=c(nT,N,1)),
                  Lambda = diag(1,2),
                  pstart=c(0,0,0),ystart=r2darray[1,,],
                  mu1=rep(0,N),mu2=rep(0,N),
                  rel_tol=1e-6,abs_tol=1e-6,max_step=1e6)


stansample<-sampling(stanDso,data=stan_data,iter=1000,warmup=500,
                     chains=2,cores=2,seed=2020+trial)

print(stansample)

shinystan::launch_shinystan(stansample)

