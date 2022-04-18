#####################################################################################
## Stochastic predator-prey model with the use of symbolic Jacobian functions in R ##
#####################################################################################

##### Data generation ----
library(Sim.DiffProc)
N=1
T=50
dt=0.125
nT=T/dt

agrowth=1
aconsumed=0.8
bgrowth=0.6
bdeath=0.7

sigma1=0.1
sigma2=0.1

fx <- expression(agrowth*x-aconsumed*x*y,bgrowth*x*y-bdeath*y)
gx <- expression(sigma1,sigma2)
trial <- 1
Seed <- 1127+trial
set.seed(Seed)

repeat{
  r2dtemp<-c(5,1)
  r2d<-r2dtemp
  for (k in seq(dt,T,by=dt)){
    mod2dtemp <- snssde2d(drift=fx,diffusion=gx,M=1,t0=k-dt,x0=as.numeric(r2dtemp),T=k,N=10,type="str",method="smilstein")
    r2dtemp <- rsde2d(mod2dtemp,at=round(k,3))
    r2d<-rbind(r2d,r2dtemp)
  }
  yobs = r2d+rnorm((nT+1)*2,0,0.1)
  
  # Just to make sure that the generated numbers are sensible
if(sum(is.na(r2d))==0){break}else{print("Regenerate")}
}
plot(yobs)

save(Seed,yobs,r2d,file=paste0("SDE_PP_Data",trial,".RData"))


##### Writing the SDE function part of the Stan file with code deriving the analytic Jacobian ----

library(Deriv)

drift <- expression(agrowth*a-aconsumed*a*b,bgrowth*a*b-bdeath*b)

drift_a<-drift[1]
drift_b<-drift[2]

J11 <- Deriv(drift[1],"a")
J12 <- Deriv(drift[1],"b")
J21 <- Deriv(drift[2],"a")
J22 <- Deriv(drift[2],"b")

params_text<-"
  matrix[2,2] sigmaP;
  real agrowth;
  real aconsumed;
  real bgrowth;
  real bdeath;
  real Pvec[3];

  agrowth=params[1];
  aconsumed=params[2];
  bgrowth=params[3];
  bdeath=params[4];

  sigmaP[1,1] = params[5];
  sigmaP[2,1] = params[6];
  sigmaP[1,2] = 0;
  sigmaP[2,2] = params[7];

  a = yp[1];
  b = yp[2];
"

Jacobian_text<-paste0("Jacob[1,1]=",J11,";",
                      "Jacob[1,2]=",J12,";",
                      "Jacob[2,1]=",J21,";",
                      "Jacob[2,2]=",J22,";")

SDEmodel_text_begin<-"
real[] PPSDE(real t,
real[] yp,
real[] params,
real[] x,
int[] x_int) {
  real dyPdt[5];
  real a;
  real b;
  matrix[2,2] Jacob;
  
"


SDEmodel_text_end<-paste0("

  dyPdt[1] = ",drift_a,";
  dyPdt[2] = ",drift_b,";
  dyPdt[3:5] = P_mat_to_array(Jacob*P_array_to_mat(yp[3:5])+P_array_to_mat(yp[3:5])*Jacob'+sigmaP*sigmaP');
  return dyPdt;
}
")

PPSDE_start<-readChar("PPSDE_start.txt",file.info("PPSDE_start.txt")$size)
PPSDE_end<-readChar("PPSDE_end.txt",file.info("PPSDE_end.txt")$size)

SDEmodel_text<-paste0(PPSDE_start,
                      SDEmodel_text_begin,
                      params_text,Jacobian_text,
                      SDEmodel_text_end,
                      PPSDE_end)

writeLines( SDEmodel_text , con="PP-SDE.stan" )

library(rstan)
rstan_options(auto_write = TRUE)

# Check whether the stan script compiles
# stanDso<-stan_model(file="PP-SDE.stan")

#### Running the model ----

# If running multiple chains on a multi-core computer, uncomment this
# options(mc.cores = parallel::detectCores())

# Organize data to be input into stan
stan_data <- list(y = yobs[-1,],
                  T = nT,
                  t0 = 0,
                  ts = matrix(seq(dt,T,by=dt),ncol=1),
                  Lambda = diag(1,2),pstart=c(0,0,0),ystart=c(5,1))


# I ran two chains separately in parallel using a computational cluster
i=1
stansample<-stan(file='PP-SDE.stan',data=stan_data,iter=1000,
                 chains=1,seed=2022+i)
save(stansample,file=paste0("PPSDE_n1_Trial",trial,"_chain",i,".RData"))
i=2
stansample<-stan(file='PP-SDE.stan',data=stan_data,iter=1000,
                 chains=1,seed=2022+i)
save(stansample,file=paste0("PPSDE_n1_Trial",trial,"_chain",i,".RData"))

# Since two chains are ran separately, the following code combines the two chains
j=1
temp_res_store <-list()
for(chain in 1:2){
  resfile <- paste0("PPSDE_n1_Trial",trial,"_chain",chain,".RData")
  if(file.exists(resfile)){
    load(resfile)
    temp_res_store[[j]]<-stansample
    j=j+1}
}
combinechains<-sflist2stanfit(temp_res_store)


## If running two chains together use the following
combinechains<-stan(file='PP-SDE.stan',data=stan_data,iter=1000,
                 chains=2,cores=2,seed=2022)
save(combinechains,file=paste0("PPSDE_n1_Trial",trial,".RData"))

# Examine the result
get_elapsed_time(combinechains)
round(summary(combinechains)[[1]],3)
shinystan::launch_shinystan(combinechains)

# Compare to true data generating values
true_values <- list(agrowth=1,aconsumed=0.8,bgrowth=0.6,bdeath=0.7,
                    sigmaP=c(0.1,0.1),rho=0,sigmaY=c(0.1,0.1))


