
rm(list=ls())

library(rstan)
rstan_options(auto_write = TRUE)

load("StanDataFerrerH_reduced20_sync.RData")
stan_data <- c(stan_data,rel_tol=1e-3,abs_tol=1e-3,max_step=1e3)

# chain = 2
set.seed(2020*chain)
system.time(odesolver<-stan(file='FerrerHelm_H_diffT2.stan',data=stan_data,iter=1000,
                init=list(list(mu_a1=rnorm(1,0.08*20,0.5),mu_a2=rnorm(1,0.04*20,0.5),
                               mu_b1=rnorm(1,0.03*20,0.5),mu_b2=rnorm(1,0.08*20,0.5),
                               sigma_a1=rnorm(1,0.18*20,0.5),sigma_a2=rnorm(1,0.16*20,0.5),
                               sigma_b1=rnorm(1,0.18*20,0.5),sigma_b2=rnorm(1,0.12*20,0.5))),
                chains=1,
                # cores=2,
                seed=20200322+chain,
                verbose=TRUE,
                sample_file=paste0("Samples_Chain",chain,".csv")))

print(odesolver)

save(odesolver,file=paste0("Res_Ferrer_H_reduced20_sync_chain",chain,".RData"))

q("no")



stanlist<-list()
j=1
for(i in 1:10){
  resfile<-paste0("Res_Ferrer_H_reduced_sync_chain",i,".RData")
  if(file.exists(resfile)){
    load(resfile)
    stanlist[[j]]<-odesolver
    j=j+1
  }
}


combinedchains<-sflist2stanfit(stanlist[1:4])

round(summary(combinedchains)[[1]][,1],3)
write.table(as.data.frame(round(summary(combinedchains)[[1]],3)),
            file="Summary_FerrerDataH15_Reduced.csv",
            sep = ",")


