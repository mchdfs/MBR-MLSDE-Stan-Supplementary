N=30
T=30
dt=0.1
nT=T/dt

trial=1

library(Sim.DiffProc)
Seed=2020+trial
set.seed(Seed)

a1=0.3
a2=0.2
b1=0.8
b2=0

mu1=0
mu2=0
  
sigma1=0.2
  
fx<-expression(a1i*(mu1-x)+a2i*(y-x),b1i*(mu2-y)+b2i*(x-y))
gx <- expression(sigma1,sigma1)
  
  
generatedxstarti<-matrix(NA,nrow=N,ncol=3)
generatedparami<-matrix(NA,nrow=N,ncol=5)
  
r2darray<-array(NA,dim=c(nT+1,2,N))
r2dall<-c()
  
  
for (j in 1:N){
  r2dtemp<-c(-2,0)+rnorm(2,0,c(1,0.1))
  
  generatedxstarti[j,]<-c(j,r2dtemp)
  
  a1i <- rnorm(1,a1,0.1)
  a2i <- rnorm(1,a2,0.2)
  b1i <- rnorm(1,b1,0.2)
  b2i <- rnorm(1,b2,0.1)
    
  generatedparami[j,] <- c(j,a1i,a2i,b1i,b2i)
  
  r2d<-r2dtemp
  for (k in seq(dt,T,by=dt)){
    mod2dtemp <- snssde2d(drift=fx,diffusion=gx,M=1,t0=k-dt,x0=as.numeric(r2dtemp),T=k,N=10,type="str",method="smilstein")
    r2dtemp <- rsde2d(mod2dtemp,at=round(k,3))
    r2d<-rbind(r2d,r2dtemp)
  }
  r2darray[,1,j]<-r2d$x
  r2darray[,2,j]<-r2d$y
  r2dall<-rbind(r2dall,cbind(r2d,id=j))
}
  
r2darray[-1,,] = r2darray[-1,,]+rnorm((nT)*N*2,0,0.5)
save(Seed,r2dall,r2darray,generatedparami,generatedxstarti,file=paste0("SimData",N,"T",T,"_",trial,".RData"))
  