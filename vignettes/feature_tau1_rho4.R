source("R/tsc.setting.QFSTS.R")
source("R/QFSTS.R")
source("R/QFSTS.forecast.R")
source("vignette/3_series_eg1.R")

tau1<-c(0.1,0.1,0.1)
count=100
for (n in count){
  nam <- paste("simdata_3series_tau1_",n, sep = "")
  assign(nam, sim.data.func(n,tau1,corr=0.4))
}

Ytrain<-as.matrix(simdata_3series_tau1_100[,1:3])
Xtrain<-as.matrix(simdata_3series_tau1_100[,4:27])

STmodel<-tsc.setting(Ytrain,mu=c(1,1,1), #mu: include trend or not
                     rho=c(0.6,0.3,0.1),
                     S=c(100,70,40))

ki<- c(8,16,24)
pii<- matrix(rep(0.5,dim(Xtrain)[2]),nrow=dim(Xtrain)[2])

#beta
b<-matrix(0,dim(Xtrain)[2])
kapp<-0.01

#v0 and V0 for obs Sigma
R2<-0.8
v0<-5

#State component Sigma
v<-0.01
ss<-0.01

Phi<-diag(c(0.7,0.6,0.9))

set.seed(1)
QFSTS.model<-QFSTS.func(Ytrain,Xtrain,STmodel,ki,pii,b,kapp,R2,v0,v,ss,tau1,Phi,mc=400,burn=100)



