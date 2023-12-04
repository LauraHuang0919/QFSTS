library(LaplacesDemon)
library(tibble)
library(GeneralizedHyperbolic)
library(pscl)
library(MASS)
library(KFAS)
library(MCMCpack)
library(Matrix)
library(BBmisc)
library(reshape)
library(ggplot2)
library(backports)

##################################################################################
##################################################################################
##################################################################################
##################################################################################
sim.data.func<-function(n,tau,corr){
  
  set.seed(100)
  
  ###############Parameter Setup###########
  m<-3 #m: dimension of target series
  mu=c(1,1,1) #mu: include trend or not
  rho=c(0.6,0.3,0.1) #rho: decay value of linear trend
  D=c(0.04,0.05,0.02) #D: global trend
  S=c(100,70,40) #S: number of seasons

  ###############Regression component###########    
  beta<-t(matrix(c(2,3,-2.5,4,0,0,-3.5,2.5,-2,-2,-3,-1,0,0,3,0,-1.5,2, -1.6, 0,0, 0,2,4), nrow=3, ncol=8)) #coefficients for predictors
  
  X1<-rnorm(n,5,2)
  X4<-rnorm(n,-2,5)
  X5<-rnorm(n,-5,2)
  X8<-rnorm(n,0,10)
  X2<-rpois(n, 10)
  X6<-rpois(n, 15)
  X7<-rpois(n, 20)
  X3<-rpois(n, 5)
  X<-cbind(X1,X2,X3,X4,X5,X6,X7,X8) #Predictors
  reg<-X%*%beta #regression componenet
  
  
  ###############Trend component###########
  trend=matrix(0,n,m) #Trend component
  delta=matrix(0,n,m) #Slope
  
  for (j in 1:m){
    for(i in 2:n){
      if (mu[j]==T){
        trend[i,j]<-trend[i-1,j]+delta[i-1,j]+rnorm(1,0,1)
        if(rho[j]!=0){
          delta[i,j]<-D[j]+rho[j]*(delta[i-1,j]-D[j])+rnorm(1,0,1)
        }
      }
    }
  }
 
  ##############Seasonal component###########
  sl=matrix(0,n,m) #Seasonal component

  for (j in 1:m){
    for(i in 1:n){
      if (S[j]!=0){
        if(S[j]>i){
          sl[i,j]<-rnorm(n=1,mean=1,sd=0.5)
        } else{
          sl[i,j]<- -sum(sl[(i-S[j]+1):(i-1),j])+rnorm(n=1,mean=1,sd=0.5)
        }
      }
    }
  }

  ###############Error term###########
  Phi<-diag(c(0.7,0.6,0.9))
  tilde_phi_tau<-as.vector((1-2*tau)/(tau*(1-tau))) #before it was 1-2*tau
  mu_err<-Phi%*%tilde_phi_tau
  Psi_tau<-diag(sqrt(2/(tau*(1-tau))))
  matrix.corr<-matrix(c(1,corr,corr,corr,1,corr,corr,corr,1), nrow=3, ncol=3)
  Sig.err<-Phi%*%Psi_tau%*%matrix.corr%*%Psi_tau%*%Phi
  Sig.err<-round(Sig.err,digits = 4)
  err<-raml(n=n,mu=mu_err,Sigma=Sig.err)   
  
  
  ###############Target series###########
  Y=reg+trend+sl+err 
  colnames(Y)<-c("Y1","Y2","Y3")
  
  ###############Generate Dataset###########
  
  output_dataset<-cbind(Y,X,X,X)
  return (output_dataset)
}
