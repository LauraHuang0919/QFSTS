# install.packages("GeneralizedHyperbolic")
# install.packages("pscl")
# install.packages("MASS")
# install.packages("KFAS")
# install.packages("MCMCpack")
# install.packages("Matrix")
# install.packages("BBmisc")
# install.packages("reshape")
library(GeneralizedHyperbolic)
library(pscl)
library(MASS)
library(KFAS)
library(MCMCpack)
library(Matrix)
library(BBmisc)
library(reshape)
library(fBasics)
library(LaplacesDemon)
library(forecast)
# Y=Ytrain
# X.star=Xtrain

#' Main function for the quantile feature selection time series.
#'
#' @param Y  A (n ∗ m)-dimensional matrix containing multiple target series, where n is thenumber of observations and m
#'           is the number of target series.
#' @param X.star A (n ∗ K)-dimensional matrix containing all candidate predictor series for each target series.
#'                K = $\Simga$ k_i is the number of all candidate predictors for all target
#                 series. The first k1 variables are the set of candidate predictors for the first target
#                 series, and the next k2 variables are the set of candidate predictors for the second
#                 target series, etc. Note that, one variable can appear in the X.star several times,
#                 since different target series can contain the same candidate predictors.
#' @param STmodel A state space model of STmodel class returned by tsc.setting.
#' @param ki A vector of integer values denoting the acumulated number of predictors for
#             target series. For example, if there are three target series where the first has 8
#             predictors, the second has 6 predictors, and the third has 10 predictors, then the
#             vector is c(8, 14, 24).
#' @param pii A vector describing the prior inclusion probability of each candidate predictor.
#' @param b NULL or a vector describing the prior means of regression coefficients. The default value is NULL.
#' @param kapp A scalar value describing the number of observations worth of weight on the prior mean vector. The default value is 0.01.
#' @param R2 A numerical value taking value in [0, 1], describing the expected percentage of variation of Y to be explained by the model. The default value is 0.8.
#' @param v0 A numerical value describing the prior degree of freedom of the inverse Wishart distribution for $\sigma_\epsilon$.
#' @param v A numerical value describing the prior degree of freedom of the inverse Wishart distribution for (Σµ, Σδ, Στ , Σω). The default value is 0.01
#' @param ss A numerical value describing the prior scale matrix of the inverse Wishart distribution for (Σµ, Σδ, Στ , Σω). The default value is 0.01.
#' @param tau A vecotr containing the quantile values.
#' @param Phi A m × m-dimensional diagonal matrix of means of error terms.
#' @param mc A positive integer giving the desired number of MCMC draws. The default value is 500
#' @param burn A positive integer giving the number of initial MCMC draws to be discarded. The default value is 50.
#'
#' @return An object of QFSTS class
#' @export
#'
#' @examples
QFSTS.func<-
  function(Y,X.star=NULL,STmodel=NULL,
           ki=NULL,pii=NULL,
           b=NULL,kapp=0.1,
           R2=0.8,v0=NULL,
           v=0.01,ss=0.01,tau=NULL,Phi=NULL,
           mc=500,burn=50){
    n=nrow(Y)       #number of observations
    m=ncol(Y)       #number of response variables
    if(!is.null(X.star)){
      I=c(rep(1,dim(X.star)[2]))    #Initialization of indicator function
    }
    # mc=2
    # burn=1
    #Initialization of Sigma_eplison
    if(is.null(STmodel)){
      Sigma2<-diag(0.1,m)
    } else {
      Sigma2<-matrix(STmodel$H,m,m) #$#$#$#$ this is Sigma_eplison
    }

    ######Initialization of output results
    #####time series components
    if(!is.null(STmodel)){
      States<-array(0,c(dim(Y)[1],dim(STmodel$R)[1],mc-burn))
      st.sig2<-matrix(0,dim(STmodel$R)[2],mc-burn)
    }
    ###regression component
    if(!is.null(X.star)){
      Ind<- matrix(0,dim(X.star)[2],mc-burn)
      beta.hat<- matrix(0,dim(X.star)[2],mc-burn)
      B.hat<-array(0,c(dim(X.star)[2],m,mc-burn))
    }
    ob.sig2<-array(0,c(m,m,mc-burn))  ##Sigma_eplison
    W.record<-array(0,c(1,mc-burn))
    Phi.record<-array(0,c(m,m,mc-burn))

    W<-rexp(1,1)


    tilde_phi<-as.vector((1-2*tau)/(tau*(1-tau)))
    tilde_Phi_tau_matrix<-t(matrix(rep(tilde_phi,n),nrow=m))
    mu_err<-Phi%*%tilde_phi
    mu_err.tilde<-rep(mu_err,each=n)


    for(jj in 1:mc){
      #jj=1
      if(!is.null(STmodel)){  #check existence of time series components
        ### draw latent state####
        LS<- simulateSSM(STmodel, type = "states")    #draw latent states by smoother
        LS<-matrix(LS,ncol=dim(LS)[2])
        State.error<- simulateSSM(STmodel,type = "eta")  #draw state disturbance
        Obs.error<- matrix(simulateSSM(STmodel,type = "epsilon"),ncol=m)  #draw residuals
        State.error<-matrix(State.error,ncol=dim(State.error)[2])

        ##### Draw state component parameters ###equation 30 in mbsts###
        v.post <-  v+n     #posterior degree freedom
        ss.post<-vector()
        State.sigma2<-vector()
        for (i in 1:ncol(State.error)) {
          ss.post[i] <- ss+crossprod(State.error[,i])/W       #posterior sum of square
          State.sigma2[i]<-rigamma(1,alpha=v.post,beta=ss.post[i]) #draw state parameters
        }
      }

      #####Transformation end##############
      if(!is.null(X.star)){ #check existence of regression components
        ####SSVS for drawing gamma, sigma_eplison, beta######
        ##transform design matrix to a larger matrix for computation convenience
        K=ncol(X.star)         #number of predictors
        Xindex<-seq(1,K)  #column index for each predictor in X
        for (i in 1:K){
          X.star[,i]<-X.star[,i]-mean(X.star[,i])  #demean predictors
        }
        Xtemp1<-X.star[,1:ki[1]]
        Xtemp2<-X.star[,(ki[1]+1):ki[2]]
        X<-as.matrix(bdiag(Xtemp1,Xtemp2))
        if(m>2){
          for (j in 3:m){
            Xtemp<-X.star[,(ki[j-1]+1):ki[j]]
            X<-as.matrix(bdiag(X,Xtemp))
          }
        }

        if(is.null(STmodel)){
          Y.star<- Y
        } else {
          Y.star<- Y-LS%*%t(matrix(STmodel$Z,nrow=m))
          #substract time series component from target series
        }
        for (i in 1:m){
          Y.star[,i]<-Y.star[,i]-mean(Y.star[,i])  #demean response variable
        }
        Y.tilde<-matrix(Y.star,ncol = 1)   #transform to vector form



        #transformed system with uncorrelated errors
        U<-chol(Sigma2)    #Cholesky decomposition for observation variance parameter
        U_inv<-solve(U)
        UI<-t(U_inv)%x%diag(n)
        X.hat<-UI%*%X            #transformed X
        Y.hat<-UI%*%Y.tilde      #trasnformed y
        mu_err.hat<-UI%*%mu_err.tilde  #trasnformed mu_err

        ###Prior parameters setting
        V0=(v0-m-1)*(1-R2)*var(Y.star)    #prior scale matrix for sigma_epsilon

        A<-kapp*(t(X)%*%X)/n   #prior information matrix for beta

        ###All posterior parameters
        V<-1/W*t(X.hat)%*%X.hat+A   #posterior co-variance matrix for beta
        N=n+v0       #posterior degree freedom for sigma_epsilon
        gama<-I     #assign value to temporary gamma

        #####draw each gamma ######
        for(j in sample(Xindex))
        {
          zero<-0
          ####To make sure at least one predictor selected
          if(sum(gama[-j])==0){
            zero<-1            #except jth predictor, other predictors not selected
            Index<-sample(Xindex[-j],1)
            gama[Index]<-1     #randomly choose one predictor except jth predictor
          }

          p<-vector()
          for(value in 1:0)
          {
            gama[j]<-value
            p.gamma<- prod(pii^gama)*prod((1-pii)^(1-gama))
            #prior probability for gamma
            b.gamma<-b[which(gama==1),]
            #prior coefficients for selected predictors
            X.gamma<-X.hat[,which(gama==1)] #design matrix with selected predictors
            ####prior parameters
            A.gamma<-A[which(gama==1),which(gama==1)]
            #prior information matrix with selected predictors

            ####posterior parameters
            V.gamma<-V[which(gama==1),which(gama==1)]
            #posterior co-variance matrix for beta with selected predictors
            Z.gamma<-1/W*t(X.gamma)%*%Y.hat+A.gamma%*%b.gamma-t(X.gamma)%*%mu_err.hat

            exp.par<-as.vector(t(b.gamma)%*%A.gamma%*%b.gamma-
                                 t(Z.gamma)%*%solve(V.gamma)%*%Z.gamma)
            ##Set constant value to shrink value for exponent##
            if(value==1){
              temp<-exp.par
            }
            ###posterior pmf for gamma
            p[value+1]<-(det(as.matrix(A.gamma))/det(as.matrix(V.gamma)))^(1/2)*
              p.gamma*exp(-0.5*(exp.par-temp))
          }
          ##debug due to computation round off
          if (p[2]/(p[1]+p[2])>1){
            gama[j]<-rbinom(1, 1, prob=1)
          } else {
            gama[j]<-rbinom(1, 1, prob=p[2]/(p[1]+p[2]))    #update the gamma[j]
          }
          if(zero==1){
            gama[Index]<-0
          }
        }
        I<-gama    #assign updated gamma to indicator function

        ###All predictors not selected (not checked and changed)
        if(sum(I)==0){
          Obs.sigma2<-cov(Obs.error)     #assign value to sigma^2_eplison
          beta<-matrix(0,nrow=K,ncol=1)  #assign zero to all coefficients
          B<-matrix(0,nrow=K,ncol=m)
        }  else{     ###At least one predictor selected
          ######draw sigma and beta############
          b.gamma<-b[which(I==1),]
          ####prior coefficients for selected predictors
          X.gamma<-X.hat[,which(I==1)]   #design matrix with selected predictors
          Xstar.gamma<-X.star[,which(I==1)]   #design matrix with selected predictors
          ####prior parameters
          A.gamma<-A[which(I==1),which(I==1)]
          ####prior information matrix with selected predictors
          ####posterior parameters
          V.gamma<-V[which(I==1),which(I==1)]
          beta.sigma<-solve(V.gamma)
          #posterior co-variance matrix for beta with selected predictors
          beta.tilde<-beta.sigma%*%(1/W*t(X.gamma)%*%Y.hat-t(X.gamma)%*%mu_err.hat+A.gamma%*%b.gamma)

          beta<-matrix(0,nrow=K,ncol=1)     #initialization of beta
          beta.temp<-mvrnorm(n=1,mu=beta.tilde,Sigma=beta.sigma)
          #draw beta from multivariate normal distribution
          beta[which(I==1),]<-beta.temp
          ##beta is generated with new gamma


          ####transfrom beta to matrix form
          betai1<-beta[1:ki[1],]
          betai2<-beta[(ki[1]+1):ki[2],]
          B<-as.matrix(bdiag(betai1,betai2))
          if(m>2){
            for (j in 3:m){
              betai<-beta[(ki[j-1]+1):ki[j],]
              B<-as.matrix(bdiag(B,betai))
            }
          }
          B.gamma<-matrix(B[which(I==1),],ncol=m)  #estimated coefficients with selected predictors
          E.gamma.temp<-Y.star-Xstar.gamma%*%B.gamma
          Phi.inv<-solve(Phi)
          mu_err.star<-matrix(mu_err.tilde,nrow = n, ncol = m, byrow = TRUE)
          E.gamma<-(E.gamma.temp-mu_err.star*W)%*%Phi.inv   #observed error terms
          #E.gamma<-(E.gamma.temp-mu_err.star*W)  #observed error terms
          SS.post<-1/W*t(E.gamma)%*%E.gamma+V0    #posterior sum of squares matrix
          Obs.sigma2.temp<-riwish(N,SS.post)  #draw observation variance parameter
          Obs.sigma2<-Phi%*%Obs.sigma2.temp%*%Phi
          #Obs.sigma2<-riwish(N,SS.post)


          U2<-chol(Obs.sigma2.temp)
          U_inv2<-solve(U2)
          #V.scale<-diag(1,m)
            x.record<-array(0,c(m,m,1000))
            for(i_ind in 2:1000){
              #i_ind=2
              x.record[,,1] = Phi
              currentx = x.record[,,i_ind-1]
              #proposedx=diag(c(0.1,0.2,0.3))
              proposedx = currentx + diag(mvrnorm(1,mu=rep(0,m),Sigma=diag(rep(0.01,m))))
              #proposedx = currentx + rtmvnorm(n = 1, mu=-Phi[i,i], sigma=1, lb = 0, ub=10000)
              proposedx.inv<-solve(proposedx)
              currentx.inv<-solve(currentx)
              #resultA=det(proposedx)^{-n}*exp(-1/2*1/W*tr(t(E.gamma.temp%*%proposedx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)
                                                          #%*%(E.gamma.temp%*%proposedx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W) ))
              #resultB=det(currentx)^{-n}*exp(-1/2*1/W*tr(t(E.gamma.temp%*%currentx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)
                                                          #%*%(E.gamma.temp%*%currentx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W) ))
              resultA=det(currentx)^{n}*exp(-1/2*tr(1/W*t(E.gamma.temp%*%proposedx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)
                                                          %*%(E.gamma.temp%*%proposedx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)))
              resultB=det(proposedx)^{n}*exp(-1/2*tr(1/W*t(E.gamma.temp%*%currentx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)
                                                         %*%(E.gamma.temp%*%currentx.inv%*%U_inv2-tilde_Phi_tau_matrix%*%U_inv2*W)))


              #AAA = (/)*(/resultB)
              #if(AAA==Inf){AAA=1}
              if(resultA>resultB){
                x.record[,,i_ind] = proposedx       # accept move with probabily min(1,A)
              } else {
                x.record[,,i_ind] = currentx        # otherwise "reject" move, and stay where we are
              }
            }
            Phi=x.record[,,1000]



          mu_err<-Phi%*%tilde_phi
          mu_err.tilde<-rep(mu_err,each=n)
          UI<-t(U_inv)%x%diag(n)
          mu_err.hat<-UI%*%mu_err.tilde  #trasnformed mu_err
          W.a<-2+t(mu_err.hat)%*%mu_err.hat #need to use the new generated mu_err.hat
          NewResidue<-Y.hat-X.hat%*%beta
          W.b<-t(NewResidue)%*%NewResidue
          W.p<-1-n/2
          W<-rgig(1, W.b, W.a, W.p, param = c(W.b, W.a, W.p))
          #https://www.rdocumentation.org/packages/GeneralizedHyperbolic/versions/0.8-4/topics/Generalized%20Inverse%20Gaussian
          #W is generated with new Simga_epsilon, mu_err, Y-X_{gamma}beta_{gamma}
        }
      }


      ####assign values to output results during each iteration
      if(jj>burn){
        if(!is.null(STmodel)){
          States[,,(jj-burn)]<-LS
          st.sig2[,(jj-burn)]<-State.sigma2
          if(is.null(X.star)){
            Obs.sigma2<-cov(Obs.error)     #assign value to sigma^2_eplison
            ob.sig2[,,jj-burn]<-Obs.sigma2
          }
        }
        if(!is.null(X.star)){
          Ind[,(jj-burn)]<- matrix(I,nrow=dim(X.star)[2])
          beta.hat[,(jj-burn)]<- matrix(beta,nrow=dim(X.star)[2])
          B.hat[,,(jj-burn)]<-B
          ob.sig2[,,jj-burn]<-Obs.sigma2
        }
        W.record[,(jj-burn)]<-W
        Phi.record[,,(jj-burn)]<-Phi
      }

      #####updating state space model parameters
      if(!is.null(STmodel)){
        if(!is.null(X.star)){
          rowname<-row.names(STmodel$y)
          STmodel$y[]<- ts(Y-X.star%*%B,names = colnames(STmodel$y))
          ####target series
          row.names(STmodel$y)<-rowname
          ####variance-covariance matrix for observation errors
          Sigma2 <- matrix(Obs.sigma2,m,m)
          STmodel$H[]<-array(Sigma2,c(m,m,1))
        } else {
          ####variance-covariance matrix for observation errors
          Sigma2 <- matrix(cov(Obs.error),m,m)
          STmodel$H[]<-array(Sigma2,c(m,m,1))
        }
        ####variance-covariance matrix for time series components
        STmodel$Q[]<- array(diag(State.sigma2,nrow = dim(STmodel$R)[2])
                            ,c(dim(STmodel$R)[2],dim(STmodel$R)[2],1))
      } else{
        Sigma2 <- matrix(Obs.sigma2,m,m)
      }
    }

    #####return output results
    if(!is.null(X.star) & !is.null(STmodel)){
      return (list(Ind=Ind, beta.hat=beta.hat, B.hat=B.hat,
                   ob.sig2=ob.sig2, States=States, st.sig2=st.sig2, W=W.record, Phi=Phi.record))
    } else {
      if(!is.null(STmodel)){
        return (list(ob.sig2=ob.sig2, States=States, st.sig2=st.sig2, W=W.record, Phi=Phi.record))
      } else {
        return (list(Ind=Ind,beta.hat=beta.hat,B.hat=B.hat,ob.sig2=ob.sig2, W=W.record, Phi=Phi.record))
      }
    }
  }
