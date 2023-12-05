#' Forecasting for the QFSTS mode
#'
#' @param mqbsts An object of the mqbsts class created by a call to the QFSTS_function function.
#' @param STmodel An object of the STModel class created by a call to the tsc.setting function.
#' @param newdata A vector or matrix containing the predictor variables to use in making the prediction. This is only required if the mbsts model has a regression component.
#' @param tau A vecotr containing the quantile values to use in making the prediction.
#' @param steps An integer value describing the number of time steps ahead to be forecasted. If it is greater than the number of new observations in the newdata, zero values will fill in missing new observations.
#'
#' @return Predicting Distribution and mean
#' @export
#'
#' @examples
QFSTS.forecast <-
  function(mqbsts,STmodel=NULL,newdata=NULL,tau=NULL,steps=1){
    ##extract values from mqbsts object
    if(!is.null(STmodel)){
      States<-mqbsts$States[dim(mqbsts$States)[1],,]
      st.sig2<-mqbsts$st.sig2
    }
    if(!is.null(newdata)){
      B.hat<-mqbsts$B.hat
    }
    ob.sig2<-round(mqbsts$ob.sig2, digits = 4)
    #tilde.phi.star<-mqbsts$tilde.phi.star
    W<-mqbsts$W


    if(!is.null(STmodel)){
      ##Initialization of time series components
      ls<-array(0,c(steps,dim(ob.sig2)[1],dim(mqbsts$States)[3]))
      st.err<-array(0,c(dim(st.sig2)[1],steps,dim(mqbsts$States)[3]))
      ####Predict time series components
      for(i in 1:dim(mqbsts$States)[3]){
        for(j in 1:steps){
          if(j==1){
            st.err[,j,i]<-t(mvrnorm(n=1,mu=c(rep(0,dim(st.sig2)[1])),
                                    Sigma=diag(st.sig2[,i],dim(st.sig2)[1])))
            newstate<-STmodel$T[,,1]%*%States[,i]+
              STmodel$R[,,1]%*%st.err[,j,i]
          } else{
            st.err[,j,i]<-t(mvrnorm(n=1,mu=c(rep(0,dim(st.sig2)[1])),
                                    Sigma=diag(st.sig2[,i],dim(st.sig2)[1])))
            newstate<-STmodel$T[,,1]%*%newstate+
              STmodel$R[,,1]%*%st.err[,j,i]
          }
          ls[j,,i]<-t(STmodel$Z[,,1]%*%newstate)
        }
      }
    }

    if(!is.null(newdata)){
      ###check if step is greater than number of observations for newdata
      # if(dim(newdata)[1]>=steps){
      #   newdata<-newdata[1:steps,]
      # } else{
        #newdata<-rbind(newdata,matrix(0,steps-dim(newdata)[1],dim(newdata)[2]))
      #}
      ##Initialization of regression components
      reg<-array(0,c(steps,dim(ob.sig2)[1],dim(B.hat)[3]))
      ##Predict regression components
      for (i in 1:dim(B.hat)[3]){
        reg[,,i]<- newdata%*%B.hat[,,i]
      }
    }
    ####Observation errors
    err<-array(0,c(steps,dim(ob.sig2)[1],dim(ob.sig2)[3]))

    mu_err<-array(0,dim(mqbsts$Phi)[2:3])
    Phi<-mqbsts$Phi
    tilde_phi_tau<-as.vector((1-2*tau)/(tau*(1-tau)))

    for(j in 1:steps){
      for(i in 1:dim(ob.sig2)[3]){
        # err[j,,i]=t(mvrnorm(n=1,mu=W[,i]*Phi.matrix%*%tilde_phi[,i],
        #                     Sigma=W[,i]*ob.sig2[,,i]))
        #err[j,,i]=raml(n=1,mu=mu_err,Sigma=ob.sig2[,,i])
        mu_err[,i]<-Phi[,,i]%*%tilde_phi_tau
        err[j,,i]=raml(n=1,mu=mu_err[,i],Sigma=ob.sig2[,,i])
      }
    }
    if(!is.null(STmodel) & !is.null(newdata)){
      pred.distribution<-ls+reg+err
    } else{
      if(!is.null(STmodel)){
        pred.distribution<-ls+err
      } else{
        pred.distribution<-reg+err
      }
    }
    pred.mean<-apply(pred.distribution,c(1,2),mean)

    return(list(pred.distribution=pred.distribution,pred.mean=pred.mean))
  }
