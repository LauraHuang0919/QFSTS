# QFSTS: Quantile Feature Selection Time Series Model

Quantile feature selection over correlated multivariate time series data has always been a methodological challenge and is an open problem. In this package, we propose a general Bayesian dimension reduction methodology for feature selection in high-dimensional joint quantile time series analysis, under the name of the quantile feature selection time series (QFSTS) model. This QFSTS model is an extension from Dr. Ning's previous [mbsts](https://cran.r-project.org/web/packages/mbsts/index.html) package with similar functionality but adding more settings on the quantile feature selection. The QFSTS model is a general structural time series model, where each component yields an additive contribution to the time series modeling with direct interpretations. Feature selection is conducted in the quantile regression component, where each time series has its own pool of contemporaneous external predictors allowing nowcasting.

## Installation 

```{r install, tidy='formatR',eval=FALSE, echo=TRUE}
devtools::install_github("LauraHuang0919/QFSTS", build_vignettes = TRUE)
```

## Usage

```{r attach, echo=T, results='hide', message=F, warning=F, tidy='formatR'}
library(QFSTS)
```
### 1. Specification of time series components

Generation of the initial time series components in the mbsts package is through the tsc.setting function. The other input parameters of the tsc.setting function are the following: The trend inclusion parameter mu and the learning rate parameter rho for the trend component; The seasonality parameter S for the seasonal component; The damping factor parameter vrho and the frequency parameter lambda for the cycle component.

```{r conversion, tidy='formatR', tidy.opts=list(width.cutoff = 70),cache=T}
STmodel<-tsc.setting(Ytrain,mu=c(1,1,1), #mu: include trend or not
                     rho=c(0.6,0.3,0.1),
                     S=c(100,70,40))
```
### 2. Model Training

Model training with the QFSTS model is performed through the QFSTS_func function in the QFSTS package. Model training uses the MCMC approach, which is to sample from a probability distribution based on constructing a Markov chain that has the desired distribution as its equilibrium distribution. That is, the states of the chain after discarding some steps as “burn-in” data, are used as samples from the desired distribution. Specifically, the MBSTS model uses the Gibbs sampler for feature selection utilizing the classical “spike and slab” prior setup (George and McCulloch, 1997). Gibbs sampler can be seen as a special case of the Metropolis–Hastings algorithm. The point of the Gibbs sampler is that given a multivariate distribution it is simpler to sample from a conditional distribution than to marginalize by integrating over a joint distribution. To implement Gibbs sampler in model training, all necessary conditional probabilities were derived. In the illustration example, for saving time, we run 10 MCMC iterations while discarding the results of the first 5 iterations. However, logically speaking, to get a good MCMC result, we need to run much more than that. 


```{r conversion, tidy='formatR', tidy.opts=list(width.cutoff = 70),cache=T}
QFSTS.model<-QFSTS.func(Ytrain,Xtrain,STmodel,ki,pii,b,kapp,R2,v0,v,ss,tau1,Phi,mc=40,burn=10)

```

## Details

For more information on QFSTS Package, please access the package documentations or [vignettes](https://github.com/LauraHuang0919/QFSTS/tree/main/vignettes). Please feel free to contact the author.
