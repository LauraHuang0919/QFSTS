# QFSTS: Quantile Feature Selection Time Series Model

Quantile feature selection over correlated multivariate time series data has always been a methodological challenge and is an open problem. In this package, we propose a general Bayesian dimension reduction methodology for feature selection in high-dimensional joint quantile time series analysis, under the name of the quantile feature selection time series (QFSTS) model. This QFSTS model is an extension from Dr. Ning's previous MBSTS package with similar functionality but adding more settings on the quantile feature selection. The QFSTS model is a general structural time series model, where each component yields an additive contribution to the time series modeling with direct interpretations. Feature selection is conducted in the quantile regression component, where each time series has its own pool of contemporaneous external predictors allowing nowcasting.

## Installation 

devtools::install_github("LauraHuang0919/QFSTS")
