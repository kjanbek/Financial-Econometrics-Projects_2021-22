#Loading libraries
library(sandwich)
library(lmtest)
library(MASS)
library(zoo)
library(strucchange)
library(urca)
library(tseries)
library(xtable)
library(tsDyn)
library(fGarch)
library(forecast)
library(pracma)
library(xts)
library(seastests) 
library(VAR.etp)
library(vars)
library(fArma)
library(tidyverse)
library(mFilter)

#Data import
setwd("C:/~/FE - EMP APP 1")
data <-read.table("FRED_data.txt", fill=TRUE, header=TRUE, dec=".")

#Declare variables as time series
Dates <- as.Date(data[,1])
SP500 <- ts(data$SP500, frequency=12,start=c(2011,10))
WTI <- ts(data$Oil_prices, frequency=12,start=c(2011,10))
Breakeven5Y <- ts(data$X5.year_breakeven_inflation_rate, frequency=12,start=c(2011,10))
FRED_data <- ts(data[,-1],start=c(2011,10),frequency=12)

#Plot the time series
plot.ts(FRED_data,plot.type="multiple")
plot(Dates,SP500,type="l")
plot(Dates,WTI,type="l")
plot(Dates,Breakeven5Y,type="l")

    ############################################
    ## APPLICATION 1 - FINANCIAL ECONOMETRICS ##
    ############################################

#Log transformation of time series (except for the 5Y breakeven infl., expressed in %)
lSP500 <-log(SP500)
lWTI <-log(WTI)

plot(Dates,lSP500,type="l", main = "log(SP500)")  
plot(Dates,lWTI,type="l", main = "log WTI")

# Decomposition of time series (to get a first glance)
decomp_lSP500 <- decompose(lSP500)
decomp_lWTI <- decompose(lWTI)
decomp_Breakeven5Y <- decompose(Breakeven5Y) 

plot(decomp_lSP500)
plot(decomp_lWTI)
plot(decomp_Breakeven5Y)

#Testing for a deterministic trend (and drift)
tsmodel_lSP500 <- tslm(lSP500~trend)
summary(tsmodel_lSP500)

tsmodel_lWTI <- tslm(lWTI~trend)
summary(tsmodel_lWTI)

tsmodel_Breakeven5Y <- tslm(Breakeven5Y~trend)
summary(tsmodel_Breakeven5Y)

#Testing seasonality (ACF/PACF and WO-test):
Acf(Breakeven5Y, lag=36)
Acf(lSP500, lag=36)
Acf(lWTI, lag=36)

pacf(Breakeven5Y, lag=36)
pacf(lSP500)
pacf(lWTI, lag=36)

  #WO tests
isSeasonal(Breakeven5Y)
summary(wo(Breakeven5Y))

isSeasonal(lSP500)
summary(wo(lSP500))

isSeasonal(lWTI)
summary(wo(lWTI))

#Removing the time trend from lSP500 and lWTI:
lSP500_num <- as.numeric(lSP500)
lSP500_detrended <- detrend(lSP500_num,'linear')
lSP500_detrended_ts <- ts(lSP500_detrended, frequency=12, start=c(2011,10))
plot(Dates, lSP500_detrended_ts, type="l", main = "De-trended log(S&P500)")


lWTI_num <- as.numeric(lWTI)
lWTI_detrended <- detrend(lWTI_num,'linear')
lWTI_detrended_ts <- ts(lWTI_detrended, frequency=12, start=c(2011,10))
plot(Dates, lWTI_detrended_ts, type="l", main = "De-trended log(WTI)")


#Testing for a seasonal variations (see above)
tsmodel_SP500 <- tslm(lSP500_detrended_ts~season)
summary(tsmodel_SP500)

tsmodel_WTI <- tslm(lWTI_detrended_ts~season)
summary(tsmodel_WTI)

tsmodel_Breakeven5Y <- tslm(Breakeven5Y~trend + season)
summary(tsmodel_Breakeven5Y)

#Testing for a stochastic trend through a succinct ADF/PP tests

  #ADF
adf.test(Breakeven5Y, k=1)
adf.test(lSP500_detrended_ts, k=1)
adf.test(lWTI_detrended_ts, k=1)

adf.test(Breakeven5Y_DIFF_ts, k=1)
adf.test(lSP500_detrended_DIFF_ts, k=1)
adf.test(lWTI_detrended_DIFF_ts, k=1)

  #PP
pp.test(Breakeven5Y)
pp.test(lSP500_detrended_ts)
pp.test(lWTI_detrended_ts)

pp.test(Breakeven5Y_DIFF_ts)
pp.test(lSP500_detrended_DIFF_ts)
pp.test(lWTI_detrended_DIFF_ts)


#Testing for a stochastic trend through an extended ADF test
Breakeven5Y_ADF <- ur.df(Breakeven5Y,type="trend",selectlags="AIC")
summary(Breakeven5Y_ADF)

SP500_ADF <- ur.df(lSP500_detrended_ts,type="trend",selectlags="AIC")
summary(SP500_ADF)

WTI_ADF <- ur.df(lWTI_detrended_ts,type="trend",selectlags="AIC")
summary(WTI_ADF)

#First-differencing the series, to make them stationary
lSP500_detrended_DIFF <- diff(lSP500_detrended_ts, lag=1, differences=1)
lWTI_detrended_DIFF <- diff(lWTI_detrended_ts, lag=1, differences=1)
Breakeven5Y_DIFF <- diff(Breakeven5Y, lag=1, differences=1)

lSP500_detrended_DIFF_ts <- ts(lSP500_detrended_DIFF, frequency=12, start=c(2011,10))
lWTI_detrended_DIFF_ts <- ts(lWTI_detrended_DIFF, frequency=12, start=c(2011,10))
Breakeven5Y_DIFF_ts <- ts(Breakeven5Y_DIFF, frequency=12, start=c(2011,10))

plot(lSP500_detrended_DIFF_ts, type="l", main = "log(S&P500) stationarized")
plot(lWTI_detrended_DIFF_ts, type="l", main = "log(WTI) stationarized")
plot(Breakeven5Y_DIFF_ts, type="l", main = "5-year breakeven stationarized")

#Checking if the first-differenced series are now stationary
SP500_DIFF_ADF <- ur.df(lSP500_detrended_DIFF_ts,type="trend",selectlags="AIC")
summary(SP500_DIFF_ADF)

WTI_DIFF_ADF <- ur.df(lWTI_detrended_DIFF_ts,type="trend",selectlags="AIC")
summary(WTI_DIFF_ADF)

Breakeven5Y_DIFF_ADF <- ur.df(Breakeven5Y_DIFF_ts,type="trend",selectlags="AIC")
summary(Breakeven5Y_DIFF_ADF)


#Estimating the cyclical component with an ARMA model(Using de-trended and first-differenced Yt)
  # Step 1: get an auto fit for p,d,q components for each ARMA model
auto.arima(Breakeven5Y_DIFF_ts, seasonal=FALSE)
auto.arima(lSP500_detrended_DIFF_ts, seasonal=FALSE)
auto.arima(lWTI_detrended_DIFF_ts, seasonal=FALSE)


  # Step 2: we fit an ARMA for each series comparing the optimal model given by the auto.arima function and the benchmark (default) model

# Optimal model (auto.arima)
Breakeven5Y_ARMA <- arima(Breakeven5Y_DIFF_ts, order=c(2,0,2))
summary(Breakeven5Y_ARMA)
tsdisplay(residuals(Breakeven5Y_ARMA),lag.max = 40)
# Benchmark model (default)
Breakeven5Y_ARMA_bench <- arima(Breakeven5Y_DIFF_ts, order=c(1,0,1))
summary(Breakeven5Y_ARMA_bench)
tsdisplay(residuals(Breakeven5Y_ARMA_bench),lag.max = 40)


# Optimal model (auto.arima)
SP500_ARMA <- arima(lSP500_detrended_DIFF_ts, order=c(2,0,2))
summary(SP500_ARMA)
tsdisplay(residuals(SP500_ARMA),lag.max = 40)
# Benchmark model (default)
SP500_ARMA_bench <- arima(lSP500_detrended_DIFF_ts, order=c(1,0,1))
summary(SP500_ARMA_bench)
tsdisplay(residuals(SP500_ARMA_bench),lag.max = 40)


# Optimal model (auto.arima)
WTI_ARMA <- arima(lWTI_detrended_DIFF_ts, order=c(0,0,1))
summary(WTI_ARMA)
tsdisplay(residuals(WTI_ARMA),lag.max = 40)
# Benchmark model (default)
WTI_ARMA_bench <- arima(lWTI_detrended_DIFF_ts, order=c(1,0,1))
summary(WTI_ARMA_bench)
tsdisplay(residuals(WTI_ARMA_bench),lag.max = 40)


  #Step 3 We test for the normality of the residuals using a Breusch-Godfrey test

bgtest(Breakeven5Y_DIFF_ts~lSP500_detrended_DIFF_ts+lSP500_detrended_DIFF_ts, data = data)

bgtest(lSP500_detrended_DIFF_ts~Breakeven5Y_DIFF_ts+lWTI_detrended_DIFF_ts, data = data)

bgtest(lWTI_detrended_DIFF_ts~lSP500_detrended_DIFF_ts+Breakeven5Y_DIFF_ts, data = data)

#Saving the Empirical Application 1 dataset
save.image("FE_App_1.RData")

    ############################################
    ## APPLICATION 2 - FINANCIAL ECONOMETRICS ##
    ############################################

#Loading libraries
library(sandwich)
library(lmtest)
library(MASS)
library(zoo)
library(strucchange)
library(urca)
library(tseries)
library(xtable)
library(tsDyn)
library(fGarch)
library(forecast)
library(pracma)
library(xts)
library(seastests) 
library(VAR.etp)
library(vars)
library(fArma)
library(tidyverse)
library(mFilter)


#Data import
setwd("C:/~/FE - EMP APP 1")
load(file = 'FE_App_1.RData' ) 

#Estimating a VAR model 
var1 <- cbind(Breakeven5Y_DIFF_ts, lSP500_detrended_DIFF_ts,lWTI_detrended_DIFF_ts)
colnames(var1) <- cbind("Breakeven5Y", "SP500", "WTI")

#Let's determine the lag order
lagselect <- VARselect(var1, lag.max = 10, type = "const")
lagselect$selection

#Optimal lag is 2. Let's now estimate our model :
Model1 <- VAR(var1, p = 2, type = "const", season = NULL, exog = NULL) 
summary(Model1)

#We check for the robustness of our model though many tests 
  
  # Serial correlation
Serial1 <- serial.test(Model1, lags.pt = 12, type = "PT.asymptotic")
Serial1

  # Heteroscedasticity
arch1 <- arch.test(Model1, lags.multi = 12, multivariate.only = T)
arch1

  # Normal distribution of residuals
shapiro.test(residuals(Model1))
hist(residuals(Model1))
boxplot(residuals(Model1))


  # Model stability - Structural breaks
struct1 <- stability(Model1, type = "OLS-CUSUM")
plot(struct1)


# Forecasting each variable at horizon 2
forecastVAR <- predict(Model1, n.ahead = 2, ci = 0.95) 
fcast_model <- forecast(Model1, h=2)
plot(fcast_model)
fanchart(forecastVAR, names = 'Breakeven5Y', main = 'Forecast of the Inflation Expectations variable at horizon 2')
fanchart(forecastVAR, names = 'SP500', main = 'Forecast of the S&P500 Stock Index variable at horizon 2')
fanchart(forecastVAR, names = 'WTI', main = 'Forecast of the WTI oil prices variable at horizon 2')


# ARIMA forecast for benchmark
forecastARMA_BRK <- forecast(Breakeven5Y_ARMA, h=5)
plot(forecastARMA_BRK)

forecastARMA_SP500 <- forecast(SP500_ARMA, h=5)
plot(forecastARMA_SP500)

forecastARMA_WTI <- forecast(WTI_ARMA, h=5)
plot(forecastARMA_WTI)


#Granger Causality tests
GrangerBRK <- causality(Model1, cause = 'Breakeven5Y')
GrangerBRK$Granger

GrangerSP500 <- causality(Model1, cause = 'SP500')
GrangerSP500$Granger

GrangerWTI <- causality(Model1, cause = 'WTI')
GrangerWTI$Granger

#Impulse Reaction Function of the three variables

  # Inflation Expectations Shock
irfBRK_SP <- irf(Model1, impulse = 'Breakeven5Y', response = 'SP500', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfBRK_SP, ylab ='S&P500 Stock Index', main = 'Shock on Inflation Expectations' )

irfBRK_WTI <- irf(Model1, impulse = 'Breakeven5Y', response = 'WTI', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfBRK_WTI, ylab ='WTI oil prices', main = 'Shock on Inflation Expectations' )


# S&P500 Shock
irfSP_BRK <- irf(Model1, impulse = 'SP500', response = 'Breakeven5Y', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfSP_BRK, ylab ='5 years Inflation Expectations', main = 'Shock on the S&P500 Stock Index' )

irfSP_WTI <- irf(Model1, impulse = 'SP500', response = '', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfSP_WTI, ylab ='WTI oil prices', main = 'Shock on the S&P500 Stock Index' )


# WTI Shock
irfWTI_BRK <- irf(Model1, impulse = 'WTI', response = 'Breakeven5Y', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfWTI_BRK, ylab ='5 years Inflation Expectations', main = 'Shock on WTI oil prices' )

irfWTI_SP <- irf(Model1, impulse = 'WTI', response = 'SP500', n.ahead = 6, boot = T, runs = 200, ci = 0.95)
plot(irfWTI_SP, ylab ='S&P500 Stock Index', main = 'Shock on WTI oil prices' )


#Save data from Application 2
save.image("FE_App_2.RData")

############################################
## APPLICATION 3 - FINANCIAL ECONOMETRICS ##
############################################





