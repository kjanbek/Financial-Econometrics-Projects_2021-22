
    ##################################################
    ## FINANCIAL ECONOMETRIC EMPIRICAL APPLICATIONS ##
    ##################################################
    
  # Authors : 
      # BOUILLOT Roland 
      # JANBEK Khalil
      # LOUAFI Mehdi
    
  # Summary 
      # Empirical Application 1 - Line 18
      # Empirical Application 2 - Line 233
      # Empirical Application 3 - Line 356
      # Empirical Application 4 - Line 640
    
    
    ############################################
    ## APPLICATION 1 - FINANCIAL ECONOMETRICS ##
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
library(pastecs)
library(mFilter)

    # PART I : Fitting the data
    

#Data import
setwd("C:/~/FE - EMP APP 1 - Copie")
data <-read.table("FRED_Data_2.txt", fill=TRUE, header=TRUE, dec=".")

#Declare variables as time series
Dates <- as.Date(data[,1])
SP500 <- ts(data$SP500, frequency=12,start = c(2011, 11))
WTI <- ts(data$WTI, frequency=12,start=c(2011,11))
CPI <- ts(data$CPI, frequency=12,start=c(2011,11))
IPROD <- ts(data$IPROD, frequency = 12, start = c(2011,11))
FRED_data <- ts(data[,-1],start=c(2011,11),frequency=12)

# Log-transforming the variables 
lSP500 <-log(SP500)
lWTI <-log(WTI)
lCPI <-log(CPI)
lIPROD <-log(IPROD)
log_FRED_data <- cbind(lSP500, lWTI, lCPI, lIPROD)

#Descriptive statistics of the variables
summary(FRED_data)
StatDesc <-stat.desc(FRED_data)
round(StatDesc, 1)

summary(log_FRED_data)
log_StatDesc <-stat.desc(log_FRED_data)
round(log_StatDesc, 2)

#Plot the time series
plot.ts(FRED_data,plot.type="multiple")
plot(Dates,SP500,type="l", main = "SP500")
plot(Dates,WTI,type="l", main = "WTI")
plot(Dates,CPI,type="l", main = "CPI")
plot(Dates,IPROD,type="l", main = "IPROD")

plot.ts(log_FRED_data,plot.type="multiple")
plot(Dates,lSP500,type="l", main = "log(SP500)")
plot(Dates,lWTI,type="l", main = "log(WTI)")
plot(Dates,lCPI,type="l", main = "log(CPI)")
plot(Dates,lIPROD,type="l", main = "log(IPROD)")

# Decomposition of the time series
decomp_lSP500 <- decompose(lSP500)
decomp_lWTI <- decompose(lWTI)
decomp_lCPI <- decompose(lCPI)
decomp_lIPROD <- decompose(lIPROD) 

plot(decomp_lSP500)
plot(decomp_lWTI)
plot(decomp_lCPI)
plot(decomp_lIPROD)

     # i. Deterministic Trend

#Testing for a deterministic trend (and drift)
trend_lSP500 <- tslm(lSP500~trend)
summary(trend_lSP500)

trend_lWTI <- tslm(lWTI~trend)
summary(trend_lWTI)

trend_lCPI <- tslm(lCPI~trend)
summary(trend_lCPI)

trend_lIPROD <- tslm(lIPROD~trend)
summary(trend_lIPROD)


    # ii. Seasonal Component

#Testing seasonality (ACF/PACF and WO-test):
Acf(lSP500, lag=36)
Acf(lCPI, lag=36)
Acf(lWTI, lag=36)
Acf(lIPROD, lag=36)

pacf(lSP500, lag=36)
pacf(lCPI, lag=36)
pacf(lWTI, lag=36)
pacf(lIPROD, lag=36)

#WO seasonality tests
isSeasonal(lSP500)
summary(wo(lSP500))

isSeasonal(lCPI)
summary(wo(lCPI))

isSeasonal(lWTI)
summary(wo(lWTI))

isSeasonal(lIPROD)
summary(wo(lIPROD))


#Removing the time trend from lSP500, lWTI and lCPI:
lSP500_num <- as.numeric(lSP500)
lSP500_detrended <- detrend(lSP500_num,'linear')
lSP500_detrended_ts <- ts(lSP500_detrended, frequency=12, start=c(2011,11))
plot(Dates, lSP500_detrended_ts, type="l", main = "De-trended log(S&P500)")


lWTI_num <- as.numeric(lWTI)
lWTI_detrended <- detrend(lWTI_num,'linear')
lWTI_detrended_ts <- ts(lWTI_detrended, frequency=12, start=c(2011,11))
plot(Dates, lWTI_detrended_ts, type="l", main = "De-trended log(WTI)")

lCPI_num <- as.numeric(lCPI)
lCPI_detrended <- detrend(lCPI_num,'linear')
lCPI_detrended_ts <- ts(lCPI_detrended, frequency=12, start=c(2011,11))
plot(Dates, lCPI_detrended_ts, type="l", main = "De-trended log(CPI)")

#First-differencing the series, to make them stationary
lSP500_DIFF <- diff(lSP500, lag=1, differences=1)
lWTI_DIFF <- diff(lWTI, lag=1, differences=1)
lCPI_DIFF <- diff(lCPI, lag=1, differences=1)
lIPROD_DIFF <- diff(lIPROD, lag=1, differences=1)

lSP500_DIFF_ts <- ts(lSP500_DIFF, frequency=12, start=c(2011,11))
lWTI_DIFF_ts <- ts(lWTI_DIFF, frequency=12, start=c(2011,11))
lCPI_DIFF_ts <- ts(lCPI_DIFF, frequency=12, start=c(2011,11))
lIPROD_DIFF_ts <- ts(lIPROD_DIFF, frequency=12, start=c(2011,11))

plot(lSP500_DIFF_ts, type="l", main = "log(S&P500) stationarized")
plot(lWTI_DIFF_ts, type="l", main = "log(WTI) stationarized")
plot(lCPI_DIFF_ts, type="l", main = "log(CPI) stationarized")
plot(lIPROD_DIFF_ts, type="l", main = "log(IPROD) stationarized")

#Testing for stationarity through ADF/PP unit root tests

#ADF
adf.test(lSP500, k=1)
adf.test(lWTI, k=1)
adf.test(lCPI, k=1)
adf.test(lIPROD, k=1)

adf.test(lSP500_DIFF_ts, k=1)
adf.test(lWTI_DIFF_ts, k=1)
adf.test(lCPI_DIFF_ts, k=1)
adf.test(lIPROD_DIFF_ts, k=1)

#PP
pp.test(lSP500)
pp.test(lWTI)
pp.test(lCPI)
pp.test(lIPROD)

pp.test(lSP500_DIFF_ts)
pp.test(lWTI_DIFF_ts)
pp.test(lCPI_DIFF_ts)
pp.test(lIPROD_DIFF_ts)

# More detailed ADF tests controlling for trend
ur_lSP500_trend <- ur.df(lSP500, type =  "trend", selectlags =  "AIC")
summary(ur_lSP500_trend)

ur_lWTI_trend <- ur.df(lWTI, type =  "trend", selectlags =  "AIC")
summary(ur_lWTI_trend)

ur_lCPI_trend <- ur.df(lCPI, type =  "trend", selectlags =  "AIC")
summary(ur_lCPI_trend)

ur_lIPROD_trend <- ur.df(lIPROD, type =  "trend", selectlags =  "AIC")
summary(ur_lIPROD_trend)

    #First Differenced variables
ur_lSP500_DIFF_trend <- ur.df(lSP500_DIFF_ts, type =  "trend", selectlags =  "AIC")
summary(ur_lSP500_DIFF_trend)

ur_lWTI_DIFF_trend <- ur.df(lWTI_DIFF_ts, type =  "trend", selectlags =  "AIC")
summary(ur_lWTI_DIFF_trend)

ur_lCPI_DIFF_trend <- ur.df(lCPI_DIFF_ts, type =  "trend", selectlags =  "AIC")
summary(ur_lCPI_DIFF_trend)

ur_lIPROD_DIFF_trend <- ur.df(lIPROD_DIFF_ts, type =  "trend", selectlags =  "AIC")
summary(ur_lIPROD_DIFF_trend)

model1 <- cbind(lSP500, lWTI, lCPI, lIPROD)
model2 <- cbind(lSP500, lWTI, lCPI)
model3 <- cbind(lSP500, lWTI, lIPROD)
model4 <- cbind(lSP500, lCPI, lIPROD)


#Saving Dataset for Application 3
save.image("FE_App_3.RData")


    # PART II : Estimating our model
    

# We load the data
load(file = 'FE_App_3.RData' ) 

# We determine the optimal number of lags 
optimal_lag <- VARselect(model1, lag.max = 10, type = "const")
optimal_lag$selection

Trace <- ca.jo(model1, type = "trace", ecdet = "const", K = 2, spec="longrun")
summary(Trace)

Eigen <- ca.jo(model1, type = "eigen", ecdet = "const", K = 2, spec="longrun")
summary(Eigen)

# Building the VEC model
VEC1 <- VECM(model1, 2, r=1, include =("const"), estim=("2OLS"))
summary(VEC1)

# Robustness checks for our VEC model

  # VEC to VAR models 
VAR1 <- vec2var(Trace, r=1)

  # Serial Correlation
SerCorr <- serial.test(VAR1, lags.pt = 2, type = "PT.asymptotic")
SerCorr

  #ARCH effects 
Arch <- arch.test(VAR1, lags.multi = 15, multivariate.only = TRUE)
Arch

  #normality
Normal <- normality.test(VAR1, multivariate.only = TRUE)
Normal
?normality.test

#Impulse reaction functions

CPI_irf <- irf(VAR1, impulse = "lCPI", response="lSP500", n.ahead = 20, boot=TRUE, runs = 500)
plot(CPI_irf, ylab ="SP500", main="CPI shock to SP500")

WTI_irf <- irf(VAR1, impulse = "lWTI", response="lSP500", n.ahead = 20, boot=TRUE, runs = 500)
plot(WTI_irf, ylab ="SP500", main="WTI shock to SP500")

IPROD_irf <- irf(VAR1, impulse = "lIPROD", response="lSP500", n.ahead = 20, boot=TRUE, runs = 500)
plot(IPROD_irf, ylab ="SP500", main="IPROD shock to SP500")

# Variance Decomposition
VarDecomp <- fevd(VAR1, n.ahead = 10)
plot(VarDecomp)

#Saving the final Dataset for Application 3
save.image("FE_App_3.RData")

    # PART III : Forecasting our model 


# VEC Model Forecasting
VECM_forecast <- predict(VAR1, n.ahead=10, ci=0.95)
fanchart(VECM_forecast, names="lSP500", main ="Forecast S&P500 Index", xlab="Horizon", ylab="lSP500")
fanchart(VECM_forecast, names="lWTI", main ="Forecast WTI Oil prices", xlab="Horizon", ylab="lWTI")
fanchart(VECM_forecast, names="lCPI", main ="Forecast CPI", xlab="Horizon", ylab="lCPI")
fanchart(VECM_forecast, names="lIPROD", main ="Forecast Industrial Production Index", xlab="Horizon", ylab="lIPROD")


    ############################################
    ## APPLICATION 4 - FINANCIAL ECONOMETRICS ##
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
library(rugarch)
library(quantmod)
library(rmgarch)
library(PerformanceAnalytics)

#___________________________
# /!\ MONTHLY DATA MODEL /!\
#___________________________

# Load the data
load(file = 'FE_App_3.RData') 


# S&P500 monthly returns
Returns_SP500 = diff(log(SP500))
plot(Returns_SP500, type="l")
hist(Returns_SP500)


# PART I : FITTING THE DATA
#__________________________

  # i. Decomposition of the time series
decomp_Returns_SP500 <- decompose(Returns_SP500)
plot(decomp_Returns_SP500)


  # ii. Deterministic Trend and Drift
trend_Returns_SP500 <- tslm(Returns_SP500~trend)
summary(trend_Returns_SP500)

  # iii. Seasonal Component

#ACF/PACF
Acf(Returns_SP500, lag=36)
pacf(Returns_SP500, lag=36)

#WO seasonality test
isSeasonal(Returns_SP500)
summary(wo(Returns_SP500))

  # iv. Stochastic trend (unit root test)

#ADF
adf.test(Returns_SP500, k=1)

#PP
pp.test(Returns_SP500)


# PART II : Standard GARCH MODEL (1,1)
#_____________________________________

  # Model 1 : ARIMA (2,0,2), sGARCH (1,1)

      #Estimation
      #__________

# We define the optimal (p,q) parameters for the ARIMA model specification
auto.arima(Returns_SP500, seasonal=FALSE, stationary= TRUE)

# An ARIMA (2,0,2) is the optimal fit for our standard GARCH specification
Spec1 <- ugarchspec(mean.model = list(armaOrder = c(2,2)),
                variance.model = list(model = "sGARCH"),
                distribution.model = 'norm')

Spec1 <- ugarchspec(mean.model = list(armaOrder = c(1,1), archm = TRUE),
                    variance.model = list(model = "sGARCH", garchOrder=c(1,0)),
                    distribution.model = 'norm')
Spec1

# We estimate the sGARCH model
Garch8_est <- ugarchfit(data = Returns_SP500, spec = Spec1)
Garch8_est
Garch7_est@fit$coef
plot(Garch1_est)

Garch1_est_var <- Garch1_est@fit$var
Garch1_est_res <- (Garch1_est@fit$residuals)^2
plot(Garch1_est_res, type = "l")
lines(Garch1_est_var, col = "red")

      #Forecast
      #________

#We forecast our GARCH model over a 20 periods horizon
Garch1_for <- ugarchforecast(fitORspec = Garch1_est, n.ahead = 20)
Garch1_for
plot(fitted(Garch1_for))
plot(sigma(Garch1_for))

Garch1_for_sig <- Garch1_for@forecast$sigmaFor
plot(Garch1_for_sig, type = "l")

# For a better display of the forecast's results, we only take the last 20 observations 
Garch1_for_var_tail <- c(tail(Garch1_est_var,20),rep(NA,10))
Garch1_for_res_tail <- c(tail(Garch1_est_res,20),rep(NA,10))
Garch1_forecast <- c(rep(NA,20),(Garch1_for_sig)^2)

plot(Garch1_for_res_tail, type = "l")
lines(Garch1_for_var_tail, col = "red")
lines(Garch1_forecast, col = "green")


#-----------------------------------------------------------------
#-----------------------------------------------------------------


#___________________________
# /!\ DAILY DATA MODEL /!\
#___________________________


# PART I : FITTING THE DATA
#__________________________

# Loading the S&P500, Google and Microsoft daily prices series
getSymbols("^GSPC",from = "2011-11-01",to = "2021-08-01")
getSymbols("GOOGL",from = "2011-11-01",to = "2021-08-01")
getSymbols("MSFT",from = "2011-11-01",to = "2021-08-01")

# We plot the time series
plot(GSPC$GSPC.Close)
plot(GOOGL$GOOGL.Close)
plot(MSFT$MSFT.Close)

# Daily returns on close prices 
    # S&P500
SPreturns <- CalculateReturns(GSPC$GSPC.Close)
SPreturns <- ts(CalculateReturns(GSPC$GSPC.Close)[-1], frequency=252,start=c(2011,11,1))
hist(SPreturns)
plot(SPreturns, type = 'l')

    # Google
GOOGreturns <- CalculateReturns(GOOGL$GOOGL.Close)
GOOGreturns <- ts(CalculateReturns(GOOGL$GOOGL.Close)[-1], frequency=252,start=c(2011,11,1))
hist(GOOGreturns)
plot(GOOGreturns, type = 'l')

    # Microsoft
MSFTreturns <- CalculateReturns(MSFT$MSFT.Close)
MSFTreturns <- ts(CalculateReturns(MSFT$MSFT.Close)[-1], frequency=252,start=c(2011,11,1))
hist(MSFTreturns)
plot(MSFTreturns, type = 'l')

# We gather the Google and Microsoft into a dataframe that will be later used for the multivariate analysis
Techreturns <- data.frame(GOOGreturns, MSFTreturns)
names(Techreturns)[1] <- "Google returns"
names(Techreturns)[2] <- "Microsoft returns"

# Deterministic trends
    # i. Decomposition of the time series
decomp_SPreturns <- decompose(SPreturns)
plot(decomp_SPreturns)

decomp_GOOGreturns <- decompose(GOOGreturns)
plot(decomp_GOOGreturns)

decomp_MSFTreturns <- decompose(MSFTreturns)
plot(decomp_MSFTreturns)

    # ii. Deterministic Trend and Drift
trend_SPreturns <- tslm(SPreturns~trend)
summary(trend_SPreturns)

trend_GOOGreturns <- tslm(GOOGreturns~trend)
summary(trend_GOOGreturns)

trend_MSFTreturns <- tslm(MSFTreturns~trend)
summary(trend_MSFTreturns)

    # iii. Seasonal Component

# ACF/PACF
Acf(na.omit(SPreturns))
pacf(na.omit(SPreturns))

Acf(na.omit(GOOGreturns))
pacf(na.omit(GOOGreturns))

Acf(na.omit(MSFTreturns))
pacf(na.omit(MSFTreturns))

#WO seasonality test
isSeasonal(SPreturns)
summary(wo(SPreturns))

isSeasonal(GOOGreturns)
summary(wo(GOOGreturns))

isSeasonal(MSFTreturns)
summary(wo(MSFTreturns))

    # iv. Stochastic trend (unit root test)

#ADF
adf.test(na.omit(SPreturns, k=1))
adf.test(na.omit(GOOGreturns, k=1))
adf.test(na.omit(MSFTreturns, k=1))

#PP
pp.test(na.omit(SPreturns))
pp.test(na.omit(GOOGreturns))
pp.test(na.omit(MSFTreturns))

#We check for ARCH effects in our data sample
SPreturnsArchTest <- ArchTest(SPreturns, lags=1, demean=TRUE)
SPreturnsArchTest

GOOGreturnsArchTest <- ArchTest(GOOGreturns, lags=1, demean=TRUE)
GOOGreturnsArchTest

MSFTreturnsArchTest <- ArchTest(MSFTreturns, lags=1, demean=TRUE)
MSFTreturnsArchTest


# PART II : Standard GARCH MODEL (1,1)
#_____________________________________

  # Model 2 : SP TIME SERIE  -  ARIMA (4,0,5), sGARCH (1,1)

      #Estimation
      #__________

# We define the optimal (p,q) parameters for the ARIMA model specification
auto.arima(SPreturns, seasonal=FALSE, stationary = TRUE)

# An ARIMA (4,0,5) is the optimal fit for our standard GARCH specification
Spec2 <- ugarchspec(mean.model = list(armaOrder = c(4,5)),
                    variance.model = list(model = "sGARCH"),
                    distribution.model = 'norm')
Spec2

# We estimate the sGARCH model
Garch2_est <- ugarchfit(data = na.omit(SPreturns), spec = Spec2)
Garch2_est
Garch2_est@fit$coef
plot(Garch2_est, which = 'all')

Garch2_est_var <- Garch2_est@fit$var
Garch2_est_res <- (Garch2_est@fit$residuals)^2
plot(Garch2_est_res, type = "l")
lines(Garch2_est_var, col = "red")

      #Forecast
      #________

#We forecast our GARCH model over a 20 periods horizon
Garch2_for <- ugarchforecast(fitORspec = Garch2_est, n.ahead = 20)
Garch2_for
plot(fitted(Garch2_for))
plot(sigma(Garch2_for))

# We forecast ou sigma (that is the h-hat for the conditional variance)
Garch2_for_sig <- Garch2_for@forecast$sigmaFor
plot(Garch2_for_sig, type = "l")

# For a better display of the forecast's results, we only take the last 50 observations 
Garch2_for_var_tail <- c(tail(Garch2_est_var,50),rep(NA,10))
Garch2_for_res_tail <- c(tail(Garch2_est_res,50),rep(NA,10))
Garch2_forecast <- c(rep(NA,50),(Garch2_for_sig)^2)

plot(Garch2_for_res_tail, type = "l")
lines(Garch2_for_var_tail, col = "red")
lines(Garch2_forecast, col = "green")

    # Model 3 : GOOGLE TIME SERIE  -  ARIMA (0,0,1), sGARCH (1,1)

#Estimation
#__________

# We define the optimal (p,q) parameters for the ARIMA model specification
auto.arima(GOOGreturns, seasonal=FALSE, stationary = TRUE)

# An ARIMA (4,0,5) is the optimal fit for our standard GARCH specification
Spec3 <- ugarchspec(mean.model = list(armaOrder = c(0,1)),
                    variance.model = list(model = "sGARCH"),
                    distribution.model = 'norm')
Spec3

# We estimate the sGARCH model
Garch3_est <- ugarchfit(data = na.omit(GOOGreturns), spec = Spec3)
Garch3_est
Garch3_est@fit$coef
plot(Garch3_est, which = 'all')

Garch3_est_var <- Garch3_est@fit$var
Garch3_est_res <- (Garch3_est@fit$residuals)^2
plot(Garch3_est_res, type = "l")
lines(Garch3_est_var, col = "red")

#Forecast
#________

#We forecast our GARCH model over a 20 periods horizon
Garch3_for <- ugarchforecast(fitORspec = Garch3_est, n.ahead = 20)
Garch3_for
plot(fitted(Garch3_for))
plot(sigma(Garch3_for))

# We forecast ou sigma (that is the h-hat for the conditional variance)
Garch3_for_sig <- Garch3_for@forecast$sigmaFor
plot(Garch3_for_sig, type = "l")

# For a better display of the forecast's results, we only take the last 50 observations 
Garch3_for_var_tail <- c(tail(Garch3_est_var,50),rep(NA,10))
Garch3_for_res_tail <- c(tail(Garch3_est_res,50),rep(NA,10))
Garch3_forecast <- c(rep(NA,50),(Garch3_for_sig)^2)

plot(Garch3_for_res_tail, type = "l")
lines(Garch3_for_var_tail, col = "red")
lines(Garch3_forecast, col = "green")

    # Model 4 : MICROSOFT TIME SERIE  -  ARIMA (4,0,4), sGARCH (1,1)

#Estimation
#__________

# We define the optimal (p,q) parameters for the ARIMA model specification
auto.arima(MSFTreturns, seasonal=FALSE, stationary = TRUE)

# An ARIMA (4,0,5) is the optimal fit for our standard GARCH specification
Spec4 <- ugarchspec(mean.model = list(armaOrder = c(4,4)),
                    variance.model = list(model = "sGARCH"),
                    distribution.model = 'norm')
Spec4

# We estimate the sGARCH model
Garch4_est <- ugarchfit(data = na.omit(MSFTreturns), spec = Spec4)
Garch4_est
Garch4_est@fit$coef
plot(Garch4_est, which = 'all')

Garch4_est_var <- Garch4_est@fit$var
Garch4_est_res <- (Garch4_est@fit$residuals)^2
plot(Garch4_est_res, type = "l")
lines(Garch4_est_var, col = "red")

#Forecast
#________

#We forecast our GARCH model over a 20 periods horizon
Garch4_for <- ugarchforecast(fitORspec = Garch4_est, n.ahead = 20)
Garch4_for
plot(fitted(Garch4_for))
plot(sigma(Garch4_for))

# We forecast ou sigma (that is the h-hat for the conditional variance)
Garch4_for_sig <- Garch4_for@forecast$sigmaFor
plot(Garch4_for_sig, type = "l")

# For a better display of the forecast's results, we only take the last 50 observations 
Garch4_for_var_tail <- c(tail(Garch4_est_var,50),rep(NA,10))
Garch4_for_res_tail <- c(tail(Garch4_est_res,50),rep(NA,10))
Garch4_forecast <- c(rep(NA,50),(Garch4_for_sig)^2)

plot(Garch4_for_res_tail, type = "l")
lines(Garch4_for_var_tail, col = "red")
lines(Garch4_forecast, col = "green")


          # PART III : ARCH-M model
          #________________________

    # Model 5 : MICROSOFT ARIMA (4,0,4), ARCH-M (1)

# We define the optimal (p,q) parameters for the ARIMA model specification
auto.arima(MSFTreturns, seasonal=FALSE, stationary = TRUE)

# An ARIMA (4,0,4) is the optimal fit for our ARCH-M specification
# We use a ARCH-M option as well as an ARCH model by defining the GARCH order as (1,0)
Spec6 <- ugarchspec(mean.model = list(armaOrder = c(4,4), archm = TRUE),
                    variance.model = list(model = "sGARCH", garchOrder=c(1,0)),
                    distribution.model = 'norm')
Spec6

# We estimate the ARCH-M model
Archm6_est <- ugarchfit(data = MSFTreturns, spec = Spec6, solver = 'hybrid')
Archm6_est
Archm6_est@fit$coef
plot(Archm6_est, which = 'all')


Archm6_est_var <- Archm6_est@fit$var
Archm6_est_res <- (Archm6_est@fit$residuals)^2
plot(Archm6_est_res, type = "l")
lines(Archm6_est_var, col = "red")

#Forecast
#________

#We forecast our ARCH-M model over a 20 periods horizon
Archm6_for <- ugarchforecast(fitORspec = Archm6_est, n.ahead = 20)
Archm6_for
plot(fitted(Archm6_for))
plot(sigma(Archm6_for))

# We forecast ou sigma (that is the h-hat for the conditional variance)
Archm6_for_sig <- Archm6_for@forecast$sigmaFor
plot(Archm6_for_sig, type = "l")

# For a better display of the forecast's results, we only take the last 50 observations 
Archm6_for_var_tail <- c(tail(Archm6_est_var,50),rep(NA,10))
Archm6_for_res_tail <- c(tail(Archm6_est_res,50),rep(NA,10))
Archm6_forecast <- c(rep(NA,50),(Archm6_for_sig)^2)

plot(Archm6_for_res_tail, type = "l")
lines(Archm6_for_var_tail, col = "red")
lines(Archm6_forecast, col = "green")


#Saving Dataset for Application 4
save.image("FE_App_4.RData")

