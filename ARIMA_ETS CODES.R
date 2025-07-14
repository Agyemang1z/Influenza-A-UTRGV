## ARIMA models and Intervention Analysis

# Combine minitab with R studio for ARIMA graphs
library(tseries)
library(ggplot2)
library(reshape2)
library(fUnitRoots)
library(forecast)
library(sarima)
library(lmtest)
library(hwwntest)
library(nortest)
library(nortestARMA)
library(bbmle)
library(seasonal)
library(car)
library(sarima)
library(stats)
library(stats4)
library(TTR)
library(vars)
library(stats)
library(seasonal)
library(forecast)
library(seastests)
library(forecastHybrid)

RTA <- read.csv("C:\\Users\\edmun\\Desktop\\UTRGV Thesis\\ARIMA_ETS\\Data and Codes\\USA.csv")
tail(RTA)
summary(RTA)

#make the data a time series data
Tseries<-ts(RTA$FluAcases, frequency=12, start=c(2009,1))
class(Tseries)
seasonplot(Tseries, col=rainbow(7), year.labels=TRUE,season.labels=TRUE,main="",
           xlab = NULL, ylab = "Number of Influenza A Cases",)

## Plotting the time series of Accident data
autoplot(Tseries, col="green", xlab = "Year", ylab = "Number of Influenza A Cases")


###Test for seasonality
kw(Tseries, freq = 12, diff = T, residuals = F, autoarima = T)

# Splitting the data into two parts
FluATrain <- head(RTA, 168)
FluTest <- tail(RTA, 12)
sum(FluTest$FluAcases)
### Convert Train cases to time series
TseriesTrain<-ts(FluATrain$FluAcases, frequency=12, start=c(2009,1))


## Decomposing the time series
DecompTseries<-decompose(TseriesTrain)
ts.stl<-stl(TseriesTrain,"periodic")  # decompose the TS
ts.sa<-seasadj(ts.stl) #de-seasonalize
acf(DecompTseries$seasonal)
## Plotting the time series of Accident data
autoplot(TseriesTrain, col="green", xlab = "Year", ylab = "Number of Influenza A Cases")

# both acf() and pacf() generates plots by default for Airpassengers
ACFED<- acf(TseriesTrain) # autocorrelation
PACFED<- pacf(TseriesTrain)  # partial autocorrelation

### test for stationarity
kpss.test(TseriesTrain)

# Seasonal Difference
ndiffs(TseriesTrain)  # number for seasonal difference needed

## Difference to make it stationary 
FluAcases_seasdiff <- diff(TseriesTrain, differences=1)  # seasonal differencing
plot(FluAcases_seasdiff, type="l", main="Seasonally Differenced",col="red",xlab="Year")  # still not stationary!
autoplot(FluAcases_seasdiff, col="red")
### test for stationary using Difference series

### KPSS and ADF have opposite hypotheses
kpss.test(FluAcases_seasdiff)
adf.test(FluAcases_seasdiff)

# both acf() and pacf() generates plots by default for Difference series
ACFSEA<- acf(FluAcases_seasdiff) # ACF plot
PACFSEA<- pacf(FluAcases_seasdiff)  # PACF plot


## Arima Function to select the best model
BestArima<-auto.arima(TseriesTrain, stepwise = FALSE, trace = TRUE)

## Checking various ARIMA models

# ARIMA 1
fit1<-Arima(FluAcases_seasdiff, order = c(0,0,1))
summary(fit1)
?arima

# if there is a seasonal component, then the code used is
SARIMA1<-Arima(TseriesTrain, order = c(0,1,3), 
          seasonal = list(order = c(0,0,1), period = 12), 
          include.mean = TRUE,include.drift = FALSE)
# Coefficient test
coeftest(SARIMA1)
## Goodness of fit
summary(SARIMA1) # accuracy test

SARIMA2<-Arima(TseriesTrain, order = c(1,1,1), 
          seasonal = list(order = c(1,0,0), period = 12), 
          include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA2)
## Goodness of fit
summary(SARIMA2) # accuracy test

SARIMA3<-Arima(TseriesTrain, order = c(1,1,0), 
          seasonal = list(order = c(0,0,1), period = 12), 
          include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA3)
## Goodness of fit
summary(SARIMA3) # accuracy test

SARIMA4<-Arima(TseriesTrain, order = c(0,1,1), 
          seasonal = list(order = c(1,0,0), period = 12), 
          include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA4)
## Goodness of fit
summary(SARIMA4) # accuracy test

SARIMA5<-Arima(TseriesTrain, order = c(0,1,2), 
          seasonal = list(order = c(1,0,0), period = 12), 
          include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA5)
## Goodness of fit
summary(SARIMA5) # accuracy test

SARIMA6<-Arima(TseriesTrain, order = c(1,1,2), 
          seasonal = list(order = c(1,0,0), period = 12), 
          include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA6)
## Goodness of fit
summary(SARIMA6) # accuracy test

SARIMA7<-Arima(TseriesTrain, order = c(0,1,3), 
               seasonal = list(order = c(1,0,2), period = 12), 
               include.mean = TRUE,include.drift = FALSE)

# Coefficient test
coeftest(SARIMA7)
## Goodness of fit
summary(SARIMA7) # accuracy test


## Auto generate of ACF and PCF plot
RTA$FluAcases%>%
Arima(order=c(0,1,3), seasonal=c(0,0,1))%>%
residuals() %>% ggtsdisplay()

# Model residual Analysis
checkresiduals(SARIMA1)

###Test for Independence
Box.test(SARIMA1$resid,type="Ljung-Box",lag=12) 

# test for normality of residuals (do for diff lags)
qqPlot(SARIMA1$resid) # Informal test of normality
lillie.test(SARIMA1$resid) # Formal test of normality

######Test for Homoscedasticity (constant variance)
# Obtain the residuals from the SARIMA model
residuals <- SARIMA1$resid

# Perform the Goldfeld-Quandt test on the residuals
gq_test <- gqtest(residuals ~ time(residuals))

# Print the test result
print(gq_test)

# Modified Ljung-Box test for larger lags (model accuracy)
Box.test(SARIMA1$resid,type="Box-Pierce",lag =12) 
Box.test(SARIMA1$resid,type="Box-Pierce",lag=24) 
Box.test(SARIMA1$resid,type="Box-Pierce",lag=36) 
Box.test(SARIMA1$resid,type="Box-Pierce",lag=48)
Box.test(SARIMA1$resid,type="Box-Pierce",lag=60)
Box.test(SARIMA1$resid,type="Box-Pierce",lag=72)

tsdisplay(residuals(SARIMA1))
tsdiag(SARIMA1)
?Box.test

pacf(SARIMA1$resid,col="green")
acf(SARIMA1$resid,col="green")


## Forecasting ARIMA for SARIMA4
Dataforecast<-forecast(SARIMA1,level=0.95, h=12) # for 12 months
plot(Dataforecast,include=24,col="green")
autoplot(Dataforecast,include=12,xlab="Year",ylab="Forecasted Flu Cases")

## Prediction ARIMA for FIT1
?predict
paw<-predict(SARIMA1,n.ahead = 12)

#### Making Plot of Forecast Errors
plot.ts(Dataforecast$residuals) # make a time plot
plotForecastErrors(Dataforecast$residuals) # make a histogram


################### Predictions
# Forecast for 12 months (test period)
Dataforecast <- forecast(SARIMA1, h=12)

# Create full time series: training + test actual
full_actual <- ts(c(as.numeric(TseriesTrain), FluTest$FluAcases), 
                  start = start(TseriesTrain), frequency = 12)

# Forecast values as ts starting after training period
forecast_ts <- ts(Dataforecast$mean, start = c(2023, 1), frequency = 12)

# Define x-axis limits and breaks for years from 2009 to 2024 every 1 year
xlims <- c(2009, 2024)
years_breaks <- seq(2009, 2024, by = 1)

# Plot full actual data (training + test) with no x-axis
plot(full_actual, col = "green", lwd=2, ylab = "Number of Influenza A Cases",
     xlab = "Year",
     xlim = xlims, xaxt = "n")

# Add training data points in green
lines(window(full_actual, end = c(2022,12)), col = "green", lwd = 2)

# Add test data points in red
lines(window(full_actual, start = c(2023,1)), col = "red", lwd = 2)

# Add forecasted values in blue dashed line
lines(forecast_ts, col = "blue", lwd = 2, lty = 2)

# Add custom x-axis with years breaks
axis(1, at = years_breaks, labels = years_breaks)

# Add legend
legend("topleft", legend = c("Training Data", "Test Data", "ARIMA Forecast"),
       col = c("green", "red", "blue"), lwd = 2, lty = c(1,1,2), bty = "n")


######### Prediction errors for SARIMA
ActualFlucases<-c(24955,3947,2376,1378,1341,1156,1377,1168,1703,5535,20126,75147)
ArimaForecast<-c(73487,30077,32007,33424,33753,26981,25501,25431,25930,36679,68314,66776)

######### Calculating Arima error
ArimaError<-ActualFlucases-ArimaForecast
SA<-ArimaError/ActualFlucases
sum(SA)
sum(ArimaError^2)

######### Prediction errors
# Given data
ActualFlucases <- c(24955, 3947, 2376, 1378, 1341, 1156, 1377, 1168, 1703, 5535, 20126, 75147)
ArimaForecast <- c(73487, 30077, 32007, 33424, 33753, 26981, 25501, 25431, 25930, 36679, 68314, 66776)

# Mean Absolute Error (MAE)
mae <- mean(abs(ActualFlucases - ArimaForecast))

# Root Mean Square Error (RMSE)
rmse <- sqrt(mean((ActualFlucases - ArimaForecast)^2))

# Mean Absolute Percentage Error (MAPE)
mape <- mean(abs((ActualFlucases - ArimaForecast) / ActualFlucases)) * 100

# Symmetric Mean Absolute Percentage Error (SMAPE)
smape <- mean(2 * abs(ActualFlucases - ArimaForecast) / (abs(ActualFlucases) + abs(ArimaForecast))) * 100

# Median Absolute Percentage Error (MdAPE)
mdape <- median(abs((ActualFlucases - ArimaForecast) / ActualFlucases)) * 100

# Geometric Mean Relative Absolute Error (GMRAE)
# First, calculate the benchmark forecast using naive method
naive_forecast <- c(NA, ActualFlucases[-length(ActualFlucases)])  # lagged actual values
naive_forecast[1] <- ActualFlucases[1]  # Handle the first value for the benchmark
gmrae <- exp(mean(log(abs((ActualFlucases[-1] - ArimaForecast[-1]) / (ActualFlucases[-1] - naive_forecast[-1])))))

# Theil U1 statistic
numerator <- sqrt(mean((ArimaForecast - ActualFlucases)^2))
denominator <- sqrt(mean(ArimaForecast^2)) + sqrt(mean(ActualFlucases^2))
theil_u1 <- numerator / denominator

# Print results
cat("MAE:", mae, "\n")
cat("RMSE:", rmse, "\n")
cat("MAPE:", mape, "%\n")
cat("SMAPE:", smape, "%\n")
cat("MdAPE:", mdape, "%\n")
cat("GMRAE:", gmrae, "\n")
cat("Theil U1:", theil_u1, "\n")

###########EXPONENTIAL SMOOTHING WITH HOLT WINTER'S METHOD

### Whole Data (up to 2023)
RTAtimeseriesforecasts <- HoltWinters(Tseries)

plot(RTAtimeseriesforecasts)

# Forecast the time series
RTAtimeseriesforecasts_forecast <- forecast(RTAtimeseriesforecasts, h=12)

# Plot the forecast
plot(RTAtimeseriesforecasts_forecast)

Box.test(RTAtimeseriesforecasts_forecast$residuals, lag=20, type="Ljung-Box")

plot.ts(RTAtimeseriesforecasts_forecast$residuals) # make a time plot
plotForecastErrors(RTAtimeseriesforecasts_forecast$residuals)

# test for normality of residuals (do for diff lags)
qqPlot(RTAtimeseriesforecasts_forecast$residuals) # Informal test of normality
lillie.test(RTAtimeseriesforecasts_forecast$residuals) # Formal test of normality


####### Train Data (up to 2022)
RTA1timeseriesforecasts <- HoltWinters(TseriesTrain)

plot(RTA1timeseriesforecasts)

# Forecast the time series
RTA1timeseriesforecasts_forecast <- forecast(RTA1timeseriesforecasts, h=12)

# Plot the forecast
plot(RTA1timeseriesforecasts_forecast)



#########Prediction Plot
# Fit Holt-Winters model on training data
RTA1timeseriesforecasts <- HoltWinters(TseriesTrain)

# Forecast for 12 months (test period)
RTA1_forecast <- forecast(RTA1timeseriesforecasts, h=12)

# Create full time series: training + test actual
full_actual <- ts(c(as.numeric(TseriesTrain), FluTest$FluAcases), 
                  start = start(TseriesTrain), frequency = 12)

# Forecast values as ts starting after training period
forecast_ts <- ts(RTA1_forecast$mean, start = c(2023, 1), frequency = 12)

# Define x-axis limits and breaks for years from 2009 to 2025 every 2 years
xlims <- c(2009, 2024)
years_breaks <- seq(2009, 2024, by = 1)

# Plot full actual data (training + test) with no x-axis
plot(full_actual, col = "green", lwd=2, ylab = "Number of Influenza A Cases",
     xlab = "Year",
     xlim = xlims, xaxt = "n")

# Add training data points in green
lines(window(full_actual, end = c(2022,12)), col = "green", lwd = 2)

# Add test data points in red
lines(window(full_actual, start = c(2023,1)), col = "red", lwd = 2)

# Add forecasted values in blue dashed line
lines(forecast_ts, col = "blue", lwd = 2, lty = 2)

# Add custom x-axis with year labels
axis(1, at = years_breaks, labels = years_breaks)

# Add legend
legend("topleft", legend = c("Training Data", "Test Data", "ETS Forecast"),
       col = c("green", "red", "blue"), lwd = 2, lty = c(1,1,2), bty = "n")

######### Prediction errors
# Given data
ActualFlucases <- c(24955, 3947, 2376, 1378, 1341, 1156, 1377, 1168, 1703, 5535, 20126, 75147)
ETSForecast <- c(141720, 155289,157984 , 144208,129951 , 112367, 113712, 118415,123645 , 143343,187241 , 170518)

# Mean Absolute Error (MAE)
mae <- mean(abs(ActualFlucases - ETSForecast))

# Root Mean Square Error (RMSE)
rmse <- sqrt(mean((ActualFlucases - ETSForecast)^2))

# Mean Absolute Percentage Error (MAPE)
mape <- mean(abs((ActualFlucases - ETSForecast) / ActualFlucases)) * 100

# Symmetric Mean Absolute Percentage Error (SMAPE)
smape <- mean(2 * abs(ActualFlucases - ETSForecast) / (abs(ActualFlucases) + abs(ETSForecast))) * 100

# Median Absolute Percentage Error (MdAPE)
mdape <- median(abs((ActualFlucases - ETSForecast) / ActualFlucases)) * 100

# Geometric Mean Relative Absolute Error (GMRAE)
# First, calculate the benchmark forecast using naive method
naive_forecast <- c(NA, ActualFlucases[-length(ActualFlucases)])  # lagged actual values
naive_forecast[1] <- ActualFlucases[1]  # Handle the first value for the benchmark
gmrae <- exp(mean(log(abs((ActualFlucases[-1] - ETSForecast[-1]) / (ActualFlucases[-1] - naive_forecast[-1])))))

# Theil U1 statistic
numerator <- sqrt(mean((ETSForecast - ActualFlucases)^2))
denominator <- sqrt(mean(ETSForecast^2)) + sqrt(mean(ActualFlucases^2))
theil_u1 <- numerator / denominator

# Print results
cat("MAE:", mae, "\n")
cat("RMSE:", rmse, "\n")
cat("MAPE:", mape, "%\n")
cat("SMAPE:", smape, "%\n")
cat("MdAPE:", mdape, "%\n")
cat("GMRAE:", gmrae, "\n")
cat("Theil U1:", theil_u1, "\n")
