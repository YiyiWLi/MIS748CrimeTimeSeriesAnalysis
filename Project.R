rm(list=ls())
library(tseries) 
library(TSA)
library(lmtest)
library(tsoutliers)
source("rolling.forecast.R")

setwd("/Users/yiyi/MIS_748/FinalProject")
cr = read.csv("data/crimebyweekly.csv")

data = ts(cr[1:269,2], start=c(2017,1), frequency = 52)
log_data = log(data)
plot(data, main="Total Number of Crime from 2017 to 2022", ylab="Number")

train_crime = ts(cr[1:167,2], start=c(2017,1), frequency = 52)
test_crime = ts(cr[168:269,2], start=c(2020,11), frequency = 52)

# log transformation
log_cr = log(train_crime)
acf(log_cr, lag.max = 150)
pacf(log_cr, lag.max = 150)

# seasonal difference
dlog_cr = diff(log_cr,52)
acf(dlog_cr, lag.max = 150)
pacf(dlog_cr, lag.max = 150)

# regular difference
ddlog_cr = diff(dlog_cr)
par(mfrow = c(1,2))
acf(ddlog_cr, main="ACF for Transformed Data", lag.max = 150)
pacf(ddlog_cr, main="PACF for Transformed Data", lag.max = 150)

# check stationary
adf.test(ddlog_cr)

# build model without p,q,d orders
model_try = arima(ddlog_cr, order = c(0, 0, 0), seasonal = list(order=c(0, 0, 0), period=52))
model_try
# AIC = -253.78
coeftest(model_try)
# intercept is not significant

# build model with d order
model_try1 = arima(log_cr, order = c(0, 1, 0), seasonal = list(order=c(0, 1, 0), period=52))
model_try1
# AIC = -255.75
acf(model_try1$residuals, main="ACF for Residual of ARIMA(0,1,0)*(0,1,0)")
# MA(1)
pacf(model_try1$residuals, main="PACF for Residual of ARIMA(0,1,0)*(0,1,0)")
# AR(3)
eacf(model_try1$residuals)
# ARMA(0,1)

# MA(1)
ma = arima(log_cr, order = c(0, 1, 1), seasonal = list(order=c(0, 1, 0), period=52))
ma
# AIC = -315.92
coeftest(ma)
Box.test(ma$residuals, type = 'Ljung',lag = 60)
# p-value = 0.9793

# AR(3)
ar = arima(log_cr, order = c(3, 1, 0), seasonal = list(order=c(0, 1, 0), period=52))
# AIC = -299.29
coeftest(ar)
Box.test(ar$residuals, type = 'Ljung',lag = 60)
# p-value = 0.8025
abs(polyroot(c(1, -ar$coef[1:3])))
# 1.588750 1.817045 1.588750 out of unit 

# Comparision fitted values and actual values
par(mfrow = c(1,1))
# MA(1)
# fit = ts(exp(fitted.values(ma)), start = c(2017,1), frequency = 52)
# Error in exp(fitted.values(ar)) : 
#   non-numeric argument to mathematical function
ma_fit = ts(exp(log_cr-ma$residuals), start = c(2017,1), frequency = 52)
plot(train_crime, ylim=c(min(train_crime)-200, max(ma_fit)+200), xlim = c(2017,2021),xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Fitted Values and Actual Values with ARIMA(0,1,1)*(0,1,0)',lwd=2)
points(ma_fit, col = 'red', pch = 2,cex=0.5)
lines(ma_fit, col = 'red')
legend.text=c("Actual values", "Fitted values")
legend("bottomright", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# AR(3)
# fit = ts(exp(fitted.values(ar)), start = c(2017,1), frequency = 52)
ar_fit = ts(exp(log_cr-ar$residuals), start = c(2017,1), frequency = 52)
plot(train_crime, ylim=c(min(train_crime)-200, max(ar_fit)+200), xlim=c(2017,2021),xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Fitted Values and Actual Values with ARIMA(3,1,0)*(0,1,0)',lwd=2)
points(ar_fit, col = 'red', pch = 2,cex=0.5)
lines(ar_fit, col = 'red')
legend.text=c("Actual values", "Fitted values")
legend("bottomright", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# rolling forecast
error1 = rolling.forecast(log_cr, 7, 160, c(0,1,1), seasonal = list(order=c(0,1,0),period=52),include.mean=F)
# 0.001437879 0.001504153 0.000618940 0.002150202 0.004890202 0.006775619 0.003485259
error2 = rolling.forecast(log_cr, 7, 160, c(3,1,0), seasonal = list(order=c(0,1,0),period=52),include.mean=F)
# 0.0011721909 0.0003950582 0.0008873423 0.0008372706 0.0023424955 0.0035833035 0.0015633521

# Prediction for 2020 to 2022
# MA(1)
ma_pred = predict(ma, 102)
ma_pred_ts = ts(exp(ma_pred$pred), start = c(2020,11), frequency=52)
plot(data, xlim=c(2018,2023), ylim=c(500,2500), xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Prediction and Actual Values from 2020 to 2022 with ARIMA(0,1,1)*(0,1,0)',lwd=2)
points(ma_pred_ts, col = 'red', pch = 2,cex=0.5)
lines(ma_pred_ts, col = 'red',cex=1.5)
legend.text=c("Actual values", "Predicted values")
legend("bottomleft", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# AR(3)
ar_pred = predict(ar, 102)
ar_pred_ts = ts(exp(ar_pred$pred), start = c(2020,11), frequency=52)
plot(data, xlim=c(2018,2023), ylim=c(500,2500), xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Prediction and Actual Values from 2020 to 2022 with ARIMA(3,1,0)*(0,1,0)',lwd=2)
points(ar_pred_ts, col = 'red', pch = 2,cex=0.5)
lines(ar_pred_ts, col = 'red',cex=1.5)
legend.text=c("Actual values", "Predicted values")
legend("bottomleft", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# check outliers
# MA(1)
model_ma = arima(log_data, order = c(0, 1, 1), seasonal = list(order=c(0, 1, 0), period=52))
locate.outliers(model_ma$residuals, pars=coefs2poly(model_ma))
# AR(3)
model_ar = arima(log_data, order = c(3, 1, 0), seasonal = list(order=c(0, 1, 0), period=52))
locate.outliers(model_ar$residuals, pars=coefs2poly(model_ar))

# add intervention variables
n = length(data)
n1 = 167
pulse = ts(c(rep(0, n1), 1, rep(0, n-n1-1)), start=2017, frequency=52) 
step=ts(c(rep(0, n1), rep(1, n-n1)), start=2017, frequency=52) 
# MA(1)
order=c(0,1,1)
seasonal=list(order=c(0,1,0), period=52)
# no intervention
ma00= arimax(log_data, order=order, seasonal=seasonal, method="ML")
ma00
# AIC = -473.73
# with step function
ma01= arimax(log_data, order=order, seasonal=seasonal, xtransf=data.frame(step), transfer=list(c(0,0)), method="ML")
ma01
# AIC = -548.13
# with pulse function
ma10= arimax(log_data, order=order, seasonal=seasonal, xtransf=data.frame(pulse), transfer=list(c(1,0)), method="ML")
ma10
# AIC = -587.77
# with step and pulse function
ma11= arimax(log_data, order=order, seasonal=seasonal, xtransf=data.frame(pulse, step), transfer=list(c(1,0), c(0,0)), method="ML")
ma11
# AIC = -586.6

# AR(3)
order1=c(3,1,0)
seasonal1=list(order=c(0,1,0), period=52)
# no intervention
ar00= arimax(log_data, order=order1, seasonal=seasonal1, method="ML")
ar00
# AIC = -471.71
# with step function
ar01= arimax(log_data, order=order1, seasonal=seasonal1, xtransf=data.frame(step), transfer=list(c(0,0)), method="ML")
ar01
# AIC = -542.64
# with pulse function
ar10= arimax(log_data, order=order1, seasonal=seasonal1, xtransf=data.frame(pulse), transfer=list(c(1,0)), method="ML")
ar10
# AIC = -556.91
# with step and pulse function
ar11= arimax(log_data, order=order1, seasonal=seasonal1, xtransf=data.frame(pulse, step), transfer=list(c(1,0), c(0,0)), method="ML")
ar11
# AIC = -554.91

# intervention variables
# MA(1)
ma_tc = filter(pulse, filter=ma10$coef[2],method='recursive',side=1)*ma10$coef[3]
plot(ma_tc, type='l', xlim=c(2017,2023), ylab='impact',main="Intervention Change with ARIMA(0,1,1)*(0,1,0)")

ma_pred1 = arima(log_data, order=order, seasonal=seasonal, xreg=data.frame(ma_tc))
ma_pred1 
# AIC = -589.77
ma_fit1 = ts(exp(log_data-ma10$residuals), start = c(2017,1), frequency=52)
plot(data, xlim=c(2017,2023), ylim=c(0,2500), xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Fitted Values and Actual Values from 2020 to 2022 with ARIMA(0,1,1)*(0,1,0)',lwd=2)
points(ma_fit1, col = 'red', pch = 2,cex=0.5)
lines(ma_fit1, col = 'red',cex=1.5)
legend.text=c("Actual values", "Predicted values")
legend("bottomleft", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# AR(3)
ar_tc = filter(pulse, filter=ar10$coef[4],method='recursive',side=1)*ar10$coef[5]
plot(ar_tc, type='l', xlim=c(2017,2023), ylab='impact',main="Intervention Change with ARIMA(3,1,0)*(0,1,0)")

ar_pred1 = arima(log_data, order=order, seasonal=seasonal, xreg=data.frame(ar_tc))
ar_pred1 
# AIC = -589.77
ar_fit1 = ts(exp(log_data-ar10$residuals), start = c(2017,1), frequency=52)
plot(data, xlim=c(2017,2023), ylim=c(0,2500), xlab = 'Year', ylab = 'Crime Number', main = 'Comparision between Fitted Values and Actual Values from 2020 to 2022 with ARIMA(3,1,0)*(0,1,0)',lwd=2)
points(ar_fit1, col = 'red', pch = 2,cex=0.5)
lines(ar_fit1, col = 'red',cex=1.5)
legend.text=c("Actual values", "Predicted values")
legend("bottomleft", legend.text, lty = rep(1,2), col = 1:2, pch = 1:2, cex=0.6)

# predition for future
h=44
pulse_new = ts(c(rep(0, n1), 1, rep(0, n+h-n1-1)), start=2017, frequency=52) 
step_new = ts(c(rep(0, n1), rep(1, n+h-n1)), start=2017, frequency=52) 
# MA(1)
ma_tc_new = filter(pulse_new, filter=ma10$coef[2],method='recursive',side=1)*ma10$coef[3]
ma_newdata=data.frame(tc = tc_new[(n+1):(n+h)])
ma_pred_full=predict(ma_pred1, h, newxreg=ma_newdata)

# AR(3)
ar_tc_new = filter(pulse_new, filter=ar10$coef[4],method='recursive',side=1)*ma10$coef[5]
ar_newdata=data.frame(tc = tc_new[(n+1):(n+h)])
ar_pred_full=predict(ar_pred1, h, newxreg=ar_newdata)

time = seq(2017, 2022+8/52, by=1/52)
time_new = seq(2022+9/52, 2023, by=1/52)
par(mfrow = c(1,1))
# plot for MA(1)
plot(time, data, xlab='Time', ylab='', xlim=c(2017, 2023), ylim=c(500,2500), type='o', main='Prediction of Crime in 2022 with ARIMA(0,1,1)*(0,1,0)',cex=0.5)
points(time, exp(log_data-ma10$residuals), type='o', pch=2, lty=1, col='red',cex=0.5)
points(time_new, exp(ma_pred_full$pred), type='o', pch=3, lty=1, col='blue',cex=0.5)
legend.txt=c("Actual", "Fitted", "Predicted")
legend("bottomright", legend.txt, col=c("black", "red", "blue"), pch=c(1,2,3), lty=c(1,2,1), cex=0.5)

# plot for AR(3)
plot(time, data, xlab='Time', ylab='', xlim=c(2017, 2023), ylim=c(500,2500), type='o', main='Prediction of Crime in 2022 with ARIMA(3,1,0)*(0,1,0)',cex=0.5)
points(time, exp(log_data-ar10$residuals), type='o', pch=2, lty=1, col='red',cex=0.5)
points(time_new, exp(ar_pred_full$pred), type='o', pch=3, lty=1, col='blue',cex=0.5)
legend.txt=c("Actual", "Fitted", "Predicted")
legend("bottomright", legend.txt, col=c("black", "red", "blue"), pch=c(1,2,3), lty=c(1,2,1), cex=0.5)

