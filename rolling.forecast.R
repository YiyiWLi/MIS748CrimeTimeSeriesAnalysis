## require library(tseries)
#
# perform rolling forecast and obtain sum of out-of-sample forecasting error
# input
## x: time series
## max.step: maximum forecasting horizon
## start: starting point for out-of-sample forecast
## order: ARIMA order (as called in function arima)
## seasonal: seasonal ARIMA order: (as called in function arima
#
rolling.forecast=function(x,max.step=1,start=50,order=c(0,0,0),
       seasonal=list(order=c(0,0,0),period=NA),include.mean=T,fixed=NULL,xreg=NULL){
nn=length(x)
kk=max.step
ss=(1:kk)*0
for (i in start:(nn-kk-1)) {
xx=x[1:i]
out1=arima(xx,order,seasonal,include.mean=include.mean,fixed=fixed,xreg=xreg)
ff=predict(out1,kk)
for (k in 1:kk) ss[k]=ss[k]+(ff$pred[k]-x[i+k])**2
}
return(ss)
}

