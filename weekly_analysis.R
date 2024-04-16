#necessary libraries
library(TSA)
library(tseries)
library(forecast)
library(MASS)
library(pracma)
#loading the data
weekly_data = read.csv("C:/Users/91938/Desktop/TSA/Final Project/weekly_data.csv", 
                        head = T)$ERCOT
#actual values to compare with forecasted values
actual = window(weekly_data, start = 992, end=1095)
#data for analysis


############################    Data   #########################################
series = window(weekly_data, end = 991)
#converting into time series data
data = ts(series, start=c(2002,1), freq=52)
#raising the plot so that xlab will be visible
par(mar = c(5, 4, 4, 2) + 0.1)
plot(data, xlab='Year(Weeks)', ylab='Load on grid', 
     main = 'ERCOT load curve from 2002 to 2020', col='blue1')
grid(lwd = 0.5, lty = "dotted")
################################################################################


###################  check the normality of data   #############################
qqnorm(data);qqline(data)
hist(data, main='Histogram of Weekly Data', xlab='Load(MW)')
shapiro.test(data)
################################################################################

###########  Analysis(Transformation + Deterministic Modelling)     ############
#boxcox transformation
b = boxcox(lm(data~time(data)))
lambda = b$x[which.max(b$y)]
lambda

#Data after box-cox transformation
bc_data <- (data^lambda - 1) / lambda
plot(bc_data, main = bquote(paste(
  'Box-Cox Transformation with lamba: ', .(lambda), sep="")), 
     ylab='Box-Cox Data')

#checking normality of box-cox data
hist(bc_data, main='Histogram of Box-Cox Data', xlab = 'Box-Cox Data')
shapiro.test(bc_data)

#Applying the Deterministic Trend Modelling
#seasonal means model + Linear Trend
moenth=season(bc_data)
model1=lm(bc_data~month+time(bc_data))
summary(model1)
#residuals
resid_linear = resid(model1)

#converting residuals into timeseries
res = ts(resid_linear, start=c(2002,1), freq=52)
#plotting residuals
plot(res, main='Residuals of Seasonal + Linear Model for Weekly Data', 
     ylab = 'Residuals', xlab= 'Years')
################################################################################

######################   Selection of Candidate Models #########################
#ACF plot of residuals
acf(resid_linear, lag.max = 60, main= 'ACF plot of Residuals') #MA(8)
#PACF of residuals
pacf(resid_linear, lag.max = 60, main= 'PACF plot of Residuals') #AR(1) AR(3)

#Unit root Tests
adf.test(resid_linear)
pp.test(resid_linear)
kpss.test(resid_linear)


#EACF plot of residuals
eacf(resid_linear) # MA(8), ARMA(1,2), ARMA(2,3), ARMA(3,3)
auto.arima(resid_linear) #AR(4)
res=armasubsets(y=rsd4,nar=10,nma=10,y.name='test',ar.method='ols')
par(mfrow=c(1,1))
plot(res)	# default is BIC
#ARMA(3,3)
#AR(4)
#ARMA(4,4)
plot(res,scale='AIC')#ARMA(9,8), ARMA(4,4)
plot(res,scale='AICc')#ARMA(9,8), ARMA(4,4)

#checking AIC, BIC, AICc score
ma8 = Arima(res, order = c(0,0,8), include.mean=F)
ar4 = Arima(res, order = c(4,0,0), include.mean=F)
arma33 = Arima(res, order = c(3,0,3), include.mean=F)
arma44 = Arima(res, order = c(4,0,4), include.mean=F)
#Selected AR(4) as a final model. Now, checking for Overfit
arma41 = Arima(resid_linear, order = c(4,0,1), include.mean=F)
ar5 = Arima(resid_linear, order = c(5,0,0), include.mean=F)
ma8
ar4  #We conclude this as it has best AIC and BIC scores
arma33
arma44
arma41
ar5
#Diagnostic plot of AR(4)
tsdiag(ar4)
################################################################################


################################    Forecast  ##################################
#forecasting for 2 years using AR(4)
ar4_xreg=Arima(bc_data,order=c(4,0,0),include.mean=F,xreg=model.matrix(model1))
pred_years = seq(from=2021,to=2023,length=104)
newpreddata=data.frame(month=as.factor(month[1:104]),tm=pred_years)
predx=predict(ar4_xreg,n.ahead=104,newxreg=
                cbind(1,rbind(diag(1,52)[,-1],diag(1,52)[,-1]),newpreddata[,2]))
pr=predx$pred
#95% Confidence Interval
uci=pr+2*predx$se
lci=pr-2*predx$se
#Converting prediction values and 95% CI into time series data
pr=ts(pr,start=2021,freq=52)
uci=ts(uci,start=2021,freq=52)
lci=ts(lci,start=2021,freq=52)

#calculating RMSE and MAPE
library(Metrics)
rmse(actual, InvBoxCox(pr, lambda = lambda))
mean(abs((actual- InvBoxCox(pr, lambda = lambda))/actual)) * 100

#Actual values which were kept aside to compare with predicted values
actual = ts(actual, start = c(2021,1), freq = 52)


#prediction on the original scale
plot(InvBoxCox(bc_data, lambda = lambda), ylab = "Load on Grid (MW)",
     main = "Forecast from Seasonal + Linear Trend - AR(4)",
     xlim = c(2002, 2023), ylim=c(0, 80000))
lines(InvBoxCox(pr, lambda = lambda),col=2)
lines(InvBoxCox(uci, lambda = lambda),col=3, lty=2)
lines(InvBoxCox(lci, lambda = lambda),col=3, lty=2)


#Zoomed plot of forecast while comparing with actual values
plot(InvBoxCox(bc_data, lambda = lambda), ylab = "Load on Grid (MW)",
     main = "Forecast from Seasonal + Linear Trend - AR(4)",
     xlim = c(2021, 2023), ylim=c(0, 80000),lwd=3) 
lines(InvBoxCox(pr, lambda = lambda),col=2, lwd=3)
lines(InvBoxCox(uci, lambda = lambda),col=3, lty=2, lwd=3)
lines(InvBoxCox(lci, lambda = lambda),col=3, lty=2, lwd=3)
lines(x=seq(from=2021,to=2023,length=104),y=actual,pch=3, lwd=3)
################################################################################


