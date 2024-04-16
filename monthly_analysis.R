library(TSA)
library(tseries)
library(forecast)
library(MASS)
library(pracma)

monthly_data = read.csv("C:/Users/udayr/Documents/MA 5781/Final Project/ERCOT data/monthly_data.csv",head=T)
data = ts(monthly_data$ERCOT[1:216], start=c(2002,1), freq=12) # Monthly frequency

# Fit the model on data from 2002 to 2020 - 1:228
# Use remaining data to compare the forecast - 229:252

plot(data, xlab="Years", ylab="Load on Grid (MW)", main="ERCOT Load Curve data from 2002 to 2020")
time(data)
par(mfrow=c(1,2))
hist(data)
qqnorm(data,main="QQ plot - 2002 to 2020")
qqline(data)
par(mfrow=c(1,1))

shapiro.test(data)
# Distribution seems right skewed

# Attempting transforms - Box-Cox
# Parameter estimation for using Box-Cox transformation
b = boxcox(lm(data~time(data)))
# Exact lambda
lambda = b$x[which.max(b$y)]
lambda
# Estimated parameter is somewhere -1.23, which is near -1 corresponding to 1/x
tx_data = (data ^ lambda - 1) / lambda

par(mfrow=c(3,1))
plot(tx_data, xlab="Years", ylab="log(MW)", main="Box-Cox tranformation with lamda = -1.23")
hist(tx_data)
qqnorm(tx_data,main="QQ plot of Box-Cox transform ")
qqline(tx_data)
par(mfrow=c(1,1))

shapiro.test(tx_data)
# Test says non-normal but at least the skew of the data has been corrected.

# Seasonal means plus linear trend on transformed series --------------------------------------------------------------------------------
tm = time(tx_data)
trend = lm(tx_data~season(tx_data) + tm)
summary(trend)
rsd = ts(resid(trend))

par(mfrow=c(3,1))
plot(rsd, ylab="Residuals",xlab="Years", main = sprintf("Residuals of Seasonal + Linear model"))
acf(rsd, main="ACF of residuals", lag.max = 100)
pacf(rsd, main="PACF of residuals", lag.max = 100)
par(mfrow=c(1,1))
# Maybe MA(2) or AR(1) models based on acf and pacf

par(mfrow=c(1,2))
hist(rsd, main="Histogram - Seasonal + Linear")
qqnorm(rsd,main="QQ plot - Seasonal + Linear")
qqline(rsd)
par(mfrow=c(1,1))

adf.test(rsd)
pp.test(rsd)
kpss.test(rsd)
# No unit roots found

auto.arima(rsd)
# Suggests ARIMA(0,0,0) with zero mean.

eacf(rsd, ar.max = 13)
# Vertex found at --- ARMA(4,0), ARMA(4,2)

res=armasubsets(y=rsd,nar=10,nma=10,y.name='test',ar.method='ols')
par(mfrow=c(1,1))
plot(res)	# default is BIC
# BIC rankings ---
# ARMA(1,0)
# ARMA(0,10)

plot(res,scale='AIC')
plot(res,scale='AICc')
# Both AIC and AICc rankings ---
# ARMA(4,10)
# ARMA(1,0) 
par(mfrow=c(1,1))

# Using AR(1)
month = season(tx_data)
arma10_xreg=Arima(tx_data, order=c(1,0,0),include.mean=F,xreg=model.matrix(trend))
pred_years = seq(from=2020,to=2021.917,length=24)
newpreddata=data.frame(month=as.factor(month[1:12]),tm=pred_years)
predx=predict(arma10_xreg,n.ahead=24,newxreg=cbind(1,rbind(diag(1,12)[,-1],diag(1,12)[,-1]),newpreddata[,2]))
pr=predx$pred
uci=pr+2*predx$se
lci=pr-2*predx$se


pr=ts(pr,start=2020,freq=12)
uci=ts(uci,start=2020,freq=12)
lci=ts(lci,start=2020,freq=12)

plot(tx_data,xlim=c(2002,2022), xlab="Years", main="Transformed Predictions - ARMA(1,0)")
lines(pr,col=2)
lines(uci,col=3)
lines(lci,col=3)

plot(InvBoxCox(tx_data,lambda = lambda), ylab="Load (in MW)", xlab="Years", main="Original Scale retreived - ARMA(1,0)",xlim=c(2002,2022))
lines(InvBoxCox(pr,lambda = lambda),col=2)
lines(InvBoxCox(uci,lambda = lambda),col=3)
lines(InvBoxCox(lci,lambda = lambda),col=3)

# Model as white noise
month = season(tx_data)
arma00_xreg=Arima(tx_data, order=c(0,0,0),include.mean=F,xreg=model.matrix(trend))
pred_years = seq(from=2020,to=2021.917,length=24)
newpreddata=data.frame(month=as.factor(month[1:12]),tm=pred_years)
predx=predict(arma00_xreg,n.ahead=24,newxreg=cbind(1,rbind(diag(1,12)[,-1],diag(1,12)[,-1]),newpreddata[,2]))
pr=predx$pred
uci=pr+2*predx$se
lci=pr-2*predx$se


pr=ts(pr,start=2020,freq=12)
uci=ts(uci,start=2020,freq=12)
lci=ts(lci,start=2020,freq=12)

plot(tx_data,xlim=c(2002,2022), xlab="Years", main="Transformed Predictions - White Noise")
lines(pr,col=2)
lines(uci,col=3)
lines(lci,col=3)

plot(InvBoxCox(tx_data,lambda = lambda), ylab="Load (in MW)", xlab="Years", main="Original Scale retreived - White Noise",xlim=c(2002,2022))
lines(InvBoxCox(pr,lambda = lambda),col=2)
lines(InvBoxCox(uci,lambda = lambda),col=3)
lines(InvBoxCox(lci,lambda = lambda),col=3)