#資料挖礦與大數據分析============================================================

install.packages("forecast")
library(forecast)
install.packages("TSA")
library(TSA)
data(co2, package = "datasets")
tsdisplay(co2)
acf(co2)
pacf(co2)

train <- ts(co2[seq(1, length(co2)-12)], frequency = 12, start = c(1959, 1))
test <- ts(co2[seq(length(co2)-11, length(co2))], frequency = 12, start = c(1997, 1)) #留著和預測做比對

tsdisplay(diff(train), main = "First difference of co2") #季節性時間位差的相關性仍十分明顯
acf(diff(train), main = "First difference of co2")
pacf(diff(train), main = "First difference of co2")

tsdisplay(diff(diff(train), lag = 12), main = "First and seasonal difference of co2") #季節性差分
acf(diff(diff(train), lag = 12), main = "First and seasonal difference of co2")
pacf(diff(diff(train), lag = 12), main = "First and seasonal difference of co2")

m1 <- Arima(train, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 1), period = 12))
names(m1)
m1$coef #ma2參數不顯著（估計值的絕對值小於等於兩倍標準誤）
m1$var.coef #參數估計間的共變異都非常小，代表參數估計之間不會相互混淆

tsdiag(m1, gof = 36)
qqnorm(residuals(m1))
qqline(residuals(m1))
legend("topleft", legend = paste("p-value =", round(shapiro.test(residuals(m1))$p.val, 4)))


fix.par <- rep(NA, 10)
fix.par[c(2, 4:8)] = 0
m1.1 <- Arima(train, order = c(0, 1, 9), seasonal = list(order = c(0, 1, 1), period = 12), fixed = fix.par)

m1.1$coef
m1.1$var.coef
tsdiag(m1.1, gof = 36)
qqnorm(residuals(m1.1))
qqline(residuals(m1.1))
legend("topleft", legend = paste("p-value =", round(shapiro.test(residuals(m1.1))$p.val, 4)))

m1 #觀察AIC, BIC, sigma^2
m1.1 #觀察AIC, BIC, sigma^2

pred1 <- predict(m1, n.ahead = 12, se.fit = TRUE)$pred
MAE <- round(sum(abs(pred1 - test) / 12), 2)
MSE <- round(sum((pred1 - test)^2) / 12, 2)
MAPE <- round(sum(abs(pred1 - test) / pred1) / 12 * 100, 2)


pred2 <- predict(m1.1, n.ahead = 12, se.fit = TRUE)$pred
MAE <- round(sum(abs(pred2 - test) / 12), 2)
MSE <- round(sum((pred2 - test)^2) / 12, 2)
MAPE <- round(sum(abs(pred2 - test) / pred2) / 12 * 100, 2)


#R軟體資料分析基礎與應用============================================================

##自迴歸移動平均模型（Autoregressive Moving Average）============================================================

require(WDI)
gdp <- WDI(country = c("US", "CA", "GB", "DE", "CN", "JP", "SG", "IL"), indicator = c("NY.GDP.PCAP.CD", "NY.GDP.MKTP.CD"), start = 1960, end = 2011)
names(gdp) <- c("iso2c", "Country", "Year", "PerCapGDP", "GDP")
head(gdp)

require(ggplot2)
require(scales)
ggplot(gdp, aes(Year, PerCapGDP, color = Country, linetye = Country)) + geom_line() + scale_y_continuous(label = dollar)

require(useful)
ggplot(gdp, aes(Year, GDP, color = Country, linetype = Country)) + geom_line() + scale_y_continuous(label = multiple_format(extra = dollar, multiple = "M"))


us <- gdp$PerCapGDP[gdp$Country == "United States"]
us <- ts(us, start = min(gdp$Year), end = max(gdp$Year)) #轉換為時間序列
us
plot(us, ylab = "Per Capita GDP", xlab = "Year")
acf(us) #自相關函數（autocorrelation function; ACF）##ACF和PACF用於看平穩性
pacf(us) #偏自相關函數（partial autocorrelation function; PACF）

x <- c(1, 4, 8, 2, 6, 6, 5, 3) #需要對時間序列做轉換才能建模，對該序列做差分或做轉換
diff(x, differences = 1) #一階差分
diff(x, differences = 2) #二階差分
diff(x, lag = 1) #元素和其之前第一個元素的差
diff(x, lag = 2) #元素和其之前第二個元素的差

require(forecast)
ndiffs(x = us) #找出做幾階層差分才是最合適的
plot(diff(us, 2))

usBest <- auto.arima(x = us,2) #auto.arima()挑出AR(2)和MA(1)組成最佳模型，挑選基準為最小AICC
usBest #好的模型殘差應呈現白噪音（white noise）的特性

acf(usBest$residuals)
pacf(usBest$residuals)
coef(usBest)

predict(usBest, n.ahead = 5, se.fit = TRUE)
theForecast <- forecast(object = usBest, h = 5)
plot(theForecast) #陰影為信賴區間


##VAR向量自我迴歸============================================================

require(reshape2)
gdpCast <- dcast(Year ~ Country, data = gdp[, c("Country", "Year", "PerCapGDP")], value.var = "PerCapGDP")
head(gdpCast)
gdpTS <- ts(data = gdpCast[, -1], start = min(gdpCast$Year), end = max(gdpCast$Year))
gdpTS
plot(gdpTS, plot.type = "single", col = 1:8)
legend("topleft", legend = colnames(gdpTS), ncol = 2, lty = 1, col = 1:8, cex = .9)

gdpTS <- gdpTS[, which(colnames(gdpTS) != "Germany")]
numDiffs <- ndiffs(gdpTS)
numDiffs

gdpDiffed <- diff(gdpTS, differences = numDiffs)
plot(gdpDiffed, plot.type = "single", col = 1:7)
legend("bottomleft", legend = colnames(gdpDiffed), ncol = 2, lty = 1, col = 1:7, cex = .9)


install.packages("vars")
require(vars)
gdpVar <- VAR(gdpDiffed, lag.max = 12) #建立模型
gdpVar$p #挑選的位階
names(gdpVar$varresult)
class(gdpVar$varresult$Canada) #每個模型都是lm物件
class(gdpVar$varresult$Japan)
head(coef(gdpVar$varresult$Canada)) #每個模型都有各自的係數
head(coef(gdpVar$varresult$Japan))

require(coefplot)
coefplot(gdpVar$varresult$Canada) #VAR模型係數圖
coefplot(gdpVar$varresult$Japan)

predict(gdpVar, n.ahead = 5)


##GARCH============================================================

install.packages("quantmod")
require(quantmod)
att <- getSymbols("T", auto.assign = FALSE) #AT&T股價
head(att)
plot(att)
chartSeries(att)
addBBands()
addMACD(32, 50, 12)


attClose <- att$T.Close #我們感興趣的是收盤價
class(attClose)
head(attClose)


install.packages("rugarch")
require(rugarch) #以ugarchspec()設定模型規格
attSpec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), #指定以GARCH(1, 1)對波動度建模
                      mean.model = list(armaOrder = c(1, 1)), #以ARMA(1, 1)對平均值建模
                      distribution.model = "std") #指定創新分佈為t分佈

attGarch <- ugarchfit(spec = attSpec, data = attClose) #以ugarchfit()建立模型
attGarch

plot(attGarch@fit$residuals, type = "l") #殘差 #attGarch是一個S4物件，裡面的資料要用@來套用
plot(attGarch, which = 10) #殘差ACF


attSpec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(1, 1)), 
                       distribution.model = "std") #以不同規格設定平均數，比較AIC

attSpec2 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                       mean.model = list(armaOrder = c(0, 0)),
                       distribution.model = "std")

attSpec3 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(0, 2)),
                       distribution.model = "std")

attSpec4 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = c(1, 2)),
                       distribution.model = "std")

attGarch1 <- ugarchfit(spec = attSpec1, data = attClose)
attGarch2 <- ugarchfit(spec = attSpec2, data = attClose)
attGarch3 <- ugarchfit(spec = attSpec3, data = attClose)
attGarch4 <- ugarchfit(spec = attSpec4, data = attClose)

infocriteria(attGarch1) #Akaike(AIC); Bayes(BIC) 
infocriteria(attGarch2)
infocriteria(attGarch3)
infocriteria(attGarch4) #attGarch1 & 4為最佳模型

attPred <- ugarchboot(attGarch, n.ahead = 50, method = c("Partial", "Full")[1])
plot(attPred, which = 2)


attLog <- diff(log(attClose))[-1]
attLogSpec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(1, 1)),
                         distribution.model = "std")
attLogGarch <- ugarchfit(spec = attLogSpec, data = attLog)
infocriteria(attLogGarch) #GARCH模型目的不是要對訊號（signal）建立更好的模型，而是要更好地捕捉波動度（volatility）的行為


#Practical Machine Learning==================================================

#Forecasting==================================================

library(quantmod)
from.dat <- as.Date("01/01/08", format = "%m/%d/%y")
to.dat <- as.Date("12/31/13", format = "%m/%d/%y")
getSymbols("GOOG", src = "google", from = from.dat, to = to.dat)

head(GOOG)

mGoog <- to.monthly(GOOG) #ERROR，待查
googOpen <- Op(mGoog)
ts1 <- ts(googOpen, frequency = 12)
plt(ts1, xlab = "Years+1", ylab = "GOOG")

ts1Train <- window(ts1, start = 1, end = 5)
ts1Test <- window(ts1, start = 5, end = (7 - 0.01))
ts1Train

plot(ts1Train)
lines(ma(ts1Train, order = 3), col = "red")

ets1 <- ets(ts1Train, model = "MMM")
fcast <- forecast(ets1)
plot(fcast)
lines(ts1Test, col = "red")
accuracy(fcast, ts1Test)