#R軟體資料分析基礎與應用==================================================

#17廣義線性模型==================================================

##羅吉斯迴歸（Logistic Regression）==================================================

acs <- read.table("http://jaredlander.com/data/acs_ny.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
acs$Income <- with(acs, FamilyIncome >= 150000)
require(ggplot2)
require(useful)
ggplot(acs, aes(x = FamilyIncome)) + geom_density(fill = "grey", color = "grey") + geom_vline(xintercept = 150000) + scale_x_continuous(label = multiple.dollar, limits = c(0, 1000000))

head(acs)

income1 <- glm(Income ~ HouseCosts + NumWorkers + OwnRent + NumBedrooms + FamilyType, 
               data = acs, family = binomial(link = "logit"))

summary(income1)

invlogit <- function(x){
        1/(1 + exp(-x))
}
invlogit(income1$coefficients) #解讀羅吉斯迴歸中的係數要先用羅吉反函數對係數作轉換！


##泊松迴歸模型==================================================

#適合運用在計數資料上

ggplot(acs, aes(x = NumChildren)) + geom_histogram(binwidth = 1) #家庭小孩數量直方圖

children1 <- glm(NumChildren ~ FamilyIncome + FamilyType + OwnRent, data = acs, 
                 family = poisson(link = "log"))
summary(children1)

install.packages("coefplot")
library(coefplot)
coefplot(children1)

#泊松迴歸可能存在過度離散（overdispersion）的問題

z <- (acs$NumChildren - children1$fitted.values) / sqrt(children1$fitted.values)
sum(z ^ 2) / children1$df.residual #計算過度離散率（OD），>= 2表示過度離散，或是< 2但p值為1

pchisq(sum(z ^ 2), children1$df.residual) #過度離散p值

#過度離散問題顯著，用準泊松分佈（quasi poisson）或負二項分佈（negative binomial）重新建模

children2 <- glm(NumChildren ~ FamilyIncome + FamilyType + OwnRent, data = acs, 
                 family = quasipoisson(link = "log"))
multiplot(children1, children2)


##其他廣義線性模型==================================================

#glm還支援伽瑪(Gamma)、反高斯（inverse gaussian）、準二項（quasibinomial）回歸
#對他們使用不同的連結函數（link functions）
#Gamma：inverse、identity、log
#Poisson：log、identity、sqrt
#inverse gaussian：1/mu^2、inverse、identity、log

#多項回歸模型：可以用來對幾種類別進行分類，可以執行好幾個羅吉斯回歸得到相同結果，或使用nnet套件的polr()函數或multinom()函數


##倖存分析==================================================

require(survival)
head(bladder) #stop（事件發生或病人離開研究的時間）、event（在該時間是否有發生事件）

bladder[100:105, ]

survObject <- with(bladder[100:105, ], Surv(stop, event)) #Surv()建立反應變數
survObject
survObject[, 1:2] #首兩列有事件發生，發生時間為12，最後兩列沒事件發生，但可能在該時間後發生，因此該資料的發生時間被設限了（censored）


#倖存分析最常使用的模型為Cox比例風險模型（Proportional Hazard Model）

cox1 <- coxph(Surv(stop, event) ~ rx + number + size + enum, data = bladder)
summary(cox1)

#倖存曲線顯示「在某個時間點有多少比例的受試者存活」

plot(survfit(cox1), xlab = "Days", ylab= "Survival Rate", conf.int = TRUE)


#此處的rx變數是病人接受治療或安慰劑的指標，可傳遞到strata機贓料分成兩群來分析，產生兩條倖存曲線

cox2 <- coxph(Surv(stop, event) ~ strata(rx) + number + size + enum, data = bladder)
summary(cox2)

plot(survfit(cox2), xlab = "Days", ylab = "Survival Rate", onf.int = TRUE, col = 1:2)
legend("bottomleft", legend = c(1, 2), lty = 1, col = 1:2, text.col = 1:2, title = "rx")

cox.zph(cox1) #檢測比例風險模型的假設
cox.zph(cox2)


#Andersen-Gill分析和倖存分析相似，但處理的是區間資料，而且可以處理多個事件
#例如不只可以處理一間急診室是否有人求診，還能計算出急診室求診個數
#同樣用coxph()，要加一個附加變數到Surv，且必須根據用來識別資料的欄位（id）對資料分群

head(bladder2)
ag1 <- coxph(Surv(start, stop, event) ~ rx + number + size + enum + cluster(id), data = bladder2)

ag2 <- coxph(Surv(start, stop, event) ~ strata(rx) + number + size + enum + cluster(id), data = bladder2)

plot(survfit(ag1), conf.int = TRUE)
plot(survfit(ag2), conf.int = TRUE, col = 1:2)
legend("topright", legend = c(1, 2), lty = 1, col = 1:2, text.col = 1:2, title = "rx")


#Practical Machine Learning==================================================

#Predicting with regression==================================================

library(caret)
data(faithful)
set.seed(333)

inTrain <- createDataPartition(y = faithful$waiting, p = .5, list = FALSE)
trainFaith <- faithful[inTrain, ]
testFaith <- faithful[-inTrain, ]
head(trainFaith)

plot(trainFaith$waiting, trainFaith$eruptions, pch = 19, col = "blue", xlab = "Waiting", ylab = "Duration")

lm1 <- lm(eruptions ~ waiting, data = trainFaith)
summary(lm1)

plot(trainFaith$waiting, trainFaith$eruptions, pch = 19, col = "blue", xlab = "Waiting", ylab = "Duration")
lines(trainFaith$waiting, lm1$fitted, lwd = 3)

coef(lm1)[1] + coef(lm1)[2] *80
newdata <- data.frame(waiting = 80)
predict(lm1, newdata)

par(mfrow = c(1, 2))
plot(trainFaith$waiting, trainFaith$eruptions, pch = 19, col = "blue", xlab = "Waiting", ylab = "Duration")
lines(trainFaith$waiting, predict(lm1), lwd = 3)
plot(testFaith$waiting, testFaith$eruptions, pch = 19, col = "blue", xlab = "Waiting", ylab = "Duration")
lines(testFaith$waiting, predict(lm1, newdata = testFaith), lwd = 3)

sqrt(sum((lm1$fitted - trainFaith$eruptions)^2)) #RMSE
sqrt(sum((predict(lm1, newdata = testFaith) - testFaith$eruptions)^2))

pred1 <- predict(lm1, newdata = testFaith, interval = "prediction")
ord <- order(testFaith$waiting)
plot(testFaith$waiting, testFaith$eruptions, pch = 19, col = "blue")
matlines(testFaith$waiting[ord], pred1[ord, ], type = "l", col = c(1, 2, 2), lty = c(1, 1, 1), lwd = 3) #Prediction intervals預測區間

modFit <- train(eruptions ~ waiting, data = trainFaith, method = "lm")
summary(modFit$finalModel)


#Predicting with regression, mutiple covariates==================================================

library(ISLR)
library(ggplot2)
library(caret)

summary(Wage)

inTrain <- createDataPartition(y = Wage$wage, p = .7, list = FALSE)
training <- Wage[inTrain, ]
testing <- Wage[-inTrain, ]
dim(training); dim(testing)

featurePlot(x = training[, c("age", "education", "jobclass")], 
            y = training$wage,
            plot = "pairs")

qplot(age, wage, data = training)
qplot(age, wage, colour = jobclass, data = training)
qplot(age, wage, colour = education, data = training)

modFit <- train(wage ~ age + jobclass + education, method = "lm", data = training)
finMod <- modFit$finalModel
print(modFit)

plot(finMod, 1, pch = 19, cex = .5, col = "#00000010")

qplot(finMod$fitted, finMod$residuals, colour = race, data = training)

plot(finMod$residuals, pch = 19)

pred <- predict(modFit, testing)
qplot(wage, pred, colour = year, data = testing)

modFitAll <-train(wage ~., data = training, method = "lm")
pred <- predict(modFitAll, testing)
qplot(wage, pred, data = testing)


#Regularized Regression==================================================

library(ElemStatLearn)
data(prostate)
str(prostate)

small <- prostate[1:5, ]
lm(lpsa ~., data = small)

#in caret methods are:
##ridge
##lasso
##relaxo


#R軟體資料分析基礎與應用==================================================

##非線性最小平方法==================================================

fileUrl <- "http://jaredlander.com/data/wifi.rdata"
download.file(fileUrl, destfile = "wifi.rdata")
load("wifi.rdata")
head(wifi)

require(ggplot2)
ggplot(wifi, aes(x = x, y = y, color = Distance)) + geom_point() + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(wifi$Distance))


#nls函數常用來計算非線性最小平方
wifiMod1 <- nls(Distance ~ sqrt((betaX - x) ^ 2 + (betaY - y) ^ 2), #指定用根號模型
                data = wifi, start = list(betaX = 50, betaY = 50)) #網格的中心作為起始值
summary(wifiMod1)

ggplot(wifi, aes(x = x, y = y, color = Distance)) + geom_point() + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = mean(wifi$Distance)) + 
        geom_point(data = as.data.frame(t(coef(wifiMod1))), aes(x = betaX, y = betaY), size = 5, color = "green")


##樣條（Splines）==================================================

#樣條讓一些非線性資料有較平滑的分佈，甚至可用來對新資料做預測
#樣條實際上在每個資料點都有專屬的轉換函數f，我們要找出f的極小化
#lambda越大，曲線越平滑

data(diamonds) #以不同的自由度做平滑
diaSpline1 <- smooth.spline(x = diamonds$carat, y = diamonds$price)
diaSpline2 <- smooth.spline(x = diamonds$carat, y = diamonds$price, df = 2)
diaSpline3 <- smooth.spline(x = diamonds$carat, y = diamonds$price, df = 10)
diaSpline4 <- smooth.spline(x = diamonds$carat, y = diamonds$price, df = 20)
diaSpline5 <- smooth.spline(x = diamonds$carat, y = diamonds$price, df = 50)
diaSpline6 <- smooth.spline(x = diamonds$carat, y = diamonds$price, df = 100)

get.spline.info <- function(object) {
        data.frame(x = object$x, y = object$y, df = object$df)
}

require(plyr)

splineDF <- ldply(list(diaSpline1, diaSpline2, diaSpline3, diaSpline4, diaSpline5, diaSpline6), get.spline.info)
head(splineDF)

g <- ggplot(diamonds, aes(x = carat, y = price)) + geom_point()
g + geom_line(data = splineDF, aes(x = x, y = y, color = factor(round(df, 0)), group = df)) + 
        scale_color_discrete("Degrees of \nFreedom")


#最好的樣條為三次自然樣條，因為它在切點可以製造平滑的轉折，並在輸入資料的端點後面製造線性的現象
#ns()函數

require(splines)
head(ns(diamonds$carat, df = 1))
head(ns(diamonds$carat, df = 2))
head(ns(diamonds$carat, df = 3))
head(ns(diamonds$carat, df = 4))

g <- ggplot(diamonds, aes(x = carat, y = price)) + geom_point()
g + stat_smooth(method = "lm", formula = y ~ ns(x, 6), color = "blue") #6個切點的三次自然樣條
g + stat_smooth(method = "lm", formula = y ~ ns(x, 3), color = "red") #3個切點的三次自然樣條


##廣義加性模型（GAMs）==================================================

creditNames <- c("Checking", "Duration", "CreditHistory", "Purpose", "CreditAmount", "Savings", "Employment", 
                 "InstallmentRate", "GenderMarital", "OtherDebtors", "YearsAtResidence", "RealEstate", "Age", 
                 "OtherInstallment", "Housing", "ExistingCredits", "Job", "NumLiable", "Phone", "Foreign", "Credit")

theURL <- "http://archive.ics.uci.edu/ml/machine-learning-databases/statlog/german/german.data"
credit <- read.table(theURL, sep = "", header = FALSE, col.names = creditNames, stringsAsFactors = FALSE)
head(credit)

head(credit[, c("CreditHistory", "Purpose", "Employment", "Credit")])
creditHistory <- c(A30 = "All Paid", A31 = "All Paid This Bank", A32 = "Up To Date", A33 = "Late Payment", A34 = "Critical Account")
purpose <- c(A40 = "car (new)", A41 = "car (used)", A42 = "furniture/equipment", A43 = "radio/television", A44 = "domestic appliances", A45 = "repairs",
             A46 = "education", A47 = "(vacation - does not exist?)", A48 = "retraining", A49 = "business", A410 = "others")
employment <- c(A71 = "unemployed", A72 = "< 1 year", A73 = "1 - 4 years", A74 = "4 - 7 years", A75 = "> = 7 years")
credit$CreditHistory <- creditHistory[credit$CreditHistory]
credit$Purpose <- purpose[credit$Purpose]
credit$Employment <- employment[credit$Employment]
credit$Credit <- ifelse(credit$Credit == 1, "Good", "Bad")
credit$Credit <- factor(credit$Credit, levels = c("Good", "Bad"))
head(credit[, c("CreditHistory", "Purpose", "Employment", "Credit")])

require(useful)
ggplot(credit, aes(x = CreditAmount, y = Credit)) + 
        geom_jitter(position = position_jitter(height = .2)) + 
        facet_grid(CreditHistory ~ Employment) + 
        xlab("Credit Amount") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
        scale_x_continuous(labels = multiple)

ggplot(credit, aes(x = CreditAmount, y = Age)) + 
        geom_point(aes(color = Credit)) + 
        facet_grid(CreditHistory ~ Employment) + 
        xlab("Credit Amount") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
        scale_x_continuous(labels = multiple)


#gam可以用無母數平滑函數，例如樣條和張量積（tensor product），將連續變數作轉換

creditGam <- gam(Credit ~ te(CreditAmount) + s(Age) + CreditHistory + Employment, data = credit, family = binomial(link = "logit"))
summary(creditGam)

plot(creditGam, select = 1, se = TRUE, shade = TRUE) #CrdeitAmount張量積的平滑曲線
plot(creditGam, select = 2, se = TRUE, shade = TRUE) #Age樣條的平滑曲線

#灰色陰影為平滑曲線的信賴區間