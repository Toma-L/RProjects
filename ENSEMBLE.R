#利用R語言打通大數據的經脈==========

#10整合學習==================================================

install.packages("adabag")
library(adabag)
bank <- read.csv("bank.csv", header = TRUE, sep = ";")
dim(bank)
head(bank)

sub <- sample(1:nrow(bank), round(nrow(bank)/4)) #建構訓練集與測試集
length(sub)
bank_train <- bank[-sub, ]
bank_test <- bank[sub, ]
dim(bank_train)
dim(bank_test)

install.packages("rpart")
library(rpart)

bag <- bagging(y~., bank_train, mfinal = 5)
names(bag)
summary(bag)

bag$formula
bag$trees[2] #5棵決策樹，總票數5票
bag$votes[105:115, ] #看各樣本的投票情況
bag$prob[105:115, ] #看各樣本屬於各類別的機率
bag$class[105:115] #模型對各樣本的預測類別
bag$samples[105:115, ] #模型bag中對第105~115個樣本在5個base learner中的抽樣情況
bag$importance #看各feature在分類過程中的相對重要性

bag1 <- bagging(y~., bank_train, mfinal = 5, control = rpart.control(maxdepth = 3)) #控制樹的複雜度
bag1$trees[2]

pre_bag <- predict(bag, bank_test) #對測試集進行預測
names(pre_bag)
pre_bag$vote[1:10, ]
pre_bag$prob[1:10, ]
pre_bag$class[1:10]
pre_bag$confusion #用混淆矩陣分析預測結果
pre_bag$error

sub_minor <- which(bank_test$y == "yes")
sub_major <- which(bank_test$y == "no")
length(sub_minor)
length(sub_major) #發現資料不平衡

err_bag <- sum(pre_bag$class != bank_test$y) / nrow(bank_test)
err_minor_bag <- sum(pre_bag$class[sub_minor] != bank_test$y[sub_minor]) / length(sub_minor)
err_major_bag <- sum(pre_bag$class[sub_major] != bank_test$y[sub_major]) / length(sub_major)
err_bag
err_minor_bag
err_major_bag #不平衡資料集用Adaboost處理具有一定優勢（但未必）

boo <- boosting(y~., bank_train, mfinal = 5)
pre_boo <- predict(boo, bank_test)
err_boo <- sum(pre_boo$class != bank_test$y) / nrow(bank_test)
err_minor_boo <- sum(pre_boo$class[sub_minor] != bank_test$y[sub_minor]) / length(sub_minor)
err_major_boo <- sum(pre_boo$class[sub_major] != bank_test$y[sub_major]) / length(sub_major)
err_boo
err_minor_boo 
err_major_boo


#11隨機森林==================================================

install.packages("randomForest")
library(randomForest)

set.seed(4)
data(mtcars)
mtcars.rf <- randomForest(mpg ~., data = mtcars, ntree = 1000, importance = TRUE)
importance(mtcars.rf) #分析模型中各個變數的重要性
importance(mtcars.rf, type = 1) #type決定用什麼衡量重要性，1是精度平均較少值，2是節點不純度，都是越大越重要


set.seed(1)
data(iris)
iris.rf <- randomForest(Species ~., iris, proximity = TRUE)
MDSplot(iris.rf, iris$Species, pallete = rep(1, 3), pch = as.numeric(iris$Species)) #對隨機森林模型進行視覺化分析


data(iris) #對遺失值做內插是隨機森林的重要用途
iris.na <- iris
iris.na[75, 2] = NA
iris.na[125, 3] = NA
set.seed(111)
iris.imputed <- rfImpute(Species ~., data = iris.na)

list("real" = iris[c(75, 125), 1:4], "have-NA" = iris.na[c(75, 125), 1:4], "disposed" = round(iris.imputed[c(75, 125), 2:5], 1)) #比較實際值與內插值


iris.rf <- randomForest(Species ~., iris)
hist(treesize(iris.rf)) #treesize()用來檢視隨機森林模型中每一棵樹所具有的節點個數


data(airquality)
set.seed(131)
ozone.rf <- randomForest(Ozone ~., data = airquality, mtry = 3, importance = TRUE, na.action = na.omit)
plot(ozone.rf) #繪製誤差與決策樹數量的關係圖


library(randomForest)
url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
wine <- read.csv(url, header = TRUE, sep = ";")
names(wine)
summary(wine)
dim(wine)

cha = 0
for(i in 1:4898) {
        if(wine[i, 12] > 6) cha[i] = "good"
        else if(wine[i, 12] > 5) cha[i] = "mid"
        else cha[i] = "bad"
}
wine[, 12] = factor(cha)
summary(wine$quality)


set.seed(71)
samp <- sample(1:4898, 3000)
set.seed(111)
wine.rf <- randomForest(quality ~., data = wine, importance = TRUE, proximity = TRUE, ntree = 500, subset = samp) #第一種方法


x <- subset(wine, select = -quality)
y <- wine$quality
set.seed(71)
samp = sample(1:4898, 3000)
xr <- x[samp, ]
yr <- y[samp]
set.seed(111)
wine.rf <- randomForest(xr, yr, importance = TRUE, proximity = TRUE, ntree = 500)

print(wine.rf)
importance(wine.rf)


n <- ncol(wine) - 1
rate = 1
for(i in 1:n) {
        set.seed(222)
        model <- randomForest(quality ~., data = wine, mtry = i, importance = TRUE, ntree = 1000)
        rate[i] <- mean(model$err.rate)
        print(model)
}
rate


set.seed(222)
model <- randomForest(quality ~., data = wine, mtry = , importance = TRUE, ntree = 1000)

plot(model, col = 1:1)
legend(800, .215, "mid", cex = .9, bty = "n")
legend(800, .28, "bad", cex = .9, bty = "n")
legend(800, .37, "good", cex = .9, bty = "n")
legend(800, .245, "total", cex = .9, bty = "n")


set.seed(222)
model <- randomForest(quality ~., data = wine, mtry = 1, proximity = TRUE, importance = TRUE, ntree = 400)
print(model)
hist(treesize(model))
MDSplot(model, wine$quality, palette = rep(1, 3), pch = as.numeric(wine$quality))


#應用R語言於資料分析==========

rm(list = ls())
gc()

library(randomForest)
data(iris)
ind <- sample(2, nrow(iris), replace = TRUE, prob = c(.8, .2))
trainData <- iris[ind == 1, ] #建立訓練集
testData <- iris[ind == 2, ]
rf <- randomForest(Species ~., data = trainData, ntree = 100)

irisPred <- predict(rf, newdata = testData)
table(irisPred, testData$Species)


library(adabag) #AdaBoost
data(iris)

ind <- sample(2, nrow(iris), replace = TRUE, prob = c(.8, .2))
trainData <- iris[ind == 1, ]
testData <- iris[ind == 2, ]

train.adaboost <- boosting(Species ~., data = trainData, boos = TRUE, mfinal = 5)

test.adaboost.pred <- predict.boosting(train.adaboost, newdata = testData)
test.adaboost.pred$confusion


#Practical Machine Learning==================================================

#Bagging==================================================

library(ElemStatLearn)
data(ozone, package = "ElemStatLearn")
ozone <- ozone[order(ozone$ozone), ]
head(ozone)

l1 <- matrix(NA, nrow = 10, ncol = 155)
for(i in 1:10) {
        ss <- sample(1:dim(ozone)[1], replace = TRUE)
        ozone0 <- ozone[ss, ]
        ozone0 <- ozone0[order(ozone0$ozone), ]
        loess0 <- loess(temperature ~ ozone, data = ozone0, span = .2)
        l1[i, ] <- predict(loess0, newdata = data.frame(ozone = 1:155))
}

plot(ozone$ozone, ozone$temperature, pch = 19, cex = .5)
for(i in 1:10) {
        lines(1:155, l1[i, ], col = "grey", lwd = 2)
}
lines(1:155, apply(l1, 2, mean), col = "red", lwd = 2)

#other options: bagEarth, treebag, bagFDA

predictors <- data.frame(ozone = ozone$ozone)
temperature <- ozone$temperature
treebag <- bag(predictors, temperature, B = 10, bagControl = bagControl(fit = ctreeBag$fit,
                                                                        predict = ctreeBag$pred,
                                                                        aggregate = ctreeBag$aggregate))

plot(ozone$ozone, temperature, col = "lightgrey", pch = 19)
points(ozone$ozone, predict(treebag$fits[[1]]$fit, predictors), pch = 19, col = "red")
points(ozone$ozone, predict(treebag, predictors), pch = 19, col = "blue")

ctreeBag$fit
ctreeBag$pred
ctreeBag$aggregate


#Random forests==================================================

data(iris)
library(ggplot2)
inTrain <- createDataPartition(y = iris$Species, p = .7, list = FALSE)

training <- iris[inTrain, ]
testing <- iris[-inTrain, ]

library(caret)
modFit <- train(Species ~., data = training, method = "rf", prox = TRUE)
modFit

getTree(modFit$finalModel, k = 2)

irisP <- classCenter(training[, c(3, 4)], training$Species, modFit$finalModel$prox) #各類別的中心
irisP <- as.data.frame(irisP)
irisP$Species <- rownames(irisP)
p <- qplot(Petal.Width, Petal.Length, col = Species, data = training)
p + geom_point(aes(x = Petal.Width, y = Petal.Length, col = Species), size = 5, shape = 4, data = irisP)

pred <- predict(modFit, testing)
testing$predRight <- pred == testing$Species
table(pred, testing$Species)

qplot(Petal.Width, Petal.Length, colour = predRight, data = testing, main = "newdata Predictions")

?rfcv #避免overfitting, 做random forest可用cross validation


#Boosting==================================================

#四種boosting
##gbm: boosting with trees
##mboost: model based boosting
##ada: statistical boosting based on additive logistic regression
##gamBoost: for boosting generalized additive models

library(ISLR)
data(Wage)
library(ggplot2)
library(caret)

Wage <- subset(Wage, select = -c(logwage))
inTrain <- createDataPartition(y = Wage$wage, p = .7, list = FALSE)
training <- Wage[inTrain, ]
testing <- Wage[-inTrain, ]
modFit <- train(wage ~ ., method = "gbm", data = training, verbose = FALSE)
print(modFit)

qplot(predict(modFit, testing), wage, data = testing)


#Combining predictors==================================================

library(ISLR)
data(Wage)
library(ggplot2)
library(caret)
Wage <- subset(Wage, select = -c(logwage))

inBuild <- createDataPartition(y = Wage$wage, p = .7, list = FALSE)
validation <- Wage[-inBuild, ]
buildData <- Wage[inBuild, ]
inTrain <- createDataPartition(y = buildData$wage, p = .7, list = FALSE)
training <- buildData[inTrain, ]
testing <- buildData[-inTrain, ]

dim(training); dim(testing); dim(validation)

mod1 <- train(wage ~., method = "glm", data = training)
mod2 <- train(wage ~., method = "rf", data = training, trControl = trainControl(method = "cv"), number = 3)

pred1 <- predict(mod1, testing)
pred2 <- predict(mod2, testing)

qplot(pred1, pred2, colour = wage, data = testing)

predDF <- data.frame(pred1, pred2, wage = testing$wage)
combModFit <- train(wage ~ ., method = "gam", data = predDF)
combPred <- predict(combModFit, predDF)

sqrt(sum((pred1 - testing$wage) ^ 2))
sqrt(sum((pred2 - testing$wage) ^2 ))
sqrt(sum((combPred - testing$wage) ^ 2))

pred1V <- predict(mod1, validation)
pred2V <- predict(mod2, validation)
predVDF <- data.frame(pred1 = pred1V, pred2 = pred2V)
combPredV <- predict(combModFit, predVDF)

sqrt(sum((pred1V - validation$wage) ^ 2))
sqrt(sum((pred2V - validation$wage) ^ 2))
sqrt(sum((combPredV - validation$wage) ^ 2))


#R軟體資料分析基礎與應用==================================================

##隨機森林（Random Forest）==================================================

require(useful)
require(randomForest)
creditFormula <- Credit ~ CreditHistory + Purpose + Employment + Duration + Age + CreditAmount
creditX <- build.x(creditFormula, data = credit)
creditY <- build.y(creditFormula, data = credit)
creditForest <- randomForest(x = creditX, y = creditY)
creditForest


# R語言與數據挖掘最佳實踐和經典案例 =====

## 4 - 決策樹與隨機森林 =====

### 4.3 randomForest套件 =====

ind <- sample(2, nrow(iris), replace = TRUE, prob = c(.7, .3))
trainData <- iris[ind == 1, ]
testData <- iris[ind == 2, ]

library(randomForest)
rf <- randomForest(Species ~ ., data = trainData, ntree = 100, proximity = TRUE)
table(predict(rf), trainData$Species)
print(rf)
attributes(rf)
plot(rf)
importance(rf) # 變數的重要性
varImpPlot(rf) # 變數的重要性（圖示）

irisPred <- predict(rf, newdata = testData)
table(irisPred, testData$Species)

plot(margin(rf, testData$Species)) # 查看預測結果