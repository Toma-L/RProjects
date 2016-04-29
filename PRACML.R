#Practical Machine Learning==================================================

#Components of a predictor:
##question ---> input data ---> features ---> algorithm ---> parameters ---> evaluation
##make your question CONCRETE

#Relative order of importance:
##question > data > features > algorithms

#Garbage in, garbage out.

#Features matter!
#Properties of good features:
##Lead to data compression
##Retain relevant information
##Are created based on expert application knowledge

#Common mistakes:
##Trying to automate feature selection
##Not paying attention to data-specific quirks
##Throwing away information unnecessarily


#What is prediction?==================================================

library(kernlab)
data(spam)
head(spam)

plot(density(spam$your[spam$type == "nonspam"]), col = "blue", main = "", 
     xlab = "Frequency of 'your'")
lines(density(spam$your[spam$type == "spam"]), col = "red")

plot(density(spam$your[spam$type == "nonspam"]), col = "blue", main = "",
     xlab = "Frequency of your'")
lines(density(spam$your[spam$type == "spam"]), col = "red")
abline(v = .5, col = "black")

prediction <- ifelse(spam$your > .5, "spam", "nonspam")
table(prediction, spam$type) / length(spam$type)


#Relative order of importance==================================================

#BEST ML Method
##Interpretable
##Simple
##Fast
##Accurate
##Scalable


#In sample and out of sample error==================================================

library(kernlab)
data(spam)
set.seed(333)
smallSpam <- spam[sample(dim(spam)[1], size = 10), ]
spamLabel <- (smallSpam$type == "spam") * 1 + 1
plot(smallSpam$capitalAve, col = spamLabel)

rule1 <- function(x) {
        prediction <- rep(NA, length(x))
        prediction[x > 2.7] <- "spam"
        prediction[x < 2.40] <- "nonspam"
        prediction[(x >= 2.40 & x <= 2.45)] <- "spam"
        prediction[(x > 2.45 & x <= 2.70)] <- "nonspam"
        return(prediction)
}
table(rule1(smallSpam$capitalAve), smallSpam$type)

rule2 <- function(x) {
        prediction <- rep(NA, length(x))
        prediction[x > 2.8] <- "spam"
        prediction[x <= 2.8] <- "nonspam"
        return(prediction)
}
table(rule2(smallSpam$capitalAve), smallSpam$type)

table(rule1(spam$capitalAve), spam$type)
table(rule2(spam$capitalAve), spam$type)
mean(rule1(spam$capitalAve) == spam$type)

sum(rule1(spam$capitalAve) == spam$type)
sum(rule2(spam$capitalAve) == spam$type)


#Types of errors:

#Sensitivity: TP / (TP + FN)
#Specificity: TN / (FP + TN)
#Positive Predictive Value: TP / (TP + FP)
#Negative Predictive Value: TN / (FN + TN)
#Accuracy: (TP + TN) / (TP + FP + FN + TN)


#Common error measures

#Mean squared error (or root mean squared error): Continuous data, sensitive to outliers
#Median absolute deviation: Continuous data, often more robust
#Sensitivity (Recall): If you want few missed positives
#Specificity: If you want few negatives called positives
#Accuracy: Weights false positives/negatives equally
#Concordance: ex: kappa


#The caret package==================================================

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
dim(training)

set.seed(32343)
modelFit <- train(type ~., data = training, method = "glm")
modelFit

modelFit <- train(type ~., data = training, method = "glm")
modelFit$finalModel

predictions <- predict(modelFit, newdata = testing)
predictions

confusionMatrix(predictions, testing$type)


#Data Slicing==================================================

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
dim(training)

set.seed(32323)
folds <- createFolds(y = spam$type, k = 10, list = TRUE, returnTrain = TRUE)
sapply(folds, length)

folds[[1]][1:10]

set.seed(32323)
folds <- createFolds(y = spam$type, k = 10, list = TRUE, returnTrain = FALSE) #return Test
sapply(folds, length)

set.seed(32323)
folds <- createResample(y = spam$type, times = 10, list = TRUE) #製造boostrap的sample
sapply(folds, length)

set.seed(32323)
tme <- 1:1000
folds <- createTimeSlices(y = tme, initialWindo = 20, horizon = 10) #Time Slices
names(folds)


#Training Options=====

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
modelFit <- train(type ~., data = training, method = "glm")

args(train.default)

#Continuous outcomes:
##RMSE: root mean squared error
##RSquared: R2 from regression models

#Categorical outcomes:
##Accuracy: Fraction correct
##Kappa: A measure of concordance

args(trainControl)

#method
##boot: bootstrapping
##boot632: bootstrapping with adjustment
##cv: cross validation
##repeatedcv: repeated cross validation
##LOOCV: leave one out cross validation

set.seed(1235)
modelFit2 <- train(type ~., data = training, method = "glm")
modelFit2

set.seed(1235)
modelFit3 <- train(type ~., data = training, method = "glm")
modelFit3


#Plotting Predictors=====

library(ISLR)
library(ggplot2)
library(caret)
data(Wage)
summary(Wage)

inTrain <- createDataPartition(y = Wage$wage, p = .7, list = FALSE)
training <- Wage[inTrain, ]
testing <- Wage[-inTrain, ]
dim(training)
dim(testing)

featurePlot(x = training[, c("age", "education", "jobclass")], y = training$wage, plot = "pairs")
qplot(age, wage, data = training)
qplot(age, wage, colour = jobclass, data = training)

qq <- qplot(age, wage, colour = education, data = training)
qq + geom_smooth(method = 'lm', formula = y ~ x)

library(Hmisc)
cutWage <- cut2(training$wage, g = 3)
table(cutWage)

p1 <- qplot(cutWage, age, data = training, fill = cutWage, geom = c("boxplot"))
p1

install.packages("gridExtra")
library(gridExtra)

p2 <- qplot(cutWage, age, data = training, fill = cutWage, geom = c("boxplot", "jitter"))
grid.arrange(p1, p2, ncol = 2)

t1 <- table(cutWage, training$jobclass)
t1

qplot(wage, colour = education, data = training, geom = "density")


#Preprocessing==================================================

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
hist(training$capitalAve, main = "", xlab = "ave. capital run length") #間隔太大

mean(training$capitalAve)
sd(training$capitalAve)

trainCapAve <- training$capitalAve
trainCapAveS <- (trainCapAve - mean(trainCapAve)) / sd(trainCapAve)
mean(trainCapAveS)
sd(trainCapAveS)

testCapAve <- testing$capitalAve
testCapAveS <- (testCapAve - mean(trainCapAve)) / sd(trainCapAve) #必須用training set得到的參數估計
mean(testCapAveS)
sd(testCapAveS)

preObj <- preProcess(training[, -58], method = c("center", "scale")) #preProcess()
trainCapAveS <- predict(preObj, training[, -58])$capitalAve
mean(trainCapAveS)
sd(trainCapAveS)

testCapAveS <- predict(preObj, testing[, -58])$capitalAve
mean(testCapAveS)
sd(testCapAveS)

set.seed(32343)
modelFit <- train(type ~., data = training, preProcess = c("center", "scale"), method = "glm")
modelFit

preObj <- preProcess(training[, -58], method = c("BoxCox")) #Box-Cox轉換使連續資料符合常態分配
trainCapAveS <- predict(preObj, training[, -58])$capitalAve
par(mfrow = c(1, 2))
hist(trainCapAveS)
qqnorm(trainCapAveS)


set.seed(13343)
training$capAve <- training$capitalAve
selectNA <- rbinom(dim(training)[1], size = 1, prob = .05) == 1
training$capAve[selectNA] <- NA #製造遺漏值

preObj <- preProcess(training[, -58], method = "knnImpute") #kNN法插補遺漏值
capAve <- predict(preObj, training[, -58])$capAve

capAveTruth <- training$capitalAve
capAveTruth <- (capAveTruth - mean(capAveTruth)) / sd(capAveTruth)

quantile(capAve - capAveTruth)
quantile((capAve - capAveTruth)[selectNA])
quantile((capAve - capAveTruth)[!selectNA])

#千萬不能用test set估計參數


#Covariate creation==================================================

library(ISLR)
library(caret)
data(Wage)

inTrain <- createDataPartition(y = Wage$wage, p = .7, list = FALSE)
training <- Wage[inTrain, ]
testing <- Wage[-inTrain, ]

table(training$jobclass)
dummies <- dummyVars(wage ~ jobclass, data = training) #dummyVars()創造虛擬變數
head(predict(dummies, newdata = training))

nsv <- nearZeroVar(training, saveMetrics = TRUE) #nearZeroVar()移除值幾乎都是0的變數
nsv

library(splines)
bsBasis <- bs(training$age, df = 3) #bs()basis function創造多項式函數
bsBasis

lm1 <- lm(wage ~ bsBasis, data = training)
plot(training$age, training$wage, pch = 19, cex = .5)
points(training$age, predict(lm1, newdata = training), col = "red", pch = 19, cex = .5)

predict(bsBasis, age = testing$age) #用training set建立的bsBasis直接用在test set上


#Preprocessing with Principal Components Analysis (PCA)==================================================

library(caret)
library(kernlab)
data(spam)

inTrain <- createDataPartition(y = spam$type, p = .75, list = FALSE)
training <- spam[inTrain, ]
testing <- spam[-inTrain, ]
M <- abs(cor(training[, -58]))
diag(M) <- 0
which(M > .8, arr.ind = TRUE) #挑出cor > .8的

names(spam)[c(34, 32)]
plot(spam[, 34], spam[, 32])

X <- .71 * training$num415 + .71 * training$num857
Y <- .71 * training$num415 - .71 * training$num857
plot(X, Y)


#SVD法：解釋強、干擾大、壓抑小
#PCA法：解釋弱、干擾小、壓抑大

smallSpam <- spam[, c(34, 32)]
prComp <- prcomp(smallSpam)
plot(prComp$x[, 1], prComp$x[, 2])

prComp$rotation #PC1: .7081 * num415 + .7061 * num857

typeColor <- ((spam$type == "spam") * 1 + 1)
prComp <- prcomp(log10(spam[, -58] + 1))
plot(prComp$x[, 1], prComp$x[, 2], col = typeColor, xlab = "PC1", ylab = "PC2")


preProc <- preProcess(log10(spam[, -58] + 1), method = "pca", pcaComp = 2) #用preProcess()做PCA
spamPC <- predict(preProc, log10(spam[, -58] + 1))
plot(spamPC[, 1], spamPC[, 2], col = typeColor)

preProc <- preProcess(log10(training[, -58] + 1), method = "pca", pcaComp = 2)
trainPC <- predict(preProc, log10(training[, -58] + 1))
modelFit <- train(training$type ~., method = "glm", data = trainPC)

testPC <- predict(preProc, log10(testing[, -58] + 1))
confusionMatrix(testing$type, predict(modelFit, testPC))

modelFit <- train(training$type ~., method = "glm", preProcess = "pca", data = training)
confusionMatrix(testing$type, predict(modelFit, testing))

#主成份分析可能會使結果更難解釋
#最好先轉換資料（logs/Box Cox）
#可以把圖畫出來確認問題在哪


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


#Predicting with trees==================================================

data(iris)
library(ggplot2)
names(iris)
table(iris$Species)

inTrain <- createDataPartition(y = iris$Species, p = .7, list = FALSE)
training <- iris[inTrain, ]
testing <- iris[-inTrain, ]
dim(training); dim(testing)

qplot(Petal.Width, Sepal.Width, colour = Species, data = training)

library(caret)
modFit <- train(Species ~., method = "rpart", data = training)
print(modFit$finalModel)

plot(modFit$finalModel, uniform = TRUE, main = "Classification Tree")
text(modFit$finalModel, use.n = TRUE, all = TRUE, cex = .8)

library(rattle)
fancyRpartPlot(modFit$finalModel) #fancyRpartPlot()

predict(modFit, newdata = testing)

#other options: party, rpart


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


#Model based prediction==================================================

#假設資料follow某機率模型，用貝氏機率去指認出來

