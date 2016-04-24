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

