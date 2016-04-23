#Practical Machine Learning==================================================

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


#Preprocessing=====

