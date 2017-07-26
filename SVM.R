#利用R語言打通大數據的經脈==========

# install.packages("e1071")
library(e1071)

data(iris)
summary(iris)

data(iris)
model = svm(Species~., data = iris)

x <- iris[, -5]
y <- iris[, 5]
model = svm(x, y, kernel = "radial", gamma = if(is.vector(x)) 1 else 1/ncol(x))

summary(model)

x <- iris[, 1:4]
pred <- predict(model, x)
pred[sample(1:150, 8)]

table(pred, y)

attach(iris)
x <- subset(iris, select = -Species)
y <- Species
type = c("C-classification", "nu-classification", "one-classification")
kernel = c("linear", "polynomial", "radial", "sigmoid")
pred = array(0, dim = c(150, 3, 4))
accuracy = matrix(0, 3, 4)
yy <- as.integer(y)

for(i in 1:3)
{
        for(j in 1:4)
        {
                pred[, i, j] = predict(svm(x, y, type = type[i], kernel = kernel[j]), x)
                if (i > 2) accuracy[i, j] = sum(pred[, i, j] != 1)
                else accuracy[i, j] = sum(pred[, i, j] != yy)
        }
}

dimnames(accuracy) <- list(type,kernel)
table(pred[, 1, 3], y)

plot(cmdscale(dist(iris[, -5])), 
     col = c("lightgray", "black", "gray")[as.integer(iris[, 5])],
     pch = c("o", "+")[1:150 %in% model$index + 1])
legend(2, -.8, c("setosa", "versicolor", "virginica"), 
       col = c("lightgray", "black", "gray"), lty = 1)

data(iris)
model = svm(Species~., data = iris)
plot(model, iris, Petal.Width ~ Petal.Length, fill = FALSE, symbolPalette = c("lightgray", "black", "gray"), svSymbol = "+")
legend(1, 2.5, c("setosa", "versicolor", "virginica"), col = c("lightgray", "black", "gray"), lty = 1)

wts <- c(1, 1, 1)
names(wts) = c("setosa", "versicolor", "virginica")
model1 <- svm(x, y, class.weights = wts)

wts <- c(1, 100, 100)
names(wts) = c("setosa", "versicolor", "virginica")
model2 = svm(x, y, class.weights = wts)
pred2 <- predict(model2, x)
table(pred2, y)

wts <- c(1, 500, 500)
names(wts) <- c("setosa", "versicolor", "virginica")
model3 <- svm(x, y, class.weights = wts)
pred3 <- predict(model3, x)
table(pred3, y)