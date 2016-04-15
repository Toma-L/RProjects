#Getting and Cleaning Data==================================================

#Using data.table

library(data.table)
DF <- data.frame(x = rnorm(9), y = rep(c("a", "b", "c"), each = 3), z = rnorm(9))
head(DF, 3)
DT <- data.table(x = rnorm(9), y = rep(c("a", "b", "c"), each = 3), z = rnorm(9))
head(DT, 3)
tables() #查看現在有哪些data.table
DT[2, ]
DT[DT$y == "a", ]
DT[c(2, 3)] #都是row
DT[, c(2, 3)] #執行c(2, 3)
{
        x = 1
        y = 2
}
k <- {print(10); 5} #10，會執行前面的
print(k) #5
DT[, list(mean(x), sum(z))]
DT[, table(y)]
DT[, w := z^2] #使用:=指定新的變數存在原本的資料裡
DT2 <- DT
DT[, y := 2] #y全變2
DT
head(DT, 3)
head(DT2, 3) #改變DT也會改變到DT2
DT[, m := {tmp <- (x+z); log2(tmp + 5)}]
head(DT)

DT[, a := x > 0]
head(DT)

DT[, b := mean(x + w), by = a] #把資料分成a為TRUE和a為FALSE，各自計算mean(x+w)，當成其b的值

set.seed(123)
DT <- data.table(x = sample(letters[1:3], 1E5, TRUE))
DT[, .N, by = x] #用.N計算次數

DT <- data.table(x = rep(c("a", "b", "c"), each = 100), y = rnorm(300))
setkey(DT, x)
DT['a']

DT1 <- data.table(x = c('a', 'a', 'b', 'dt1'), y = 1:4)
DT2 <- data.table(x = c('a', 'b', 'dt2'), z = 5:7)
setkey(DT1, x)
setkey(DT2, x)
merge(DT1, DT2) #只留下都有的

big_df <- data.frame(x = rnorm(1E6), y = rnorm(1E6))
file <- tempfile()
write.table(big_df, file = file, row.names = FALSE, col.names = TRUE, sep"\t", quote = FALSE)
system.time(fread(file)) #fread()類似read.table()
system.time(read.table(file, header = TRUE, sep = "\t"))

#Subsetting and sorting

set.seed(13435)
X <- data.frame("var1" = sample(1:5), "var2" = sample(6:10), "var3" = sample(11:15))
X <- X[sample(1:5), ]
X$var2[c(1, 3)] <- NA
X

X[, 1]
X[, "var1"]
X[1:2, "var2"]
X[(X$var1 <= 3 & X$var3 > 11), ]
X[(X$var1 <= 3 | X$var3 > 15), ]
X[which(X$var2 > 8), ]
sort(X$var1)
sort(X$var1, decreasing = TRUE)
sort(X$var2, na.last = TRUE)
X[order(X$var1), ]
X[order(X$var1, X$var3), ]

library(plyr)
arrange(X, var1)
arrange(X, desc(var1))

X$var4 <- rnorm(5)
X

Y <- cbind(X, rnorm(5))
Y

#Summarizing data

if(!file.exists("./data")) {dir.create("./data")}
fileUrl <- "http://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType = DOWNLOAD"
download.file(fileUrl, destfile = "./data/restaurants.csv")
restData <- read.csv("./data/restaurants.csv")
head(restData, n = 3)
tail(restData, n = 3)
summary(restData)
str(restData)

quantile(restData$councilDistrict, na.rm = TRUE)
table(restData$zipCode, useNA = "ifany") #如果有NA的話
table(restData$councilDistrict, restData$zipCode)

sum(is.na(restData$councilDistrict))
any(is.na(restData$councilDistrict))
all(restData$zipCode > 0)

colSums(is.na(restData))
all(colSums(is.na(restData)) == 0)

table(restData$zipCode %in% c("21212")) #%in%前面有多少符合後面
table(restData$zipCode %in% c("21212", "21213"))

restData[restData$zipCode %in% c("21212", "21213"), ]

data(UCBAdmissions)
DF <- as.data.frame(UCBAdmissions)
summary(DF)

xt <- xtabs(Freq ~ Gender + Admit, data = DF) #cross tabs

warpbreaks$replicate <- rep(1:9, len = 54)
xt <- xtabs(breaks ~., data = warpbreaks)
xt

ftable(xt)

fakeData <- rnorm(1e5)
object.size(fakeData) #資料大小

#Creating new variables

s1 <- seq(1, 10, by = 2); s1
s2 <- seq(1, 10, length = 3); s2
x <- c(1, 3, 8, 25, 100); seq(along = x)
restData$nearMe <- restData$neighborhood %in% c("Roland Park", "Homeland")
table(restData$nearMe)

restData$zipWrong <- ifelse(restData$zipCode < 0, TRUE, FALSE)
table(restData$zipWrong, restData$zipCode < 0)

restData$zipGroups <- cut(restData$zipCode, breaks = quantile(restData$zipCode)) #cut()
table(restData$zipGroups)
table(restData$zipGroups, restData$zipCode)

library(Hmisc)
restData$zipGroups <- cut2(restData$zipCode, g = 4) #easier cutting
table(restData$zipGroups)

restData$zcf <- factor(restData$zipCode)
restData$zcf[1:10]
class(restData$zcf)

yesno <- sample(c("yes", "no"), size = 10, replace = TRUE)
yesnofac <- factor(yesno, levels = c("yes", "no")) #levels為factor排序
relevel(yesnofac, ref = "yes")

as.numeric(yesnofac)

library(Hmisc)
restData$zipGroups <- cut2(restData$zipCode, g = 4)
table(restData$zipGroups)

library(Hmisc)
library(plyr)
restData2 <- mutate(restData, zipGroups = cut2(zipCode, g = 4)) #mutate()創造新變數
table(restData2$zipGroups)

#Reshaping data

library(reshape2)
head(mtcars)
mtcars$carname <- rownames(mtcars)
carMelt <- melt(mtcars, id = c("carname", "gear", "cyl"), measure.vars = c("mpg", "hp"))
head(carMelt, n = 3)
tail(carMelt, n = 3)

cylData <- dcast(carMelt, cyl ~ variable)
cylData

cylData <- dcast(carMelt, cyl ~ variable, mean) #對variable的值做mean
cylData

head(InsectSprays)
tapply(InsectSprays$count, InsectSprays$spray, sum)

spIns <- split(InsectSprays$count, InsectSprays$spray)
spIns
sprCount <- lapply(spIns, sum)
sprCount
unlist(sprCount)

sapply(spIns, sum)

ddply(InsectSprays, .(spray), summarize, sum = sum(count))

spraySums <- ddply(InsectSprays, .(spray), summarize, sum = ave(count, FUN = sum))
head(spraySums)

#Merging data

if(!file.exists("./data")) {dir.create("./data")}
fileUrl1 <- "http://dl.dropboxusercontent.com/u/7710864/data/reviews-apr29.csv"
fileUrl2 <- "http://dl.dropboxusercontent.com/u/7710864/data/solutions-apr29.csv"
download.file(fileUrl1, destfile = "./data/reviews.csv")
download.file(fileUrl2, destfile = "./data/solutions.csv")
reviews <- read.csv("./data/reviews.csv")
solutions <- read.csv("./data/solutions.csv")
head(reviews, 2)
head(solutions, 2)
names(reviews)
names(solutions)

mergedData <- merge(reviews, solutions, by.x = "solution_id", by.y = "id", all = TRUE)
head(mergedData)
intersect(names(solutions), names(reviews)) #intersect()找共同都有的變數

mergedData2 <- merge(reviews, solutions, all = TRUE) #不設定會自動merge名字一樣的變數
head(mergedData2)

df1 <- data.frame(id = sample(1:10), x = rnorm(10))
df2 <- data.frame(id = sample(1:10), y = rnorm(10))
arrange(join(df1, df2), id) #用join()合併

df1 <- data.frame(id = sample(1:10), x = rnorm(10))
df2 <- data.frame(id = sample(1:10), y = rnorm(10))
df3 <- data.frame(id = sample(1:10), z = rnorm(10))
dfList <- list(df1, df2, df3)
join_all(dfList) #很多個data.frame用join_all合併

#Editing text variables

cameraData <- read.csv("./data/cameras.csv")
names(cameraData)

splitNames <- strsplit(names(cameraData), "\\.")
splitNames[[5]]
splitNames[[6]]

mylist <- list(letters = c("A", "b", "c"), numbers = 1:3, matrix(1:25, ncol = 5))
head(mylist)
mylist[1]
mylist$letters
mylist[[1]]

splitNames[[6]][1]
firstElement <- function(x) {x[1]}
sapply(splitNames, firstElement)

names(reviews)
sub("_", "", names(reviews), ) #用""取代掉"_"，用sub()

testName <- "this_is_a_test"
sub("_", "", testName) #sub只會作用在第1個

gsub("_", "", testName) #gsub()作用在全部

grep("Alameda", cameraData$intersection) #grep()
table(grepl("Alameda", cameraData$intersection)) #grepl()變布林值

grep("Alameda", cameraData$intersection, value = TRUE) #直接給值，不要編號
grep("JeffStreet", cameraData$intersection)
length(grep("JeffStreet", cameraData$intersection))

library(stringr)
nchar("Jeffrey Leek") #計算字串長度
substr("Jeffrey Leek", 1, 7) #只取字串的某部分
paste("Jeffrey", "Leek")
paste0("Jeffrey", "Leek") #paste0()不會有空白
str_trim("Jeff       ") #清除空白

#Regular Expressions

## ^i think：以i think為首的
## morning$：以morning作結的
## [Bb][Uu][Ss][Hh]：不論大小寫
## ^[Ii] am
## ^[0-9][a-zA-Z]
## [^?.]$：[]放開頭代表排除

#Regular Expressions II

## 9.11：.代表什麼字都可以
## flood|fire：|表示「或」
## flood|earthquake|hurricane|coldfire
## ^[Gg]ood|[Bb]ad：bad或Bad就未必在句首
## ^([Gg]ood|[Bb]ad)：good和bad都只要句首
## [Gg]eorge([Ww]\.)? [Bb]ush：?代表「未必要」這個
## [Gg]eorge([Ww]\.)? [Bb]ush：加\代表.是純文字
## (.*)：.代表()裡面可以是任何字，*代表「不限字數」
## [0-9]+(.*)[0-9]+：+代表「至少要有一個」
## [Bb]ush( +[^ ]+ +){1,5}：{1,5}代表重複1~5次

## +([a-zA-Z]+) +\1 +
## ^s(.*)s
## ^s(.*?)s$

#Working with dates

d1 <- date()
d1
class(d1) #character

d2 <- Sys.Date()
d2
class(d2) #Date

format(d2, "%a %b %d")

x <- c("1jan1960", "2jan1960", "31mar1960", "30jul1960")
z <- as.Date(x, "%d%b%Y")
z
z[1] - z[2]
as.numeric(z[1] - z[2])

weekdays(d2)
months(d2)
julian(d2)

library(lubridate)
ymd("20140108")
mdy("08/04/2013")
dmy("03-04-2013")
ymd_hms("2011-08-03 10:15:03")
ymd_hms("2011-08-03 10:15:03", tz = "Pacific/Auckland")
?Sys.timezone

x <- dmy(c("1jan2013", "2jan2013", "31mar2013", "30jul2013"))
wday(x[1])
wday(x[1], label = TRUE)


#實用R程式設計==================================================

##R繪圖==================================================

data(cars)
str(cars)

###繪圖視窗之設定

windows(width = 4.5, height = 3.3, pointsize = 8)
old.par <- par(mex = .8, mar = c(5, 5, 4, 2) + .1) #mex是邊界文字縮放比
plot(cars)
par(old.par)

win.graph(width = 4.5, height = 3.3, pointsize = 8)
old.par <- par(mex = .8, mar = c(5, 5, 4, 2) + .1)
plot(cars)
par(old.par)

windows()
win.graph()

###常用的圖形參數

windows(width = 4.5, height = 3.3, pointsize = 8)
old.par <- par(mex = .8, mar = c(7, 5, 4, 2) + .1)
plot(cars, xlim = c(0, 30), ylim = c(0, 130), xlab = "xlab",
     ylab = "ylab", main = "main title", sub = "subtitle", 
     cex = .8, pch = 16, col = "red") #cex為文字及符號相對於內定值之縮放比
par(old.par)


##探索資料圖形==================================================

data(iris)
str(iris)
summary(iris)
head(iris, 5)
tail(iris, 5)
iris[iris$Species == "setosa", 1]
subset(iris, Species == "setosa", select = Sepal.Length)
iris[iris$Species == "setosa", 1:2]
subset(iris, Species == "setosa", select = c(Sepal.Length, Sepal.Width))

x <- iris[, 1]
c(min(x), max(x))
range(x)
summary(x)
fivenum(x)
mean(x)
median(x)
IQR(x)
mean(x, trim = .1)
sd(x)
cv <- sd(x) / mean(x)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(4, 4, 3, 2) + .1)
boxplot(x)
rug(x, side = 4)
boxplot(x, horizontal = TRUE)
rug(x, side = 1)
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 5, 3, 1) + .1)
hist(x, freq = TRUE, breaks = "Sturges")
rug(x, side = 1)
hist(x, prob = TRUE, breaks = "Sturges", col = "lightblue", border = "magenta")
rug(x, side = 1)
par(old.par)

hist(x, prob = TRUE, breaks = "Sturges")
hist(x, prob = TRUE, breaks = "Scott")
hist(x, prob = TRUE, breaks = "Freedman-Diaconis")

hist(x, prob = TRUE, breaks = 10)
hist(x, prbo = TRUE, breaks = seq(from = 4, to = 8, by = .25))

library(MASS)
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 5, 3, 1) + .1)
truehist(x, prob = FALSE, ylab = "Frequency", main = "Histogram")
truehist(x, prob = TRUE, ylab = "Density", main = "Histogram")
par(old.par)

stem(x, scale = .5)
sum(x == 4.4)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
stripchart(x, method = "overplot", at = .7)
text(6, .65, "overplot")
stripchart(x, method = "stack", add = TRUE, at = .85)
text(6, .8, "stack")
stripchart(x, method = "jitter", add = TRUE, at = 1.2)
text(6, 1.05, "jitter")
title(main = "strip chart")
par(old.par)

y <- cut(x, breaks = 6)
z <- table(y)
old.par <- par(mfrow = c(2, 2), mex = .2, mar = c(3, 3, 3, 2) + .1)
pie(z)
pie(z, clockwise = TRUE)
pie(z, col = terrain.colors(6))
pie(z, col = gray(seq(from = .4, to = 1.0, length = 6)))
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(density(x), col = "red", main = "Kernel density estimate")
rug(x, side = 1)
hist(x, prob = TRUE, breaks = "Sturges", main = "Histogram and KDE")
lines(density(x), col = "red")
rug(x, side = 1)
par(old.par)

f <- density(x)
class(f)
names(f)
print(f)
plot(f, type = "n")
polygon(f, col = "wheat")

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot.ecdf(x)
par(old.par)

F <- ecdf(x)
class(F)
names(F)
print(F)
summary(F)
plot(F)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
qqnorm(x)
qqline(x, col = "red", lwd = 2)
par(old.par)

data(VADeaths)
str(VADeaths)
class(VADeaths)
VADeaths
names <- c("RM", "RF", "UM", "UF")
colors <- c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk")
old.par <- par(mfrow = c(2, 2), mex = .8, mar = c(3, 3, 3, 2) + .1)
barplot(VADeaths, names.arg = names)
barplot(VADeaths, names.arg = names, horiz = TRUE)
barplot(VADeaths, names.arg = names, col = colors, border = "blue")
barplot(VADeaths, names.arg = names, col = colors, border = "blue", space = 1.5)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(VADeaths, beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan", "lavender", "cornsilk"), legend.text = rownames(VADeaths), ylim = c(0, 100))
title(main = "Death Rates in Virginia", font.main = 4)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(t(VADeaths), beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan", "lavender"), legend.text = rownames(t(VADeaths)), ylim = c(0, 80), args.legend = list(x = "topleft"))
title(main = "Death Rates in Virginia", font.main = 4)
par(old.par)

colnames(VADeaths) <- c("RM", "RF", "UM", "UF")
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
dotchart(VADeaths, xlim = c(0, 100), xlab = "Deaths per 1000", main = "Death rates")
dotchart(t(VADeaths), xlim = c(0, 100), xlab = "Deaths per 1000", main = "Death Rates")
par(old.par)

data(warpbreaks)
str(warpbreaks)
with(warpbreaks, tapply(breaks, INDEX = wool, FUN = sum))
with(warpbreaks, tapply(breaks, INDEX = tension, FUN = sum))
with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = sum))
xtabs(breaks ~ wool, data = warpbreaks)
xtabs(breaks ~ tension, data = warpbreaks)
xtabs(breaks ~ wool + tension, data = warpbreaks)
ftable(xtabs(breaks ~ wool + tension, data = warpbreaks))

t1 <- with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = sum))
t2 <- with(warpbreaks, tapply(breaks, INDEX = list(wool, tension), FUN = mean))
old.par <- par(mfrow = c(1, 2), mex = .8, mar = c(5, 4, 3, 1) + .1)
barplot(t1, beside = TRUE, col = c("lightblue", "mistyrose"), main = "counts", legend.text = rownames(t1), ylim = c(0, max(t1)))
barplot(t2, beside = TRUE, col = c("lightblue", "mistyrose"), main = "means", legend.text = rownames(t2), ylim = c(0, max(t2)))
par(old.par)

brks <- as.integer(xtabs(breaks ~ wool + tension, data = warpbreaks))
label <- c("AL", "BL", "AM", "BM", "AH", "BH")
old.par <- par(mfrow = c(1, 2), mex = .2, mar = c(3, 3, 3, 2) + .1)
pie(brks, label = label)
pie(brks, label = label, col = gray(seq(from = 0.4, to = 1.0, length = 8)))
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], panel = panel.smooth))
par(old.par)

panel.hist <- function(x, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(usr[1:2], 0, 1.5))
        h <- hist(x, plot = FALSE)
        breaks <- h$breaks
        nB <- length(breaks)
        y <- h$counts
        y <- y / max(y)
        rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
        lines(density(x, na.rm = TRUE), col = "red")
}
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], panel = panel.smooth, pch = 1, bg = "lightcyan", diag.panel = panel.hist, font.labels = 2, cex.labels = 1.2)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
pairs(iris[, 1:4], pch = c(1, 2, 4)[iris$Species], col = c("red", "green", "blue")[iris$Species])
par(old.par)

library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(Sepal.Length + Sepal.Width ~ Petal.Length + Petal.Width | Species, data = iris, layout = c(2, 2), scales = list(cex = .5, cex.lab = .5, relation = "free"), auto.key = TRUE)

setosa <- iris[iris$Species == "setosa", 1:4]
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
boxplot(setosa, names = c("sep.len", "sep.wid", "pet.len", "pet.wid"), main = "Iris setosa")
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 4, 2) + .1)
with(iris, boxplot(Sepal.Length ~ Species, main = "Sepal length"))
with(iris, boxplot(Sepal.Length ~ Species, notch = TRUE, main = "Sepal length"))
par(old.par)

sx <- with(iris, split(Sepal.Length, Species))
old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 4, 2), + .1)
boxplot(sx, main = "Sepal length")
boxplot(sx, notch = TRUE, main = "Sepal length")
par(old.par)

sx <- with(iris, split(Sepal.Length, Species))
sy <- with(iris, split(Sepal.Width, Species))

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(0, xlim = range(sx), ylim = range(sy), type = "n", xlab = "x", ylab = "y")
points(sx[[1]], sy[[1]], pch = 1, col = 1)
points(sx[[2]], sy[[2]], pch = 2, col = 2)
points(sx[[3]], sy[[3]], pch = 3, col = 3)
for(i in 1:3) abline(lm(sy[[i]] ~ sx[[i]]), col = i)
legend("topright", legend = c("setosa", "versicolor", "virginica"), lty = 1, pch = 1:3, col = 1:3)
par(old.par)

x <- iris[[1]]
y <- iris[[2]]
species <- iris[[5]]
library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(y ~ x, groups = species, type = c("g", "p", "r"), auto.key = TRUE)
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(y ~ x | species, type = c("g", "p", "r"), auto.key = TRUE)
par(old.par)

x <- iris[[1]]
y <- iris[[2]]

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
qqplot(x, y)
abline(a = 0, b = 1, col = "red")
par(old.par)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2)
        10 * sin(r) / r
}
x <- seq(from = -10, to = 10, length = 50)
y <- x
z <- outer(x, y, f)
z[is.na(z)] <- 1

old.par <- par(mfrow = c(1, 2), mex = .4, mar = c(5, 4, 3, 1) + .1, bg = "white")
persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightgreen")
persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightblue", ltheta = 120, shade = .75, ticktype = "detailed", xlab = "x", ylab = "y", zlab = "f(x, y)")
par(old.par)

persp(x, y, z, theta = 30, phi = 30, expand = .5, col = "lightgreen", box = FALSE)

library(lattice)
f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
y <- x <- seq(from = -10, to = 10, length = 50)
h <- expand.grid(x = x, y = y)
z <- f(h$x, h$y)
z[is.na(z)] <- 1
old.par <- par(mex = .8, mar = c(5, 5, 3, 1) + .1)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)")
par(old.par)

old.par <- par(mex = .8, mar = c(5, 5, 3, 1) + .1, cex.axis = .5)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), light.source = c(10, 0, 10), col = "lightblue")
par(old.par)

wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = TRUE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue")
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue", colorkey = TRUE)
wireframe(z ~ h$x * h$y, xlab = "x", ylab = "y", zlab = "f(x, y)", scales = list(arrows = FALSE), shade = TRUE, light.source = c(10, 0, 10), col = "lightblue", colorkey = TRUE, drape = TRUE, screen = list(z = 30, x = -60))

library(rgl)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
y <- x <- seq(from = -10, to = 10, length = 50)
z <- outer(x, y, f)
z[is.na(z)] <- 1

persp3d(x, y, z, xlab = "x", ylab = "y", zlab = "f(x, y)", col = "lightblue")

rgl.snapshot("/Users/thomas/Documents/R Files/persp1.png")
snapshot3d("/Users/thomas/Documents/R Files/persp2.png")

persp3d(x, y, z, xlab = "x", ylab = "y", zlab = "f(x, y)", col = "lightblue")
play3d(spin3d(axis = c(0, 0, 1), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 0), rpm = 6), duration = 10)
play3d(spin3d(axis = c(0, 1, 1), rpm = 8), duration = 5)

f <- function(x, y) {
        r <- sqrt(x^2 + y^2);
        10 * sin(r) / r
}
x <- seq(from = -10, to = 10, length = 100)
y <- x
z <- outer(x, y, f)
z[is.na(z)] <- 1

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 4, 3, 1) + .1)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, add = TRUE)
col <- gray(seq(from = .4, to = 1.0, length = 6))
image(x, y, z, col = col)
contour(x, y, z, add = TRUE)
par(old.par)

old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(6, 4, 3, 1) + .1)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, levels = c(-5, -2, 0, 2, 5), lwd = c(1, 1, 2, 1, 1), add = TRUE)
image(x, y, z, col = cm.colors(20))
contour(x, y, z, levels = 0, lwd = 2, add = TRUE)
par(old.par)

old.par <- par(mex = .8, mar = c(3, 3, 2, 1) + .1)
filled.contour(x, y, z)
par(old.par)

old.par <- par(mex = .8, mar = c(3, 3, 2, 1) + .1)
filled.contour(x, y, z, levels = c(-5, -2, 0, 2, 5), lwd = c(1, 1, 2, 1, 1))
par(old.par)

data(iris)
x <- iris[, 1]
y <- iris[, 2]
z <- iris[, 3]

library(lattice)
old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
cloud(z ~ x * y, groups = iris$Species, pch = 1:3, col = 1:3, scales = list(arrows = FALSE), light.source = c(10, 0, 10))
par(old.par)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
cloud(z ~ x * y, groups = iris$Species, pch = 1:3, col = 1:3, scales = list(arrows = FALSE), light.source = c(10, 0, 10), shade = TRUE, colorkey = TRUE)
par(old.par)

library(scatterplot3d)
old.par <- par(mfrow = c(1, 2), mex = .5, mar = c(5, 5, 3, 1) + .1)
scatterplot3d(x, y, z, xlab = "x", ylab = "y", zlab = "z", angle = 30, y.margin.add = .1, scale.y = .7, pch =c(1, 2, 3)[iris$Species], color = c("red", "green", "blue")[iris$Species])
scatterplot3d(x, y, z, xlab = "x", ylab = "y", zlab = "z", type = "h", angle = 30, y.margin.add = .1, scale.y = .7, pch = c(1, 2, 3)[iris$Species], color = c("red", "green", "blue")[iris$Species])
par(old.par)

library(rgl)
plot3d(x, y, z, pch = c(1, 2, 3)[iris$Species], col = c("red", "green", "blue")[iris$Species])
plot3d(x, y, z, pch = c(1, 2, 3)[iris$Species], col = c("red", "green", "blue")[iris$Species])
play3d(spin3d(axis = c(0, 0, 1), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 0), rpm = 8), duration = 10)
play3d(spin3d(axis = c(0, 1, 1), rpm = 8), duration = 10)

install.packages("SemiPar")
library(SemiPar)
data(pig.weights, package = "SemiPar")
str(pig.weights)
xtabs(weight ~ id.num + num.weeks, data = pig.weights)
ftable(xtabs(weight ~ id.num + num.weeks, data = pig.weights))

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
plot(weight ~ num.weeks, data = pig.weights)
for(i in 1:48) lines(weight ~ num.weeks, subset(pig.weights, id.num == i), col = i)
par(old.par)

library(lattice)
xyplot(weight ~ num.weeks, data = pig.weights, groups = id.num, type = "b")

library(nlme)
pig.growth <- groupedData(weight ~ num.weeks | id.num, data = pig.weights)
class(pig.growth)
plot(pig.growth, outer = ~1, key = FALSE)
plot(pig.growth, outer = ~1, key = TRUE)

old.par <- par(mex = .8, mar = c(5, 4, 3, 1) + .1)
xyplot(weight ~ num.weeks | id.num, data = pig.weights, type = c("b", "g", "p", "r"), scales = list(cex = .5, cex.axis = .5, cex.lab = .5), auto.key = TRUE)
par(old.par)

library(nlme)
pig.growth <- groupedData(weight ~ num.weeks | id.num, data = pig.weights)
class(pig.growth)
plot(pig.growth, pch = 16)


#R軟體資料分析基礎與應用==================================================

##07統計繪圖==================================================

require(ggplot2)
data(diamonds)
head(diamonds)
hist(diamonds$carat, main = "Carat Histogram", xlab = "Carat")
plot(price ~ carat, data = diamonds)
plot(diamonds$carat, diamonds$price)
identical(plot(price ~ carat, data = diamonds), plot(diamonds$carat, diamonds$price))
boxplot(diamonds$carat)

ggplot(data = diamonds) + geom_histogram(aes(x = carat))
ggplot(data = diamonds) + geom_point(aes(x = carat, y = price))
ggplot(data = diamonds, aes(x = carat, y = price)) + geom_point()
g <- ggplot(diamonds, aes(x = carat, y = price))
g + geom_point(aes(color = color))
g + geom_point(aes(color = color)) + facet_wrap(~color)
g + geom_point(aes(color = color)) + facet_grid(cut ~ clarity)
ggplot(diamonds, aes(x = carat)) + geom_histogram() + facet_wrap(~color)
ggplot(diamonds, aes(y = carat, x = 1)) + geom_boxplot()
ggplot(diamonds) + geom_boxplot(aes(y = carat, x = clarity))
ggplot(diamonds, aes(x = cut, y = carat)) + geom_boxplot()
ggplot(diamonds, aes(x = cut, y = carat)) + geom_violin()
ggplot(diamonds, aes(x = cut, y = carat)) + geom_violin() + geom_point()

ggplot(economics, aes(x = date, y = pop)) + geom_line()
head(economics)

require(lubridate)
economics$year <- year(economics$date)
economics$month <- month(economics$date)
head(economics)

econ2000 <- economics[which(economics$year >= 2000), ]

require(scales)
g <- ggplot(econ2000, aes(x = month, y = pop))
g <- g + geom_line(aes(color = factor(year), group = year))
g <- g + scale_color_discrete(name = "Year")
g <- g + scale_y_continuous(labels = comma)
g <- g + labs(title = "Population Growth", x = "Month", y = "Population")
g

g2 <- ggplot(diamonds, aes(x = carat, y = price)) + geom_point(aes(color = color))
g2
require(ggthemes)
g2 + theme_economist()
g2 + theme_economist() + scale_colour_economist()
g2 + theme_excel()
g2 + theme_excel() + scale_colour_excel()
g2 + theme_tufte()
g2 + theme_wsj()


#利用R語言打通大數據的經脈==================================================

##探索性資料分析==================================================

###資料集

library(MASS)
data(Insurance)
nrow(Insurance)
ncol(Insurance)
dim(Insurance)
head(Insurance)


###數位化探索

####變數概況

attributes(Insurance)
str(Insurance)
summary(Insurance)

####變數詳情

library(Hmisc)
describe(Insurance[, 1:3]) #想獲得變數詳情可用Hmisc套件中的describe()函數
describe(Insurance[, 4:5])

library(fBasics) #fBasics中有指標更豐富、面向更廣的basicStats()函數
basicStats(Insurance$Holders) #含95%信心水準的上下限、偏度峰度

####分佈指標

install.packages("timeDate")
library(timeDate)

skewness(Insurance[, 4:5]) #絕對值大於1為顯著偏倚，正值為右偏
kurtosis(Insurance[, 4:5]) #等於0為正態峰度（標準峰度）、>0為尖頂、<0為平頂

####稀疏性

library(Matrix) #Matrix套件處理高密度矩陣或稀疏矩陣
i <- sample(1:10, 10, replace = TRUE)
j <- sample(1:10, 10, replace = TRUE)
A <- sparseMatrix(i, j, x = 1) #對i行j列的元素設定為1，其他位置元素為空的稀疏矩陣A
loca <- which(A == 1, arr.ind = TRUE)
plot(loca, pch = 22)

####遺漏值

library(mice)
for(i in 1:10) {
        row <- sample(1:64, 1)
        col <- sample(1:5, 1)
        Insurance[row, col] = NA
}
md.pattern(Insurance) #檢查遺漏值分佈狀況

####相關性

cor(Insurance$Holders, Insurance$Claims)
library(rattle)
data(weather)
head(weather[, 12:21])
var <- c(12:21)
cor_matrix <- cor(weather[var], use = "pairwise")
cor_matrix

library(ellipse)
plotcorr(cor_matrix, col = rep(c("white", "black"), 5))
plotcorr(cor_matrix, diag = TRUE, type = "lower", col = rep(c("white", "black"), 5))


###視覺化探索

####長條圖

hist(Insurance$Claims, main = "Histogram of Freq of Insurance$Claims")
hist(Insurance$Claims, freq = FALSE, density = 20, 
     main = "Histogram of Density of Insurance$Claims") #各矩形面積和為1
lines(density(Insurance$Claims, na.rm = TRUE)) #線條下面積為1
str(hist(Insurance$Claims, breaks = 20, labels = TRUE, col = "black", #得相對應的輸出值
         border = "white", main = "Histogram of Insurance$Claims with 20 bars"))

####累積分佈圖

library(Hmisc)
Ecdf(Insurance$Claims, xlab = "Claims", main = "Cumulative Distribution of Claims") #繪製累積分佈圖

data_plot <- with(Insurance,
                  rbind(data.frame(var1 = Claims[Age == "<25"], var2 = "<25"),
                        data.frame(var1 = Claims[Age == "25-29"], var2 = "25-29"),
                        data.frame(var1 = Claims[Age == "30-35"], var2 = "30-35"),
                        data.frame(var1 = Claims[Age == ">35"], var2 = ">35")
                        ))
data_plot

Ecdf(data_plot$var1, lty = 2, group = data_plot$var2, label.curves = 1:4, xlab = "Claims",
     main = "Cumulative Distribution of Claims by Age")
Ecdf(Insurance$Claims, add = TRUE) #Claims變數及各Age組的累積分佈圖

####箱型圖

Claims_bp <- boxplot(Insurance$Claims, main = "Distribution of Claims")
Claims_bp$stats #取得箱型圖的5個界限值
points(x = 1, y = mean(Insurance$Claims, na.rm = TRUE), pch = 8)
Claims_points <- as.matrix(Insurance$Claims[which(Insurance$Claims > 102)], 6, 1)
Claims_text <- rbind(Claims_bp$stats, mean(Insurance$Claims, na.rm = TRUE), Claims_points)
for(i in 1:length(Claims_text)) text(x = 1.1, y = Claims_text[i, ], labels = Claims_text[i, ])

boxplot(var1 ~ var2, data = data_plot, horizontal = TRUE, 
        main = "Distribution of Claims by Age", xlab = "Claims", ylab = "Age")

with(Insurance, {
        boxplot(Holders ~ Age, boxwex = .25, at = 1:4 + .2, subset = Age == ">35")
        boxplot(Holders ~ Age, add = TRUE, boxwex = .25, at = 1:4 + .2, subset = Age == "30-35")
        boxplot(Holders ~ Age, add = TRUE, boxwex = .25, at = 1:4 + .2, subset = Age == "25-29")
        boxplot(Holders ~ Age, add = TRUE, boxwex = .25, at = 1:4 + .2, subset = Age == "<25")
})
boxplot(var1 ~ var2, data = data_plot, add = TRUE, boxwex = .25, at = 1:4 - .2, col = "lightgrey", main = "Distribution of Claims&Holders by Age", xlab = "Age", ylab = "Claims&Holders")
legend(x = "topleft", c("Claims", "Holders"), fill = c("lightgrey", "white"))


data_bp <- list(data_plot$var1[which(data_plot$var2 == "<25")],
                data_plot$var1[which(data_plot$var2 == "25-29")],
                data_plot$var1[which(data_plot$var2 == "30-35")],
                data_plot$var1[which(data_plot$var2 == ">35")]) #以list形式產生資料集
data_bp
bpplot(data_bp, name = c("<25", "25-20", "30-35", ">35"), xlab = "Age", ylab = "Claims") #繪製比例箱型圖

####橫條圖

#柱狀圖適合連續型變數，橫條圖適合離散型變數
Claims_Age <- with(Insurance,
                   c(sum(Claims[which(Age == "<25")], na.rm = TRUE),
                     sum(Claims[which(Age == "25-29")], na.rm = TRUE),
                     sum(Claims[which(Age == "30-35")], na.rm = TRUE),
                     sum(Claims[which(Age == ">35")], na.rm = TRUE)))
Claims_Age
barplot(Claims_Age, names.arg = c("<25", "25-29", "30-35", ">35"), density = rep(20, 4),
        main = "Distribution of Age by Claims", xlab = "Age", ylab = "Claims")

Holders_Age <- with(Insurance,
                    c(sum(Holders[which(Age == "<25")], na.rm = TRUE),
                      sum(Holders[which(Age == "25-29")], na.rm = TRUE),
                      sum(Holders[which(Age == "30-35")], na.rm = TRUE),
                      sum(Holders[which(Age == ">35")], na.rm = TRUE)))
Holders_Age
data_bar <- rbind(Claims_Age, Holders_Age)
data_bar

barplot(data_bar, names.arg = c("<25", "25-30", "30-35", ">35"), beside = TRUE, 
        main = "Age Distribution by Claims and Holders", xlab = "Age", ylab = "Claims&Holders",
        col = c("black", "darkgrey"))
legend(x = "topleft", rownames(data_bar), fill = c("black", "darkgrey"))

barplot(data_bar, names.arg = c("<25", "25-30", "30-35", ">35"), 
        main = "Age Distribution by Claims and Holders", xlab = "Age", ylab = "Claims&Holders",
        col = c("black", "darkgrey")) #beside = FALSE
legend(x = "topleft", rownames(data_bar), fill = c("black", "darkgrey"))

####點（陣）圖

dotchart(data_bar, xlab = "Claims&Holders", pch = 1:2, main = "Age Distribution by Claims and Holders")
legend(x = 14000, y = 15, "<25", bty = "n")
legend(x = 14000, y = 11, "25-29", bty = "n")
legend(x = 14000, y = 7, "30-35", bty = "n")
legend(x = 14000, y = 3, ">35", bty = "n")

####圓餅圖

pie(Claims_Age, labels = c("<25", "25-29", "30-35", ">35"),
    main = "Pie Chart of Claims by Age", col = c("white", "lightgrey", "darkgrey", "black"))

percent <- round(Claims_Age/sum(Claims_Age) * 100)
label <- paste(paste(c("<25", "25-29", "30-35"), ":"), percent, "%", sep = "")
pie(Claims_Age, labels = label, main = "Pie Chart of Claims by Age",
    col = c("white", "lightgrey", "darkgrey", "black"))

install.packages("plotrix")
library(plotrix)
pie3D(Claims_Age, labels = c("<25", "25-29", "30-35", ">35"), explode = .05, 
      main = "3D Pie Chart of Age by Claims", labelcex = .8, 
      col = c("white", "lightgrey", "darkgrey", "black"))