#實用R程式設計==================================================

#資料類別與基本運算==================================================

x <- 3
x
y <- 1:30
y
y[1]
y[18]
c(3, 7, 5)
x <- c(3, 7, 5)
x
y <- c("beautiful", "handsome")
y
z <- c(x, y)
z
class(z)
class(z[1])

-5:5
seq(from = -5, to = 5, by = 1)
seq(from = -5, to = 5, length = 11)
seq(from = -5, to = 5, by = 2)
seq(from = -5, to = 5, length = 6)
seq(from = 0, to = 1, by = .1)
seq(from = 0, to = 1, length = 11)
seq(from = -5, to = 5, length = 10)
seq(from = 0, to = 1, by = .3)
rep("beauty", 5)
rep("beauty", times = 5)
rep(3, times = 5)
rep(c(2, 3), times = 4)
rep(c(2, 3), length = 4)
rep(c(2, 3), each = 4)
rep(c(4, 7, 1, 8), c(3, 2, 5, 2))
rep(1:5, 1:5)
rep(1:5, rep(3, 5))
rep(1:5, c(2, 3)) #invalid

x <- c(2, 4, 6)
x
is.numeric(x)
is.integer(x)
is.character(x)
class(x)
class(c("beautiful", "handsome"))
class(c(2, "beauty"))

u1 <- c(2.2, 4.4, 6.6)
u2 <- as.numeric(c(2.2, 4.4, 6.6))
u3 <- as.vector(c(2.2, 4.4, 6.6), mode = "numeric")

v1 <- c(2, 4, 6)
v2 <- as.numeric(c(2, 4, 6))
v3 <- as.vector(c(2, 4, 6), mode = "numeric")

x1 <- as.integer(c(2, 4, 6))
x2 <- as.vector(c(2, 4, 6), mode = "integer")

y1 <- c("beautiful", "handsome")
y2 <- as.character(c("beautiful", "handsome"))
y3 <- as.vector(c("beautiful", "handsome"), mode = "character")

as.integer(c(-1.2, 0, 3.4))
as.character(c(2, 4, 6))
as.integer(c("beautiful", "cool", "handsome"))

numeric(3)
vector(mode = "numeric", length = 3)
integer(3)
vector(mode = "integer", length = 3)
character(3)
vector(mode = "character", length = 3)

numeric()
numeric(0)
length(numeric())

x <- c(1.2, 3.4, 5.6)
y <- as.numeric(c(1.2, 3.4, 5.6))
z <- numeric()
z[1:3] <- c(1.2, 3.4, 5.6)
z
w <- 1.2
w[2] <- 3.4
w[3] <- 5.6
w

aa <- NA
aa
class(aa)
length(aa)

bb <- NULL
bb
class(bb)
length(bb)

0 / 0
1 / 0
-1 / 0

x <- c(1, NA, 3)
x
x[2]
x[2] + x[2]
x[1] + x[2]
a <- Inf
b <- Inf
a + b
a * b
a / b

3.6 + 1.25
x <- 3.6
y <- 1.25
x + y
x - y
x * y
x / y
x <- 1:10
y <- 1:10 * 2
x
y
x + y
y - x
(y - x) / 2

a <- 5
b <- 2
a %/% b
a %% b

x <- 2
y <- 0:4
x
y
x^y
y^x

log(10)
log10(100)
log2(8)
exp(1.2)
x <- 0:7 * pi / 16
x

sx <- sin(x)
sx
x1 <- asin(sx)
x1 - x
cx <- cos(x)
cx
x2 <- acos(cx)
x2 - x
tx <- tan(x)
tx
x3 <- atan(tx)
x3 - x

factorial(4)
ceiling(2.55)
floor(2.55)
trunc(2.55)
round(3.555, digits = 2)
signif(3.555, digits = 2)
x <- c(1.2, 3.5, -4.7, 0)
min(x)
max(x)
range(x)
length(x)
sum(x)
mean(x)
sd(x)
var(x)
median(x)

x1 <- 1:6
x2 <- rep(4, 6)
x3 <- 6:1
x1
x2
x3
min(x1, x2, x3)
max(x1, x2, x3)
pmin(x1, x2, x3)
pmax(x1, x2, x3)

x <- c(1.22, 3.555, -4.75, 0)
ceiling(x)
floor(x)
trunc(x)
round(x, digits = 1)
signif(x, digits = 2)

x <- c(1.22, 3.55, -4.75, 0)
x[2]
x[c(1, 3)]
x[-2]
x[-c(2, 4)]
x[-length(x)]
x[c(1, 0, 2)]
x[c(-1, 0, -2)]
x[c(1, 0, -2)] #Error
x[c(1.2, 3.4)]
x[c(-1.2, -3.4)]
x <- c(1, 0, -1, 2)
y <- c(-2, 1, 0, -2)
x * y
2 * x
x * 2
c(rep(2, 4)) * x
x * (1:2)
x * c(rep(c(1, 2), 2))
sum(x * y)

look <- c("beautiful", "handsome", "cool")
substr(look, start = 1, stop = 4)
substring(look, first = 1, last = 4)
look.more <- c(look, "pretty")
look.more
paste(look, "people")
paste("They look like", look, "people", collapse = ". ")
letters[1:10]
LETTERS[1:10]
month.name[1:12]
month.abb[1:12]

A <- matrix(c(1.2, -6.5, -3.4, 4.3, 5.6, -2.1), nrow = 2, ncol = 3)
matrix(c(1.2, -3.4, 5.6, -6.5, 4.3, -2.1), nrow = 2, ncol = 3, byrow = TRUE)
class(A)
dim(A)
attributes(A)
nrow(A)
dim(A)[1]
ncol(A)
dim(A)[2]
u <- as.numeric(A)
u
v <- c(A)
v
dim(v)
length(v)
nrow(A) * ncol(A)
length(A)
matrix(0, nrow = 2, ncol = 3)
matrix(0, nrow = 3, ncol = 3)
diag(0, nrow = 3)
diag(3)
diag(2.5, nrow = 3)
diag(c(1, 2, 3), nrow = 3)
diag(c(1, 2, 3))
A <- matrix(1:9, nrow = 3, ncol = 3)
A
B <- matrix(1:9, nrow = 3, ncol = 3, byrow = TRUE)
B
A + B
A - B
A %*% B
A * B
A[2, 3]
A[1, 2]
A[4]
c(A)[4]
R1 <- A[1, ]
R1
R1[1, 2] #Error
R1 %*% B
RR1 <- A[1, , drop = FALSE]
RR1[1, 2]
RR1 %*% B
as.matrix(R1)
E <- A[c(1, 3), ]
E
class(E)
F <- A[c(1, 3), 2]
F
class(F)
C <- matrix(1:4, nrow = 2, ncol = 2)
D <- matrix(1:6, nrow = 2, ncol = 3)
cbind(C, D)
E <- matrix(1:4, nrow = 2, ncol = 2)
F <- matrix(1:6, nrow = 3, ncol = 2)
rbind(E, F)
t(A)
t(A) %*% A
diag(A)
sum(diag(A))
A <- matrix(c(1, 0, 0, 3, .5, 0, 2, 1, .25), nrow = 3, ncol = 3)
A
det(A)
Ainv <- solve(A)
Ainv %*% A
b <- c(2, 1, 3)
solve(A, b)

A <- array(1:24, dim = c(4, 3, 2))
A[, 2:3, ]
A[2:4, 2:3, ]
A[2:4, 2:3, 2]

camera <- list(c("Leica", "Pentax", "Olympus", "Nikon"), c(1.2, 3.4), c("red", "green", "blue"))
camera
camera <- list(brand = c("Leica", "Pentax", "Olympus", "Nikon"), real.number = c(1.2, 3.4), color = c("red", "green", "blue"))
camera
a1 <- camera[1]
a1
camera["brand"]
class(a1)
a2 <- camera[[1]]
a2
camera[["brand"]]
class(a2)
a3 <- camera$brand
a3
class(a3)

camera[1][1]
identical(camera[1][1], a1)
a2[c(1, 2)]
a3[2]

x1 <- c("father", "mother", "brother", "sister")
x2 <- c("Leica", "Pentax", "Olympus", "Nikon")
x3 <- c("gold", "red", "green", "blue")
x4 <- c(2, 1, 1, 2)

camera <- data.frame(member = x1, brand = x2, color = x3, amount = x4)
camera
class(camera)
names(camera)
colnames(camera)
rownames(camera)
camera$brand
camera[, "brand"]
identical(camera$brand, camera[, "brand"])
x5 <- c(8, 3, 2, 2)
camera$cost <- x5
camera
test <- camera
colnames(test)[c(4, 5)] <- c("number", "money")
test

camera[camera$brand == "Leica", ]
subset(camera, brand == "Leica")
camera[camera$cost > 2, ]
subset(camera, cost > 2)
identical(camera[camera$cost > 2, ], subset(camera, cost > 2))
transform(camera, log.cost = log(cost))
as.numeric(camera) #Error
as.matrix(camera)
class(as.matrix(camera))
is.character(as.matrix(camera))
A <- matrix(c(1.2, -6.5, -3.4, 4.3, 5.6, -2.1), nrow = 2, ncol = 3)
A
rownames(A)
colnames(A)
D <- as.data.frame(A)
D
names(D)
colnames(D)
rownames(D)
D$V1

x <- c("R", "G", "B", "R", "R", "B", "R", "G", "G")
x
class(x)
y <- factor(x)
y
class(y)
as.integer(y)
levels(y)
levels(y)[2]
nlevels(y)
levels(y)[as.integer(y)]

gl(5, 3)
gl(n = 5, k = 3)
class(gl(5, 3))
gl(5, 2, 13)
gl(n = 5, k = 2, length = 13)
is.factor(gl(5, 2, 13))

rm(list = ls())
ls()
mywd <- "/Users/thomas/Documents/R Files"
wd(mywd)
getwd()

x <- 1:10
y <- matrix(1:6, nrow = 2, ncol = 3)
x
y
ls()
dump(c("x", "y"), file = "dump.txt")
rm(x)
rm(y)
ls()

source(file = "dump.txt")
ls()
x
y

dput(y, file = "dput.txt")
newy <- dget("dput.txt")
newy

sink("test.txt")
x
y
sink("test.txt")

dataf <- iris[c(1, 2, 51, 52, 101, 102), c(1, 2, 5)]
dataf

dataf <- edit(dataf)
write.table(dataf, "dataf.txt")
list.files()
df <- read.table("dataf.txt")
df
sf0 <- read.table("dataf0.txt")
sf0
sf1 <- read.table("dataf1.txt", header = FALSE)
sf1
sf2 <- read.table("dataf2.txt", header = TRUE, row.names = NULL)
sf2
sf3 <- read.table("dataf3.txt", header = FALSE, row.names = NULL)
sf3

sf_csv2 <- read.table("dataf2.csv", sep = ",", row.names = NULL, header = TRUE)
sf_csv2
sf_csv2 <- read.csv("dataf2.csv")
sf_csv2
typhoon.data <- read.table("Typhoon-01.txt", header = TRUE)
x <- scan()
1.1
2.2
3.3

x

sv <- scan("scanvct.txt")
sv

s1 <- scan("scanlst.txt", list(Sepal.Length = 0, Sepal.Width = 0, Species = ""))

x <- c(1.2, -3.4, 5.6, -6, 0, 3)
cbind(x, sort(x), rank(x))
cbind(x, sort = sort(x), rank = rank(x))
data.frame(x, sort = sort(x), rank = rank(x))
y <- cbind(x, sort = sort(x), rank = rank(x))
z <- data.frame(x, sort = sort(x), rank = rank(x))
cat(x, "\n")
cat(y, "\n")
cat(z, "\n") #Error
print(x)
print(y)
print(z)

#邏輯運算與流程控制==================================================

as.logical(c(0, 1))
as.logical(c(-2.2, -1, 0, 1, 2.2))
as.logical(c("T", "TRUE", "True", "true"))
as.logical(c("F", "FALSE", "False", "false"))
as.logical("handsome") #NA

logical(3)
vector(mode = "logical", length = 3)

is.logical(3 < 5)
is.logical(c(TRUE, FALSE, FALSE, TRUE))
is.logical(c(-2.2, -1, 0, 1, 2.2))
is.logical("handsome")
a <- c(FALSE, TRUE, FALSE, TRUE)
b <- c(FALSE, TRUE, TRUE, FALSE)
!a
a & b
a | b
x <- 1:5
is.na(x)
y <- c(x, "NA")
is.na(y)
z <- c(x, NA)
is.na(z)
x <- c(-1.2, 0.5, 1.0, 1.3, 2.4, 5, 6.3)
any(1 < x & x < 5)
any(1 < x) & any(x < 5)
identical(any(1 < x & x < 5), any(1 < x) & any(x < 5))
all(1 < x & x < 5)
all(1 < x) & all(x < 5)

x <- c(1.2, -3.4, 5.7, -6, 0, 3)
which(x >= 1)
which((x >= 1) & (x <= 4))
which(x >= 6)

class(which(x >= 1))
x[which(x >= 1)]
length(which(x >= 1))
length(which(x >= 6))

3 < 5
3 > 5
class(3 < 5)
as.integer(3 < 5)
as.integer(3 > 5)
2.5 * (3 < 5)
2.5 * (3 > 5)

x <- c(1.2, -3.4, 5.7, -6, 0, 3)
x >= 0
as.integer(x >= 0)
sum(x >= 0)
table(x >= 0)
x[x >= 0]
y <- c(2.2, -4.4, 6.6, -8.8, 0, 3.3)
x < y
sum(x < y)
table(x < y)
x[x < y]

x <- 0.5 - .3
y <- .3 - .1
x
y
x == y
sprintf("%.20f", x)
sprintf("%.20f", y)

all.equal(x, y)
identical(all.equal(x, y), TRUE)
round(x, 10) == round(y, 10)

a <- 1:10 /16
a
sprintf("%.20f", a)
a <- 1:10 / 10
a
sprintf("%.20f", a)

x <- "handsome"
y <- 2
if(x == "beautiful") {
        y <- y + 3
}
y

x <- 2.5
y <- 4.7
if(x < y) {
        z <- x
} else {
        z <- y
}
z
min(x, y)

x <- -1.5
ifelse(x > 0, x, -x)
abs(x)

x <- c(1.2, -3.4, 5.7)
y <- c(0, 1, 2)

c1 <- 1
c2 <- 2

z <- ifelse(c1 < c2, x, y)
z

ifelse(rep(c1 < c2, length(x)), x, y)

if(c1 < c2) {
        z <- x
} else {
        z <- y
}

x <- 3
switch(x, 2 + 2, mean(1:10), 1:5)

x <- 6
switch(x, 2 + 2, mean(1:10), 1:5)

y <- "fruit"
switch(y, fruit = "banana", vegetable = "broccoli")

y <- "meat"
switch(y, fruit = "banana", "Neither")

centre <- function(x, type) {
        switch(type, 
               mean = mean(x),
               median = median(x),
               trimmed.mean = mean(x, trim = .1))
}

set.seed(1)
x <- rcauchy(10)
centre(x, "mean")
centre(x, "median")
centre(x, "trimmed.mean")

x <- .2
for(k in 2:5) {
        x[k] <- 4 * x[k-1] * (1 - x[k-1])
}
round(x, 4)

x <- 1:10
odd <- seq(from = 1, to = 9, by = 2)
even <- seq(from = 2, to = 10, by = 2)

odd.sum <- even.sum <- 0
for(i in odd) {
        odd.sum <- odd.sum + x[i]
}
for(i in even) {
        even.sum <- even.sum + x[i]
}
odd.sum
even.sum
odd.sum - even.sum

sum(x[odd]) - sum(x[even])

x <- c(1.2, 3.4, 2.1, 4.3, 3.2, 5.5, 6.7)
total <- x[1]
count <- 0
while(total <= 12){
        count <- count + 1
        total <- total + x[count + 1]
}
count
y <- cumsum(x)
sum(y < 12)

#函數與程式==================================================

array(1:12)
array(, c(3, 4))
array(1:12, c(3, 4))
array(dim = c(3, 4), data = 1:12)
args(array)
str(array)
x <- c(1, 8, 5, 2, 3, 1)
x
length(x)
diff(x)
sum(x)
prod(x)
max(x)
min(x)
which.max(x)
which.min(x)
range(x)
round(x * pi, 2)
cumsum(x)
cumprod(x)
unique(x)
mean(x)
median(x)
var(x)
sd(x)
summary(x)

x <- c(1.2, -3.4, 5.7, -6, 0, 3)
x
sort(x)
rank(x)
order(x)
x[order(x)]
identical(x[order(x)], sort(x))
order(x)[3]
which(rank(x) == 3)
sort(x)[3]
x[order(x)][3]
x[order(x)[3]]
sort(x, decreasing = TRUE)
rev(sort(x))
rev(rank(x))
order(x, decreasing = TRUE)
rev(order(x))
x <- c(2.5, 2.5, 2.5, 2.5, 2.3, 4.7, -2.2, 4.6, 4.6)
x
sort(x)
rank(x, ties.method = "average")
rank(x, ties.method = "first")
rank(x, ties.method = "random")
rank(x, ties.method = "max")
rank(x, ties.method = "min")
rank(x)

A <- matrix(1:12, nrow = 4, ncol = 3)
A
apply(A, MARGIN = 1, FUN = mean)
apply(A, MARGIN = 2, FUN = sum)
apply(A, MARGIN = 1, FUN = function(x) sd(x) / mean(x))

func <- function(x, c1, c2) c(mean(x[c1]), sum(x[c2]))
apply(A, MARGIN = 1, FUN = func, c1 = c(1, 2), c2 = c(2, 3))

A <- array(1:24, dim = c(4, 3, 2))
A
apply(A, MARGIN = 1, FUN = sum)
apply(A, MARGIN = c(1, 2), FUN = sum)
A <- matrix(1:12, nrow = 4, ncol = 3)
A
sweep(A, MARGIN = 2, STATS = 1:3, FUN = "+")
sweep(A, MARGIN = 2, STATS = 2:4, FUN = "*")
u <- apply(A, MARGIN = 2, FUN = mean)
u
func <- function(x) {
        u <- apply(x, MARGIN = 2, FUN = mean)
        v <- apply(x, MARGIN = 2, FUN = sd)
        w <- sweep(x, MARGIN = 2, STATS = u, FUN = "-")
        z <- sweep(w, MARGIN = 2, STATS = v, FUN = "/")
        return(z)
}
A <- matrix(1:12, nrow = 4, ncol = 3)
func(A)
scale(A)
attr(, "scaled:center")
attr(, "scaled:scale")
data(CO2)
str(CO2)

with(CO2, tapply(uptake, INDEX = Type, FUN = mean))
with(CO2, tapply(uptake, INDEX = Treatment, FUN = mean))
with(CO2, tapply(uptake, INDEX = list(Type, Treatment), FUN = mean))

apply(CO2[, 4:5], MARGIN = 2, FUN = mean)
sapply(CO2[, 4:5], FUN = mean)
lapply(CO2[, 4:5], FUN = mean)

x <- list(a = 1:3, beta = log(c(5, 6, 7)), logic = c(TRUE, FALSE, FALSE, TRUE))
x
sapply(x, FUN = median)
lapply(x, FUN = median)

f <- function(x) .01 * x ^ 3 * cos(x) - .2 * x ^ 2 * sin(x) + .05 * x - 1
f(-5)
f(0)
f(5)
f(c(-5, 0, 5))
class(f)
g <- f
g(c(-5, 0, 5))

f <- function(x) ifelse((-1 <= x) & (x <= 1), cos(x), 0)
f <- function(x) cos(x) * ((-1 <= x) & (x <= 1))
f(c(-2, 0, pi / 4, pi))

f <- function(x) 3 * x[1] - 4 * x[2] + x[1] * x [3]
f(c(0, 0, 0))
f(c(1, -1, 1))

sign(c(2, -2, 0))
sgn <- function(x) {
        if(x < 0) {
                value <- -1
        } else if (x == 0) {
                value <- 0
        } else {
                value <- 1
        }
        return(value)
}
sgn(2)
sgn(-2)
sgn(0)
sgn(c(2, -2, 0)) #warning

sgn.vec <- function(x) {
        n <- length(x)
        value <- integer(n)
        for(i in 1:n) {
                if(x[i] < 0) {
                        value[i] <- -1
                } else if(x[i] == 0) {
                        value[i] <- 0
                } else {
                        value[i] <- 1
                }
        }
        return(value)
}
sgn.vec(c(2, -2, 0))

f<- function(x) ifelse(x < 0, -1, ifelse(x > 0, 1, 0))
f(c(2, -2, 0))

f <- function(x) 2 * x + abs(x)
x <- c(1, 3, 2, 4)
f(x)
y <- matrix(1:6, nrow = 2, ncol = 3)
y
f(y)



#R繪圖==================================================

#探索資料圖形==================================================

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


#隨機樣本==================================================

#R程式範例==================================================

#常用進階技術==================================================


###R軟體資料分析基礎與應用==================================================

#04R軟體基礎==================================================

1 + 1
1 + 2 + 3
3 * 7 * 2
4 / 2
4 / 3
4 * 6 + 5
(4 * 6) + 5
4 * (6 + 5)

x <- 2
x
y = 5
y
3 <- z #Error

a <- b <- 7
a
b
assign("j", 4)
j
rm(j)
theVariable <- 17
theVariable
THEVARIABLE #Error

class(x)
is.numeric(x)
i <- 5L
i
is.integer(i)
is.numeric(i)
class(4L)
class(2.8)
4L * 2.8
class(4L * 2.8)
class(5L)
class(2L)
5L / 2L
class(5L / 2L)

x <- "data"
x
y <- factor("data")
y
nchar(x)
nchar("hello")
nchar(3)
nchar(452)
nchar(y) #Error

date1 <- as.Date("2012-06-28")
date1
class(date1)
as.numeric(date1)
ddate2 <- as.POSIXct("2012-06-28 17:42")
date2
class(date2)
as.numeric(date2)
class(date1)
class(as.numeric(date1))

TRUE * 5
FALSE * 5
k <- TRUE
class(k)
is.logical(k)
TRUE
T
class(T)
T <- 7
T
class(T)

2 == 3
2 != 3
2 < 3
2 <= 3
2 > 3
2 >= 3
"data" == "stats"
"data" < "stats"

x <- c(1:10)
x
x * 3
x + 2
x - 3
x / 4
x^2
sqrt(x)
1:10
10:1
-2:3
5:-7
x <- 1:10
y <- -5:4
x - y
x * y
x / y
x^y
length(x)
length(y)
length(x + y)
x + c(1, 2)
x + c(1, 2, 3)
x <= 5
x > y
x < y
x <- 10:1
y <- -4:5
any(x < y)
all(x < y)
q <- c("Hockey", "Football", "Baseball", "Curling", "Rugby", "Lacrosse", "Basketball", "Tennis", "Cricket", "Soccer")
nchar(q)
nchar(y)
x[1]
x[1:2]
x[c(1, 4)]

c(One = "a", Two = "y", Last = "r")
w <- 1:3
names(w) <- c("a", "b", "c")
w

q2 <- c(q, "Hockey", "Lacrosse", "Hockey", "Water Polo", "Hockey", "Lacrosse")
q2Factor <- as.factor(q2)
q2Factor
as.numeric(q2Factor)

factor(x = c("High School", "College", "Masters", "Doctorate"), levels = c("High School", "College", "Masters", "Doctorate"), ordered = TRUE)

mean(x)

?'+'
?'*'
?'=='
apropos("mea")

z <- c(1, 2, NA, 8, 3, NA, 3)
z
is.na(z)

zChar <- c("Hockey", NA, "Lacrosse")
zChar
is.na(zChar)

z <- c(1, NULL, 3)
z

d <- NULL
is.null(d)
is.null(7)

#05進階資料結構==================================================

x <- 10:1
y <- -4:5
q <- c("Hockey", "Football", "Baseball", "Curling", "Rugby", "Lacrosse", "Basketball", "Tennis", "Cricket", "Soccer")
theDF <- data.frame(x, y, q)
theDF
theDF <- data.frame(First = x, Second = y, Sport = q)
theDF
nrow(theDF)
ncol(theDF)
dim(theDF)

names(theDF)
names(theDF)[3]
rownames(theDF)
rownames(theDF) <- c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine", "Ten")
rownames(theDF)
rownames(theDF) <- NULL
rownames(theDF)
head(theDF)
head(theDF, n = 7)
tail(theDF)
class(theDF)
theDF$Sport
theDF["Sport"]
theDF[3, 2]
theDF[3, 2:3]
theDF[c(3, 5), 2]
theDF[c(3, 5), 2:3]
theDF[, 3]
theDF[, 2:3]
theDF[2, ]
theDF[2:4, ]
theDF[, c("First", "Sport")]
theDF[, "Sport"]
theDF["Sport"]
theDF[["Sport"]]

class(theDF[, c("First", "Sport")]) #data.frame
class(theDF[, "Sport"]) #factor vector
class(theDF[, "First"]) #integer vector
class(theDF["Sport"]) #data.frame
class(theDF[c("First", "Sport")]) #data.frame
class(theDF[["Sport"]]) #factor vector

theDF[, "Sport", drop = FALSE]
class(theDF[, "Sport", drop = FALSE])
class(theDF[, 3, drop = FALSE])

newFactor <- factor(c("Pennsylvania", "New York", "New Jersey", "New York", "Tennessee", "Massachusetts", "Pennsylvania", "New York"))
model.matrix(~newFactor - 1)

list1 <- list(1, 2, 3)
list1
list2 <-list(c(1, 2, 3))
list2
list3 <- list(c(1, 2, 3), 3:7)
list3
list4 <- list(theDF, 1:10)
list4
list5 <- list(theDF, 1:10, list3)
list5
names(list5)
names(list5) <- c("data.frame", "vector", "list")
names(list5)
list5
list6 <- list(TheDataFrame = theDF, TheVector = 1:10, TheList = list3)
names(list6)
list6

empty <- vector(mode = "list", length = 4)
empty

list5[[1]]
list5[["data.frame"]]
list5[[1]]$Sport
list5[["data.frame"]][, "Second"]

length(list5)
list5[[4]] <- 2
length(list5)
list5[["NewElement"]] <- 3:6
length(list5)
names(list5)
list5

A <- matrix(1:10, nrow = 5)
A
B <- matrix(21:30, nrow = 5)
B
C <- matrix(21:40, nrow = 2)
C
nrow(A)
ncol(A)
dim(A)
A + B
A * B
A == B
A %*% t(B)

colnames(A)
rownames(A)
colnames(A) <- c("Left", "Right")
rownames(A) <- c("1st", "2nd", "3rd", "4th", "5th")
colnames(B)
rownames(B)
colnames(B) <- c("First", "Second")
rownames(B) <- c("One", "Two", "Three", "Four", "Fifth")
colnames(C)
rownames(C)
colnames(C) <- LETTERS[1:10]
rownames(C) <- c("Top", "Bottom")

t(A)
A %*% C

theArray <- array(1:12, dim = c(2, 3, 2))
theArray

#06讀取各類資料==================================================

theUrl <- "http://www.jaredlander.com/data/Tomato%20First.csv"
tomato <- read.table(file = theUrl, header = TRUE, sep = ",")
head(tomato)


#07統計繪圖==================================================

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

#08建立R函數==================================================

#09流程控制==================================================

#10迴圈 - 迭代元素的傳統做法==================================================

#11群組操作==================================================

theMatrix <- matrix(1:9, nrow = 3)
apply(theMatrix, 1, sum)
apply(theMatrix, 2, sum)
rowSums(theMatrix)
colSums(theMatrix)
class(apply(theMatrix, 1, sum))
class(rowSums(theMatrix))
identical(as.numeric(apply(theMatrix, 1, sum)), rowSums(theMatrix))

theMatrix[2, 1] <- NA
apply(theMatrix, 1, sum)
apply(theMatrix, 1, sum, na.rm = TRUE)
rowSums(theMatrix)
rowSums(theMatrix, na.rm = TRUE)

theList <- list(A = matrix(1:9, 3), B = 1:5, C = matrix(1:4, 2), D = 2)
lapply(theList, sum)
sapply(theList, sum)
firstList <- list(A = matrix(1:16, 4), B = matrix(1:16, 2), C = 1:5)
secondList <- list(A = matrix(1:16, 4), B = matrix(1:16, 8), C = 15:1)
mapply(identical, firstList, secondList)
simpleFunc <- function(x, y){
        NROW(x) + NROW(y)
}
mapply(simpleFunc, firstList, secondList)

require(ggplot2)
data(diamonds)
head(diamonds)
aggregate(price ~ cut, diamonds, mean)
aggregate(price ~ cut + color, diamonds, mean)
aggregate(cbind(price, carat) ~ cut, diamonds, mean)

require(plyr)
head(baseball)
baseball$sf[baseball$year < 1954] <- 0
any(is.na(baseball$sf))
baseball$hbp[is.na(baseball$hbp)] <- 0
baseball <- baseball[baseball$sb >= 50, ]
baseball$OBP <- with(baseball, (h + bb + hbp) / (ab + bb + hbp + sf))
tail(baseball)
obp <- function(data){
        c(OBP = with(data, sum(h + bb + hbp) / sum(ab + bb + hbp + sf)))
}
careerOBP <- ddply(baseball, .variables = "id", .fun = obp)
careerOBP <- careerOBP[order(careerOBP$OBP, decreasing = TRUE), ]
head(careerOBP, 10)

theList <- list(A = matrix(1:9, 3), B = 1:5, C = matrix(1:4, 2), D = 2)
lapply(theList, sum)
llply(theList, sum)
identical(lapply(theList, sum), llply(theList, sum))

sapply(theList, sum)
laply(theList, sum)
identical(sapply(theList, sum), laply(theList, sum))

aggregate(price ~ cut, diamonds, each(mean, median))
system.time(dlply(baseball, "id", nrow))
iBaseball <- idata.frame(baseball)
system.time(dlply(iBaseball, "id", nrow))

require(data.table)
theDF <- data.frame(A = 1:10, B = letters[1:10], C = LETTERS[11:20], D = rep(c("One", "Two", "Three"), length.out = 10))
theDF
theDT <- data.table(A = 1:10, B = letters[1:10], C = LETTERS[11:20], D = rep(c("One", "Two", "Three"), length.out = 10))
theDT
class(theDF$B)
class(theDT$B)
diamondsDT <- data.table(diamonds)
diamondsDT
theDT[1:2, ]
theDT[theDT$A >= 7, ]
theDF[, c("A", "C")]
theDT[, list(A, C)]
theDT[, B]
theDT[, list(B)]
theDT[, "B", with = FALSE]
theDT[, c("A", "C"), with = FALSE]

tables()

setkey(theDT, D)
theDT
key(theDT)
tables()
theDT["One", ]
theDT[c("One", "Two"), ]
setkey(diamondsDT, cut, color)
diamondsDT[J("Ideal", "E"), ]

aggregate(price ~ cut, diamonds, mean)
diamondsDT[, mean(price), by = cut]
diamondsDT[, list(price = mean(price)), by = cut]
diamondsDT[, mean(price), by = list(cut, color)]
diamondsDT[, list(price = mean(price)), by = list(cut, color)]
diamondsDT[, list(price = mean(price), carat = mean(carat)), by = cut]
diamondsDT[, list(price = mean(price), carat = mean(carat), caratSum = sum(carat), by = cut)]
diamondsDT[, list(price = mean(price), carat = mean(carat)), by = list(cut, color)]

#12資料整理==================================================

#cbind和rbind資料合併

sport <- c("Hockey", "Baseball", "Football")
league <- c("NHL", "MLB", "NFL")
trophy <- c("Stanley Cup", "Commissioner's Trophy", "Vince Lombardi Trophy")
trophies1 <- cbind(sport, league, trophy)
trophies2 <- data.frame(sport = c("Basketball", "Golf"), league = c("NBA", "PGA"), trophy = c("Larry O'Brien Championship Trophy", "Wanamaker Trophy"), stringsAsFactors = FALSE)
trophies <- rbind(trophies1, trophies2)

cbind(Sport = sport, Association = league, Prize = trophy)


#資料連結

download.file(url = "http://jaredlander.com/data/US_Foreign_Aid.zip", destfile = "ForeignAid.zio")
unzip("ForeignAid.zip", exdir = "data12")

require(stringr)
theFiles <- dir("data12/", pattern = "\\.csv") #取得檔案的列表
for(a in theFiles) {
        nameToUse <- str_sub(string = a, start = 12, end = 18)
        temp <- read.table(file = file.path("data12", a), header = TRUE, #用file.path()來指定資料夾和檔名
                           sep = ",", stringsAsFactors = FALSE)
        assign(x = nameToUse, value = temp)
}


##用merge合併兩個data.frame

Aid90s00s <- merge(x = Aid_90s, y = Aid_00s, by.x = c("Country.Name", "Program.Name"),
                   by.y = c("Country.Name", "Program.Name")) #merge()可指定不同名稱的變數進行連結，但速度比其他方法慢
head(Aid90s00s)


##用plyr join合併data.frame

require(plyr)
Aid90s00sJoin <- join(x = Aid_90s, y = Aid_00s, 
                      by = c("Country.Name", "Program.Name")) #join()速度比較快，但變數名稱要一樣
head(Aid90s00sJoin)

frameNames <- str_sub(string = theFiles, start = 12, end = 18)
frameList <- vector("list", length(frameNames)) #建立空的list
names(frameList) <- frameNames
for(a in frameNames) {
        frameList[[a]] <- eval(parse(text = a)) #parse解析字元，轉換為變數
}

head(frameList[[1]]) #data.frame都在list裡了

head(frameList[["Aid_00s"]])
head(frameList[[5]])
head(frameList[["Aid_60s"]])

allAid <- Reduce(function(...){ #用Reduce()函數把所有元素連結在一起，加快運行速度
        join(..., by = c("Country.Name", "Program.Name"))
}, frameList) #Reduce的原理是把1跟2先加起來，再把3加進來
dim(allAid)
require(useful)
corner(allAid, c = 15)
bottomleft(allAid, c = 15)


##data.table中的資料合併

require(data.table)
dt90 <- data.table(Aid_90s, key = c("Country.Name", "Program.Name"))
dt00 <- data.table(Aid_00s, key = c("Country.Name", "Program.Name"))
dt0090 <- dt90[dt00] #90在左，00在右，該指令是左連結，連結時需要用到關鍵詞，在建立data.table時已指定


#用reshape2套件置換行、列資料

##melt

head(Aid_00s)
require(reshape2)
melt00 <- melt(Aid_00s, id.vars = c("Country.Name", "Program.Name"), variable.name = "Year", value.name = "Dollars")
tail(melt00, 10)

require(scales)
melt00$Year <- as.numeric(str_sub(melt00$Year, start = 3, 6)) #把Year中的FY去掉，轉成numeric
meltAgg <- aggregate(Dollars ~ Program.Name + Year, data = melt00, sum, na.rm = TRUE)
meltAgg$Program.Name <- str_sub(meltAgg$Program.Name, start = 1, end = 10)
ggplot(meltAgg, aes(x = Year, y = Dollars)) + 
        geom_line(aes(group = Program.Name)) +
        facet_wrap(~ Program.Name) +
        scale_x_continuous(breaks = seq(from = 2000, to = 2009, by = 2)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 0)) +
        scale_y_continuous(labels = multiple_format(extra = dollar, multiple = "B"))


##dcast

cast00 <- dcast(melt00, Country.Name + Program.Name ~ Year, #~左邊放要保留的，右邊放新指定的
                value.var = "Dollars")
head(cast00)


#13字串處理==================================================

##用paste建立字串

paste("Hello", "Jared", "and others")
paste("Hello", "Jared", "and others", sep = "/")
paste(c("Hello", "Hey", "Howdy"), c("Jared", "Bob", "David")) #paste()也允許向量化運算
paste("Hello", c("Jared", "Bob", "David"))
paste("Hello", c("Jared", "Bob", "David"), c("Goodbye", "Seeya"))
vectorOfText <- c("Hello", "Everyone", "out there", ".")
paste(vectorOfText, collapse = " ") #collapse參數將文字折疊


##用sprintf建立含有變數的字串

person <- "Jared"
partySize <- "eight"
waitTime <- 25

paste("Hello ", person, ", your party of ", partySize, 
      " will be seated in ", waitTime, " minutes.", sep = " ") #很麻煩

sprintf("Hello %s, your party of %s will be seated in %s minutes", person, partySize, waitTime)
sprintf("Hello %s, your party of %s will be seated in %s minutes", 
        c("Jared", "Bob"), c("eight", 16, "four", 10), waitTime)


##抽取文字

require(XML)
theURL <- "http://www.loc.gov/rr/print/list/057_chron.html"
presidents <- readHTMLTable(theURL, which = 3, as.data.frame = TRUE, 
                            skip.rows = 1, header = TRUE, stringsAsFactors = FALSE)

head(presidents)
tail(presidents$YEAR)
presidents <- presidents[1:64, ]
require(stringr)
yearList <- str_split(string = presidents$YEAR, pattern = "-") #拆開字串
head(yearList)

yearMatrix <- data.frame(Reduce(rbind, yearList)) #合併成一個matrix
head(yearMatrix)
names(yearMatrix) <- c("Start", "Stop")
presidents <- cbind(presidents, yearMatrix)
presidents$Start <- as.numeric(as.character(presidents$Start))
presidents$Stop <- as.numeric(as.character(presidents$Stop))
head(presidents)
tail(presidents)

str_sub(string = presidents$PRESIDENT, start = 1, end = 3)
str_sub(string = presidents$PRESIDENT, start = 4, end = 8)
presidents[str_sub(string = presidents$Start, start = 4, end = 4) == 1, c("YEAR", "PRESIDENT", "Start", "Stop")]


##正規表示法

johnPos <- str_detect(string = presidents$PRESIDENT, pattern = "John") #不知道字會出現在什麼地方用str_detect()
presidents[johnPos, c("YEAR", "PRESIDENT", "Start", "Stop")]

badSearch <- str_detect(presidents$PRESIDENT, "john")
goodSearch <- str_detect(presidents$PRESIDENT, ignore.case("John")) #若要忽略大小寫，ignore.case
sum(badSearch)
sum(goodSearch)

con <- url("http://jaredlander.com/data/warTimes.rdata")
load(con) #用load()載入
close(con) #用close()關閉連結

head(warTimes, 10)

warTimes[str_detect(string = warTimes, pattern = "-")]
theTimes <- str_split(string = warTimes, pattern = "(ACAEA)|-", n = 2) #維基百科用的分隔符號為ACAEA
head(theTimes)

which(str_detect(string = warTimes, pattern = "-"))
theTimes[[147]]
theTimes[[150]]

theStart <- sapply(theTimes, FUN = function(x) x[1])
head(theStart)

theStart <- str_trim(theStart) #移除不必要的空格
head(theStart)

str_extract(string = theStart, pattern = "January") #抽出含有January字串的地方，用str_extract()
theStart[str_detect(string = theStart, pattern = "January")] #找出含有January的元素並整項回傳，用str_detect()

head(str_extract(string = theStart, "[0-9][0-9][0-9][0-9]"), 20) #找四個數字一起出現的情況
head(str_extract(string = theStart, "[0-9]{4}"), 20) #搜尋由4位數組成的數字
head(str_extract(string = theStart, "\\d{4}"), 20) #\\d代表任何整數
str_extract(string = theStart, "\\d{1,3}") #搜尋任何一個出現過1次、2次或3次的數字
head(str_extract(string = theStart, pattern = "^\\d{4}"), 30) #前端用^後端用$
head(str_extract(string = theStart, pattern = "\\d{4}$"), 30)
head(str_extract(string = theStart, pattern = "^\\d{4}$"), 30)
head(str_replace(string = theStart, pattern = "\\d", replacement = "x"), 30) #將第一個數字取代為x
head(str_replace_all(string = theStart, pattern = "\\d", replacement = "x"), 30) #將所有數字取代為x
head(str_replace_all(string = theStart, pattern = "nndf1, 4g", replacement = "x"), 30) #取代任何由1位數到4位數組成的數字為x

commands <- c("<a href = index.html>The Link is here</a>", "<b>This is bold text</b>")
str_replace(string = commands, pattern = "<.+?>(.+?)<.+>", replacement = "\\1")
#.代表搜尋任何文字
#+代表至少搜尋1次
#?代表非貪婪搜尋，採用貪婪搜尋若字串中多組符合匹配結果，可能會抓到額外內容
#\\1是回溯引用的符號，將搜尋的規則取代掉，1代表用搜尋的第一組文字做取代
#完整詳見?regex


#14機率分佈==================================================

#15基本統計==================================================

#16線性模型==================================================

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


#18模型診斷==================================================

##殘差==================================================

housing <- read.table("housing.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
names(housing) <- c("Neighborhood", "Class", "Units", "YearBuilt", "SqFt", "Income", "IncomePerSqFt", "Expense", "ExpensePerSqFt", "NetIncome", "Value", "ValuePerSqFt", "Boro")
head(housing)

house1 <- lm(ValuePerSqFt ~ Units + SqFt + Boro, data = housing)
summary(house1)

require(coefplot)
coefplot(house1) #迴歸係數圖

require(ggplot2)
head(fortify(house1))

h1 <- ggplot(aes(x = .fitted, y = .resid), data = house1) + geom_point() + geom_hline(yintercept = 0) + 
        geom_smooth(se = FALSE) + labs(x = "Fitted Values", y = "Residuals")
h1

plot(house1, which = 1)
plot(house1, which = 1, col = as.numeric(factor(house1$model$Boro)))
legend("topright", legend = levels(factor(house1$model$Boro)), pch = 1, col = as.numeric(factor(levels(factor(house1$model$Boro)))),
       text.col = as.numeric(factor(levels(factor(house1$model$Boro)))), title = "Boro")

plot(house1, which = 2) #Q-Q圖

ggplot(house1, aes(sample = .stdresid)) + stat_qq() + geom_abline() #ggplot畫Q-Q圖

ggplot(house1, aes(x = .resid)) + geom_histogram()


##模型比較==================================================

house2 <- lm(ValuePerSqFt ~ Units * SqFt + Boro, data = housing)
house3 <- lm(ValuePerSqFt ~ Units + SqFt * Boro + Class, data = housing)
house4 <- lm(ValuePerSqFt ~ Units + SqFt * Boro + SqFt * Class, data = housing)
house5 <- lm(ValuePerSqFt ~ Boro + Class, data = housing)

multiplot(house1, house2, house3, house4, house5, pointSize = 2)

anova(house1, house2, house3, house4, house5) #anova()可以用來比較不同模型
anova(house1, house2, house3, house4, house5)$RSS #RSS（殘差平方和）越低越好

#變數增加，RSS本來就會降低，要小心

#也可比較AIC、BIC，越低越好

AIC(house1, house2, house3, house4, house5)
BIC(house1, house2, house3, house4, house5)


housing$HighValue <- housing$ValuePerSqFt >= 150
high1 <- glm(HighValue ~ Units + SqFt + Boro, data = housing, family = binomial(link = "logit")) #二元羅吉斯迴歸模型
high2 <- glm(HighValue ~ Units * SqFt + Boro, data = housing, family = binomial(link = "logit"))
high3 <- glm(HighValue ~ Units + SqFt * Boro + Class, data = housing, family = binomial(link = "logit"))
high4 <- glm(HighValue ~ Units + SqFt * Boro + SqFt * Class, data = housing, family = binomial(link = "logit"))
high5 <- glm(HighValue ~ Boro + Class, data = housing, family = binomial(link = "logit"))

anova(high1, high2, high3, high4, high5)
AIC(high1, high2, high3, high4, high5)
BIC(high1, high2, high3, high4, high5)


##交叉驗證（Cross Validation）

require(boot)
houseG1 <- glm(ValuePerSqFt ~ Units + SqFt + Boro, data = housing, family = gaussian(link = "identity"))
identical(coef(house1), coef(houseG1))

houseCV1 <- cv.glm(housing, houseG1, K = 5) #5折CV
houseCV1$delta #檢視誤差（原始CV誤差、調整後CV誤差），調整是基於沒有使用LOOCV的補償

houseG2 <- glm(ValuePerSqFt ~ Units * SqFt + Boro, data = housing)
houseG3 <- glm(ValuePerSqFt ~ Units + SqFt * Boro + Class, data = housing)
houseG4 <- glm(ValuePerSqFt ~ Units + SqFt * Boro + SqFt * Class, data = housing)
houseG5 <- glm(ValuePerSqFt ~ Boro + Class, data = housing)

houseCV2 <- cv.glm(housing, houseG2, K = 5)
houseCV3 <- cv.glm(housing, houseG3, K = 5)
houseCV4 <- cv.glm(housing, houseG4, K = 5)
houseCV5 <- cv.glm(housing, houseG4, K = 5)

cvResults <- as.data.frame(rbind(houseCV1$delta, houseCV2$delta, houseCV3$delta, houseCV4$delta, houseCV5$delta))
names(cvResults) <- c("Error", "Adjusted.Error")
cvResults$Model <- sprintf("houseG%s", 1:5)
cvResults

cvANOVA <- anova(houseG1, houseG2, houseG3, houseG4, houseG5)
cvResults$ANOVA <- cvANOVA$'Resid. Dev'
cvResults$AIC <- AIC(houseG1, houseG2, houseG3, houseG4, houseG5)$AIC

require(reshape2)
cvMelt <- melt(cvResults, id.vars = "Model", variable.name = "Measure", value.name = "Value")
cvMelt
ggplot(cvMelt, aes(x = Model, y = Value)) + geom_line(aes(group = Measure, color = Measure)) + 
        facet_wrap( ~ Measure, scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = .5)) + 
        guides(color = FALSE)

cv.work <- function(fun, k = 5, data, cost = function(y, yhat) mean((y - yhat) ^ 2),
                    response = "y", ...) {
        folds <- data.frame(Fold = sample(rep(x = 1:k, length.out = nrow(data))),
                            Row = nrow(data))
        error <- 0
        for(f in 1:max(folds$Fold))
        {
                theRows <- folds$Row[folds$Fold == f]
                mod <- fun(data = data[-theRows, ], ...)
                pred <- predict(mod, data[theRows, ])
                error <- error + cost(data[theRows, response], pred) *
                        (length(theRows)/nrow(data))
        }
        return(error)
}

cv1 <- cv.work(fun = lm, k = 5, data = housing, response = "ValuePerSqFt", formula = ValuePerSqFt ~ Units + SqFt + Boro)
cv2 <- cv.work(fun = lm, k = 5, data = housing, response = "ValuePerSqFt", formula = ValuePerSqFt ~ Units * SqFt + Boro)
cv3 <- cv.work(fun = lm, k = 5, data = housing, response = "ValuePerSqFt", formula = ValuePerSqFt ~ Units + SqFt * Boro + Class)
cv4 <- cv.work(fun = lm, k = 5, data = housing, response = "ValuePerSqFt", formula = ValuePerSqFt ~ Units + SqFt * Boro + SqFt * Class)
cv5 <- cv.work(fun = lm, k = 5, data = housing, response = "ValuePerSqFt", formula = ValuePerSqFt ~ Boro + Class)
cvResults <- data.frame(Model = sprintf("house%s", 1:5), Erro = c(cv1, cv2, cv3, cv4, cv5))
cvResults


##自助抽樣法（Boostrap）

require(plry)
baseball <- baseball[baseball$year >= 1990, ]
head(baseball)

bat.avg <- function(data, indices = 1:NROW(data), hits = "h", at.bats = "ab") {
        sum(data[indices, hits], na.rm = TRUE) / 
                sum(data[indices, at.bats], na.rm = TRUE)
}

bat.avg(baseball)

avgBoot <- boot(data = baseball, statistic = bat.avg, R = 1200, stype = "i") #boostrap呼叫1200次
avgBoot #原資料的測量、估計值的偏差和標準誤差

names(avgBoot)

boot.ci(avgBoot, conf = .95, type = "norm") #信賴區間

ggplot() + geom_histogram(aes(x = avgBoot$t), fill = "grey", color = "grey") + geom_vline(xintercept = avgBoot$t0 + c(-1, 1) * 2 * sqrt(var(avgBoot$t)), linetype = 2)

#boot套件還可以對時間序列和被設限的資料做boostrap
#只有罕見的情況不能用，像是要測量有偏估計量的不確定性，像是lasso取得的估計量


##逐步向前變數選取（Stepwise Variable Selection）

#step函數可以用來對所有可能模型進行迭代
##scope參數：用來指定可以接受的最小和最大模型
##direction參數：要在模型增添變數（forward）或移除變數（backward），還是雙向迭代（both）

nullModel <- lm(ValuePerSqFt ~ 1, data = housing) #最小模型基本上就是直線平均
fullModel <- lm(ValuePerSqFt ~ Units + SqFt * Boro + Boro * Class, data = housing) #最大模型

houseStep <- step(nullModel, scope = list(lower = nullModel, upper = fullModel, data = housing))

houseStep #顯示被挑選的模型

#LASSO迴歸是更好的變數選取方式
#boostrap可以用來檢視模型不確定性


#19正規化和壓縮方法==================================================

acs <- read.table("http://jaredlander.com/data/acs_ny.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(acs)
testFrame <- data.frame(First = sample(1:10, 20, replace = TRUE),
                       Second = sample(1:20, 20, replace = TRUE),
                       Third = sample(1:10, 20, replace = TRUE),
                       Fourth = factor(rep(c("Alice", "Bob", "Charlie", "David"),
                                           5)),
                       Fifth = ordered(rep(c("Edward", "Frank", "Georgia", "Hank", "Isaac"), 4)),
                       Sixth = rep(c("a", "b"), 10), stringsAsFactors = FALSE)
head(testFrame)
head(model.matrix(First ~ Second + Fourth + Fifth, testFrame)) #model.matrix()

#Fifth是ordered factor，level間有大小關係
#Fourth的level少了一個

require(useful)
head(build.x(First ~ Second + Fourth + Fifth, testFrame, contrasts = FALSE)) #build.x()
head(build.x(First ~ Second + Fourth + Fifth, testFrame, 
             contrasts = c(Fourth = FALSE, Fifth = TRUE))) #只對Fourth使用所有level

acs$Income <- with(acs, FamilyIncome >= 150000)
head(acs)
acsX <- build.x(Income ~ NumBedrooms + NumChildren + NumPeople + NumRooms + NumUnits + NumVehicles + NumWorkers + OwnRent + YearBuilt + ElectricBill + FoodStamp + HeatingFuel + Insurance + Language - 1, 
                data = acs, contrasts = FALSE) #建立預測函數矩陣
class(acsX)
dim(acsX)
topleft(acsX, c = 6)
topright(acsX, c = 6)
acsY <- build.y(Income ~ NumBedrooms + NumChildren + NumPeople + NumRooms + 
                        NumUnits + NumVehicles + NumWorkers + OwnRent + YearBuilt + 
                        ElectricBill + FoodStamp + HeatingFuel + Insurance + Language - 1, data = acs) #建立反應變數
head(acsY)
tail(acsY)


require(glmnet)
set.seed(1863561)
acsCV1 <- cv.glmnet(x = acsX, y = acsY, 
                    family = "binomial", nfold = 5) #cv.glmnet()可以自動交叉驗證的值，預設alpha為1（LASSO）

acsCV1$lambda.min
acsCV1$lambda.1se
plot(acsCV1)

coef(acsCV1, s = "lambda.1se") #點代表沒被選中的變數，LASSO會把高度相關的變數排除掉
plot(acsCV1$glmnet.fit, xvar = "lambda")
abline(v = log(c(acsCV1$lambda.min, acsCV1$lambda.1se)), lty = 2)


set.seed(71623)
acsCV2 <- cv.glmnet(x = acsX, y = acsY, 
                    family = "binomial", nfold = 5, alpha = 0) #alpha = 0建立脊迴歸模型
acsCV2$lambda.min
acsCV2$lambda.1se
coef(acsCV2, s = "lambda.1se") #每個變數都會被保留，只是會被壓縮接近0
plot(acsCV2)
plot(acsCV2$glmnet.fit, xvar = "lambda")
abline(v = log(c(acsCV2$lambda.min, acsCV2$lambda.1se)), lty = 2)

set.seed(2834673)
theFolds <- sample(rep(x = 1:5, length.out = nrow(acsX))) #建立層別，要觀測值每次執行都落在同一層
alphas <- seq(from = .5, to = 1, by = .05) #尋找最佳的alpha值要加一層交叉驗證，只考慮>0.5的alpha，因為傾向LASSO好過Ridge Reg.
set.seed(5127151)
cl <- makeCluster(2) #啟動叢集

library(doParallel) #進行平行化運算
registerDoParallel(cl)
before <- Sys.time()
acsDouble <- foreach(i = 1:length(alphas), 
                     .errorhandling = "remove", #若發生錯誤，該迭代跳過
                     .inorder = FALSE, #整合的先後次序不重要
                     .multicombine = TRUE, #能夠同時接受好幾個引數
                     .export = c("acsX", "acsY", "alphas", "theFolds"), #將幾個變數通過.export載入foreach environment
                     .packages = "glmnet") %dopar% #每個worker都載入glmnet，%dopar%讓foreach以平行運算的方式執行
{
        print(alphas[i])
        cv.glmnet(x = acsX, y = acsY, family = "binomial", nfolds = 5,
                  foldid = theFolds, alpha = alphas[i])
}
after <- Sys.time()
stopCluster(cl) #確保完成後將叢集終止
after - before
sapply(acsDouble, class) #檢測該list是一個有11個cv.glmnet的物件列表

extractGlmnetInfo <- function(object) {
        lambdaMin <- object$lambda.min 
        lambda1se <- object$lambda.1se #找出備選中的lambda
        
        whichMin <- which(object$lambda == lambdaMin) #找出lambda落在路徑何處
        which1se <- which(object$lambda == lambda1se)
        data.frame(lambda.min = lambdaMin, error.min = object$cvm[whichMin],
                   lambda.1se = lambda1se, error1se = object$cvm[which1se]) #建立data.frame含有備選中的lambda和相關錯誤訊息
}

alphaInfo <- Reduce(rbind, lapply(acsDouble, extractGlmnetInfo)) #整合到一個data.frame
alphaInfo2 <- plyr::ldply(acsDouble, extractGlmnetInfo) #也可以用ldply
identical(alphaInfo, alphaInfo2)

alphaInfo$Alpha <- alphas
alphaInfo

require(reshape2)
require(stringr)
alphaMelt <- melt(alphaInfo, id.vars = "Alpha", value.name = "Value", variable.name = "Measure")
alphaMelt$Type <- str_extract(string = alphaMelt$Measure, pattern = "(min)|(1se)")
alphaMelt$Measure <- str_replace(string = alphaMelt$Measure, pattern = "//,(min|1se)", replacement = "")
alphaCast <- dcast(alphaMelt, Alpha + Type ~ Measure, value.var = "Value")
together <- function (x) {
        for(i in 1:dim(x)[1]) {
                if(i %% 2 == 1) {
                        x[i, 3] = x[i, 4]
                }
        }
        print(x)
}
alphaCast <- together(alphaCast)
alphaCast <- alphaCast[, -4]
names(alphaCast)[3] <- "error"

ggplot(alphaCast, aes(x = Alpha, y = error)) + geom_line(aes(group = Type)) + facet_wrap(~Type, scales = "free_y", ncol = 1) + geom_point(aes(size = lambda))


#20非線性模型==================================================

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


##決策樹==================================================

require(rpart)
creditTree <- rpart(Credit ~ CreditAmount + Age + CreditHistory + Employment, data = credit)
creditTree
require(rpart.plot)
rpart.plot(creditTree, extra = 4)

#容易因為overfitting而導致很高的變異，模型會很不穩定，資料略有改變就會對模型造成很大影響


##隨機森林（Random Forest）==================================================

require(useful)
require(randomForest)
creditFormula <- Credit ~ CreditHistory + Purpose + Employment + Duration + Age + CreditAmount
creditX <- build.x(creditFormula, data = credit)
creditY <- build.y(creditFormula, data = credit)
creditForest <- randomForest(x = creditX, y = creditY)
creditForest


#21時間序列與自相關性==================================================

##自迴歸移動平均模型（Autoregressive Moving Average）==================================================

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


##VAR向量自我迴歸==================================================

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


##GARCH==================================================

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


#22資料分群==================================================

#K-means分群法==================================================

wine <- read.table("wine.csv", header = TRUE, sep = ",")
head(wine)

wineTrain <- wine[, which(names(wine) != "Cutivar")] #cultivar和組別非常相似，先排除掉

#所有資料必須是numeric！！！

set.seed(278613)
wineK3 <- kmeans(x = wineTrain, centers = 3)
wineK3

require(useful)
plot(wineK3, data = wineTrain) #將資料投影到二維空間

plot(wineK3, data = wine, class = "Cultivar")

#K-means以隨機條件作開始的，建議用不同的隨機初始條件多進行幾次，用nstart參數

set.seed(278613)
wineK3N25 <- kmeans(wineTrain, centers = 3, nstart = 25)
wineK3$size
wineK3N25$size

wineBest <- FitKMeans(wineTrain, max.clusters = 20, nstart = 25, seed = 278613) #尋找最佳分群數
wineBest
PlotHartigan(wineBest) #13群最佳（Hartigan法則，>10就可以增加）

table(wine$Cultivar, wineK3N25$cluster)
plot(table(wine$Cultivar, wineK3N25$cluster), main = "Confusion Matrix for wine Clustering", xlab = "Cultivar", ylab = "Cluster") #繪製混淆矩陣

#除了Hartigan法則，另一個可用的是Gap統計量
#Gap（差距）統計量，比較分群資料與自助抽樣法所抽出的樣本之間的群內相異度，測量觀測與預期之間的差距

require(cluster)
theGap <- clusGap(wineTrain, FUNcluster = pam, K.max = 20)
gapDF <- as.data.frame(theGap$Tab)
gapDF

require(ggplot2)
ggplot(gapDF, aes(x = 1:nrow(gapDF))) + #logW曲線
        geom_line(aes(y = logW), color = "blue") +
        geom_point(aes(y = logW), color = "blue") + 
        geom_line(aes(y = E.logW), color = "green") +
        geom_line(aes(y = E.logW), color = "green") +
        labs(x = "Number of Clusters")

ggplot(gapDF, aes(x = 1:nrow(gapDF))) + 
        geom_line(aes(y = gap), color = "red") +
        geom_point(aes(y = gap), color = "red") +
        geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), color = "red") +
        labs(x = "Number of Clusters", y = "Gap")

#K-means不能用在類別資料，而且容易受離群值影響


#PAM分割環繞物件法（Partitioning Around Medoids）==================================================

#K-medoids最常用的演算法為PAM

indicators <- c("BX.KLT.DINV.WD.GD.ZS", "NY.GDP.DEFL.KD.ZG", 
                "NY.GDP.MKTP.CD", "NY.GDP.MKTP.KD.ZG",
                "NY.GDP.PCAP.CD", "NY.GDP.PCAP.KD.ZG",
                "TG.VAL.TOTL.GD.ZS")
require(WDI)

wbInfo <- WDI(country = "all", indicator = indicators, start = 2011, end = 2011, extra = TRUE)
wbInfo <- wbInfo[wbInfo$region != "Aggregates", ] #移除Aggregates資訊
wbInfo <- wbInfo[which(rowSums(!is.na(wbInfo[, indicators])) > 0), ] #移除指標變數為NA的國家
wbInfo <- wbInfo[!is.na(wbInfo$iso2c), ]

rownames(wbInfo) <- wbInfo$iso2c
wbInfo$region <- factor(wbInfo$region) #因素化region, income, lending
wbInfo$income <- factor(wbInfo$income) #這樣他們的level有任何變化都能被考量在內
wbInfo$lending <- factor(wbInfo$lending)

keep.cols <- which(!names(wbInfo) %in% c("iso2c", "country", "year", "capital", "iso3c")) #找出要保留的直行
wbPam <- pam(x = wbInfo[, keep.cols], k = 12, keep.diss = TRUE, keep.data = TRUE) #分群
wbPam$medoids #顯示medoid觀測值


download.file(url = "http://jaredlander.com/data/worldmap.zip", destfile = "worldmap.zip")
unzip(zipfile = "worldmap.zip", exdir = "data") #解壓縮也可以用R完成

require(maptools)
world <- readShapeSpatial("data/world_country_admin_boundary_shapefile_with_fips_codes.shp")
head(world@data) #shapefile和WDI的國家編碼有差異

require(plyr)
world@data$FipsCntry <- as.character( #revalue(replace = c())
        revalue(world@data$FipsCntry, 
                replace = c(AU = "AT", AS = "AU", VM = "VN", BM = "MM", SP = "ES",
                            PO = "PT", IC = "IL", SF = "ZA", TU = "TR", IZ = "IQ", 
                            UK = "GB", EI = "IE", SU = "SD", MA = "MG", MO = "MA",
                            JA = "JP", SW = "SE", SN = "SG")))


#shapefile轉為data.frame才可以用ggplot2
world@data$id <- rownames(world@data)
require(ggplot2)
require(rgeos)
install.packages("gpclib") #一定要做這個
gpclibPermit()
world.df <- fortify(world, region = "id") #fortify()用來將shapefile轉為data.frame
head(world.df)

world.df <- join(world.df, world@data[, c("id", "CntryName", "FipsCntry")], by = "id")
head(world.df)

#開始結合分群資料和世界銀行資料

clusterMembership <- data.frame(FipsCntry = names(wbPam$clustering), Cluster = wbPam$clustering, stringsAsFactors = FALSE)
head(clusterMembership)

world.df <- join(world.df, clusterMembership, by = "FipsCntry")
world.df$Cluster <- as.character(world.df$Cluster)
world.df$Cluster <- factor(world.df$Cluster, levels = 1:12)

ggplot() + geom_polygon(data = world.df, aes(x = long, y = lat, group = group, fill = Cluster, color = Cluster)) + 
        labs(x = NULL, y = NULL) + coord_equal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_blank())

#PAM分群世界圖，灰色代表世銀沒有資訊，或者沒有結合好

wbPam$clusinfo #相異度資訊


##階層分群法==================================================

#不要求預先設定分群數
#可以用在類別資料

wineH <- hclust(d = dist(wineTrain))
plot(wineH)

keep.cols <- which(!names(wbInfo) %in% c("iso2c", "country", "year", "capital", "iso3c"))
wbDaisy <- daisy(x = wbInfo[, keep.cols]) #計算「相異度矩陣」（可用在類別資料）
wbH <- hclust(wbDaisy)
plot(wbH)

#計算距離的方法有single、complete、average、centroid，average常被認為最合適

wineH1 <- hclust(dist(wineTrain), method = "single")
wineH2 <- hclust(dist(wineTrain), method = "complete")
wineH3 <- hclust(dist(wineTrain), method = "average")
wineH4 <- hclust(dist(wineTrain), method = "centroid")

plot(wineH1, labels = FALSE, main = "Single")
plot(wineH2, labels  = FALSE, main = "Complete")
plot(wineH3, labels = FALSE, main = "average")
plot(wineH4, labels = FALSE, main = "centroid")

#可以用分群數來做修剪
plot(wineH)
rect.hclust(wineH, k = 3, border = "red")
rect.hclust(wineH, k = 13, border = "blue")

#可以用高度來做修剪
plot(wineH)
rect.hclust(wineH, h = 200, border = "red")
rect.hclust(wineH, h = 800, border = "blue")

#想要更有效率，可以考慮fastcluster套件
install.packages("fastcluster")


#R語言程序開發==================================================

library(datasets)
data(airquality)
cor(airquality)

x <- 1
print(x)
msg <- "hello"

x <- 5
x
print(x)

x <- 1:20
x
x <- c(.5, .6)
class(x)
x <- c(TRUE, FALSE)
class(x)
x <- c(T, F)
class(x)
x <- c("a", "b", "c")
class(x)
x <- 9:29
class(x)
x <- c(1 + 0i, 2 + 4i)
class(x)
x <- vector("numeric", length = 10)
class(x)

y <- c(1.7, "a")
class(y)
y <- c(TRUE, 2)
class(y)
y <- c("a", TRUE)
class(y)

x <- 0:6
class(x)
as.numeric(x)
as.logical(x)
as.character(x)

x <- c("a", "b", "c")
as.numeric(x)

m <- matrix(nrow = 2, ncol = 3)
m
dim(m)
attributes(m)
m <- matrix(1:6, nrow = 2, ncol = 3)
m
m <- 1:10
m
dim(m) <- c(2, 5)
m

x <- 1:3
y <- 10:12
cbind(x, y)
rbind(x, y)

x <- list(1, "a", TRUE, 1 + 4i)
x

x <- factor(c("yes", "yes", "no", "yes", "no"))
table(x)
unclass(x)
x <- factor(c("yes", "yes", "no", "yes", "no"), levels = c("yes", "no"))
x

x <- c(1, 2, NA, 10, 3)
is.na(x)
x <- c(1, 2, NaN, NA, 4)
is.na(x)
is.nan(x)

x <- data.frame(foo = 1:4, bar = c(TRUE, TRUE, FALSE, FALSE))
x
nrow(x)
ncol(x)

x <- 1:3
names(x)
names(x) <- c("foo", "bar", "norf")
x
names(x)

x <- list(a = 1, b = 2, c = 3)
x
m <- matrix(1:4, nrow = 2, ncol = 2)
dimnames(m) <- list(c("a", "b"), c("c", "d"))
m

x <- c("a", "b", "c", "c", "d", "a")
x[1]
x[2]
x[1:4]
x[x > "a"]
u <- x > "a"
u
x[u]

x <- matrix(1:6, 2, 3)
x[1, 2]
x[2, 1]
x[1, ]
x[, 2]
x <- matrix(1:6, 2, 3)
x[1, 2]
x[1, 2, drop = FALSE]
x <- matrix(1:6, 2, 3)
x[1, ]
x[1, ,drop = FALSE]
x <- list(foo = 1:4, bar = .6)
x[1]
x[[1]]
x$bar
x[["bar"]]
x["bar"]
x <- list(foo = 1:4, bar = .6, baz = "hello")
x[c(1, 3)]
x <- list(foo = 1:4, bar = .6, baz = "hello")
name <- "foo"
x[[name]]
x$name #NULL
x$foo
x <- list(a = list(10, 12, 14), b = c(3.14, 2.81))
x[[c(1, 3)]]
x[[1]][[3]]
x[[c(2, 1)]]

x <- list(aardvark = 1:5)
x$a
x[["a"]] #NULL
x[["a", exact = FALSE]]
x <- c(1, 2, NA, 4, NA, 5)
bad <- is.na(x)
x[!bad]

x <- c(1, 2, NA, 4, NA, 5)
y <- c("a", "b", NA, "d", NA, "f")
good <- complete.cases(x, y)
good
x[good]
y[good]

airquality[1:6, ]
good <- complete.cases(airquality)
airquality[good, ][1:6, ]

x <- 1:4
y <- 6:9
x + y
x > 2
x >= 2
y == 8
x * y
x / y
x <- matrix(1:4, 2, 2)
y <- matrix(rep(10, 4), 2, 2)
x * y
x / y
x %*% y

data <- rea.table("foo.txt")
initial <- read.table("datatable.txt", nrows = 100)
classes <- sapply(initial, class)
tabAll <- read.table("datatable.txt", colClasses = classes)
y <- data.frame(a = 1, b = "a")
dput(y)
dput(y, file = "y.R")
new.y <- dget("y.R")
new.y

x <- "foo"
y <- data.frame(a = 1, b = "a")
dump(c("x", "y"), file = "data.R")
rm(x, y)
source("data.R")
y
x

str(file)
con <- file("foo.txt", "r")
data <- read.csv(con)
close(con)
data <- read.csv("foo.txt")
con <- gzfile("words.gz")
x <- readLines(con, 10)
x

con <- url("http://www.jhsph.edu", "r")
x <- readLines(con)
head(x)

if(x > 3) {
        y <- 10
} else {
        y <- 0
}

y <- if(x > 3) {
        10
} else {
        0
}

for(i in 1:10) {
        print(i)
}

x <- c("a", "b", "c", "d")

for(i in 1:4) {
        print(x[i])
}

for(letter in x) {
        print(letter)
}

for(i in 1:4) print(x[i])

x <- matrix(1:6, 2, 3)
for(i in seq_len(nrow(x))) {
        for(j in seq_len(ncol(x))) {
                print(x[i, j])
        }
}

count <- 0
while(count < 10){
        print(count)
        count <- count + 1
}

z <- 5
while(z >= 3 && z <= 10) {
        print(z)
        coin <- rbinom(1, 1, 0.5)
        
        if(coin == 1) {
                z <- z + 1
        } else {
                z <- z - 1
        }
}

x0 <- 1
tol <- 1e-8
repeat {
        x1 <- computeEstimate() #not a real function
        if(abs(x1 - x0) < tol) {
                break
        } else {
                x0 <- x1
        }
}

for(i in 1:100) {
        if(i <= 20) {
                next
        }
}

mydata <- rnorm(100)
sd(mydata)
sd(x = mydata)
sd(x = mydata, na,rm = FALSE)
sd(na.rm = FALSE, x = mydata)
sd(na.rm = FALSE, mydata)

lm(data = mydata, y ~ x, model = FALSE, 1:100)
lm(y ~ x, mydata, 1:100, model = FALSE)

f <- function(a, b = 1, c = 2, d = NULL) {
        
}

f <- function(a, b) {
        a ^ 2
}
f(2)

f <- function(a, b) {
        print(a)
        print(b)
}
f(45) #no arg. b

x <- list(a = 1:5, b = rnorm(10))
lapply(x, mean)
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = nrorm(100, 5))
lapply(x, mean)

x <- 1:4
lapply(x, runif)

x <- 1:4
lapply(x, runif, min = 0, max = 10)

x <- list(a = matrix(1:4, 2, 2), b = matrix(1:6, 3, 2))
x

lapply(x, function(elt) elt[, 1])

x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(x, mean)

sapply(x, mean)
mean(x) #warning

str(tapply)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
f
tapply(x, f, mean)
tapply(x, f, mean, simplify = FALSE)
tapply(x, f, range)

str(split)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
split(x, f)
lapply(split(x, f), mean)

library(datasets)
head(airquality)
s <- split(airquality, airquality$Month)
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")]))
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm = TRUE))

x <- rnorm(10)
f1 <- gl(2, 5)
f2 <- gl(5, 2)
str(split(x, list(f1, f2)))
str(split(x, list(f1, f2), drop = TRUE))

str(mapply)
mapply(rep, 1:4, 4:1)
noise <- function(n, mean, sd) {
        rnorm(n, mean, sd)
}
noise(5, 1, 2)
noise(1:5, 1:5, 2)
mapply(noise, 1:5, 1:5, 2)

str(apply)

x <- matrix(rnorm(200), 20, 10)
apply(x, 2, mean)
apply(x, 1, sum)

x <- matrix(rnorm(200), 20, 10)
apply(x, 1, quantile, probs = c(.25, .75))
a <- array(rnorm(2 * 2 * 10), c(2, 2, 10))
apply(a, c(1, 2), mean)

x <- as.Date("1970-01-01")
x
unclass(x)
unclass(as.Date("1970-01-02"))

x <- Sys.time()
x
p <- as.POSIXlt(x)
p
names(unclass(p))
p$sec
x <- Sys.time()
x
unclass(x)
x$sec
p <- as.POSIXlt(x)
p$sec

Sys.setlocale("LC_TIME", "C") #IMPORTANT!!!
datestring <- c("January 10, 2012 10:40", "December 9, 2011 9:10")
x <- strptime(datestring, "%B %d, %Y %H:%M")
x
class(x)

x <- as.Date("2012-01-01")
y <- strptime("9 Jan 2011 11:34:21", "%d %b %Y %H:%M:%S")
x - y
x <- as.POSIXlt(x)
x - y

x <- as.Date("2012-03-01")
y <- as.Date("2012-02-28")
x - y

x <- as.POSIXct("2012-10-25 01:00:00")
y <- as.POSIXct("2012-10-25 06:00:00", tz = "GMT")
y - x

lm <- function(x) {
        x * x
}
lm

search()

make.power <- function(n) {
        pow <- function(x) {
                x ^ n
        }
        pow
}

cube <- make.power(3)
cube(3)
square <- make.power(2)
square(3)

ls(environment(cube))
get("n", environment(cube))

ls(environment(square))
get("n", environment(square))

y <- 10
f <- function(x) {
        y <- 2
        y ^ 2 + g(x)
}
g <- function(x) {
        x * y
}
f(3)

rm(y)
g <- function(x) {
        a <- 3
        x + a + y
}
g(2)
y <- 3
g(2)

make.NegLogLik <- function(data, fixed = c(FALSE, FALSE)) {
        params <- fixed
        function(p) {
                params[!fixed] <- p
                mu <- params[1]
                sigma <- params[2]
                a <- -.5 * length(data) * log(2 * pi * sigma ^ 2)
                b <- -.5 * sum((data - mu) ^ 2) / (sigma ^ 2)
                -(a + b)
        }
}

set.seed(1)
normals <- rnorm(100, 1, 2)
nLL <- make.NegLogLik(normals)
nLL
ls(environment(nLL))

optim(c(mu = 0, sigma = 1), nLL)$par
nLL <- make.NegLogLik(normals, c(FALSE, 2))
optimize(nLL, c(-1, 3))$minimum

nLL <- make.NegLogLik(normals, c(1, FALSE))
optimize(nLL, c(1e-6, 10))$minimum

nLL <- make.NegLogLik(normals, c(1, FALSE))
x <- seq(1.7, 1.9, len = 100)
y <- sapply(x, nLL)
plot(x, exp(-(y - min(y))), type = "l")

nLL <- make.NegLogLik(normals, c(FALSE, 2))
x <- seq(0.5, 1.5, len = 100)
y <- sapply(x, nLL)
plot(x, exp(-(y - min(y))), type = "l")

log(-1) #NaN

printmessage <- function(x) {
        if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x)
}

printmessage <- function(x) {
        if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x)
}
printmessage(1)
printmessage(NA)

printmessage2 <- function(x) {
        if(is.na(x))
                print("x is a missing value!")
        else if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x)
}

printmessages2 <- function(x) {
        if(is.na(x))
                print("x is a missing value!")
        else if(x > 0)
                print("x is greater than zero")
        else 
                print("x is less than or equal to zero")
        invisible(x)
}
x <- log(-1)
printmessage2(x)

mean(x)
traceback()
lm(y ~ x)
traceback()
debug(lm)
lm(y ~ x)

options(error = recover)
read.csv("nosuchfile")

dnorm(x, mean = 0, sd = 1, log = FALSE)
pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
rnorm(n, mean = 0, sd = 1)

x <- rnorm(10)
x
x <- rnorm(10, 20, 2)
x
summary(x)

set.seed(10)
x <- rbinom(100, 1, .5)
e <- rnorm(100, 0, 2)
y <- .5 + 2 * x + e
summary(y)

set.seed(1)
x <- rnorm(100)
log.mu <- .5 + .3 * x
y <- rpois(100, exp(log.mu))
summary(y)
plot(x, y)

set.seed(1)
sample(1:10, 4)
sample(1:10, 4)
sample(letters, 5)
sample(1:10)
sample(1:10)
sample(1:10, replace = TRUE)

system.time(readLines("http://www.jhsph.edu"))
hilbert <- function(n) {
        i <- 1:n
        1 / outer(i - 1, i, "+")
}
x <- hilbert(1000)
system.time(svd(x))

system.time({
        n <- 1000
        r <- numeric(n)
        for(i in 1:n) {
                x <- rnorm(n)
                r[i] <- mean(x)
        }
})

#獲取和整理數據==================================================


###利用R語言打通大數據的經脈==================================================

#02資料概覽==================================================

#03用R取得資料==================================================

#04探索性資料分析==================================================

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


#05資料前置處理==================================================

##資料集載入

install.packages("lattice")
install.packages("MASS")
install.packages("nnet")
library(lattice)
library(MASS)
library(nnet)

install.packages("mice")
library(mice)
data(nhanes2)
nrow(nhanes2)
ncol(nhanes2)
summary(nhanes2)
head(nhanes2)


##資料清理

pay <- c(11, 19, 14, 22, 14, 28, 13, 81, 12, 43, 11, 16, 31, 16, 23, 42, 22, 26, 17, 22, 13, 27, 180, 
         16, 43, 82, 14, 11, 51, 76, 28, 66, 29, 14, 14, 65, 37, 16, 37, 35, 39, 27, 14, 17, 13, 38, 
         28, 40, 85, 32, 25, 26, 16, 12, 54, 40, 18, 27, 16, 14, 33, 29, 77, 50, 19, 34)
par(mfrow = c(2, 2))
hist(pay)
dotchart(pay)
boxplot(pay, horizontal = TRUE)
qqnorm(pay); qqline(pay)


###遺漏值處理

sum(is.na(nhanes2)) #計算有多少樣本有遺漏值
sum(complete.cases(nhanes2)) #計算完整樣本的數量
md.pattern(nhanes2) #觀察遺漏值的情況

#若要以chl為依變數，以age、hyp、bmi為自變數建立線性迴歸模型，因nhanes2有遺漏值，
#不能直接建構，因此要用mice()函數，用下列方式建模型

imp <- mice(nhanes2, m = 4) #產生4組完整資料
#多重插補法（Multivariate Imputation）透過變數間的關係對缺失資料進行預測，
#利用蒙地卡羅方法產生多個完整資料集，再對這些資料集分別進行分析。

fit <- with(imp, lm(chl ~ age + hyp + bmi))
pooled <- pool(fit) #整理4個迴歸模型
summary(pooled)

####刪除法
####插補法

sub <- which(is.na(nhanes2[, 4] == TRUE))
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE[, 4] <- sample(dataTR[, 4], length(dataTE[, 4]), replace = TRUE) #隨機插補法
dataTE

sub <- which(is.na(nhanes2[, 4]) == TRUE)
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE[, 4] <- mean(dataTR[, 4]) #平均值法
dataTE

sub <- which(is.na(nhanes2[, 4]) == TRUE)
dataTR <- nhanes2[-sub, ]
dataTE <- nhanes2[sub, ]
dataTE
lm <- lm(chl ~ age, data = dataTR)
nhanes2[sub, 4] <- round(predict(lm, dataTE)) #迴歸模型插補法
head(nhanes2)

#熱平台插補法是在非遺漏值資料集中找一個與遺漏值所在樣本相似的樣本
accept <- nhanes2[which(apply(is.na(nhanes2), 1, sum) != 0), ] #有遺漏值的樣本
donate <- nhanes2[which(apply(is.na(nhanes2), 1, sum) == 0), ] #無遺漏值樣本
accept[1, ]
donate[1, ]
accept[2, ]
sa <- donate[which(donate[, 1] == accept[2, 1] & donate[, 3] == accept[2, 3] & donate[, 4] == accept[2, 4]), ]
#尋找跟accept第2個樣本相似的樣本
accept[2, 2] = sa[1, 2]
accept[2, ]
#熱平台插補法在變數數量很多時很難找到需要的樣本

#冷平台插補法是按照某些變數將資料分層，在層中對遺漏值使用平均值插補法
level1 <- nhanes2[which(nhanes2[, 3] == "yes"), ] #按照變數hyp分層，挑出hyp值為yes的樣本
level1

level1[4, 4] <- mean(level1[1:3, 4])
level1


###雜訊資料處理

install.packages("outliers")
library(outliers)
set.seed(1)
s1 <- .Random.seed #產生100個標準常態亂數
y <- rnorm(100)
outlier(y)
outlier(y, opposite = TRUE)
dotchart(y)
dim(y) <- c(20, 5) #將y的資料重新劃分成20*5的矩陣
outlier(y) #求每列的離群最遠值
outlier(y, opposite = TRUE) #離群最遠值的相反值

set.seed(1)
s1 <- .Random.seed
y <- rnorm(10)
outlier(y, logical = TRUE)
plot(y)

#離群值也可以用分群方法檢測

#雜訊檢查後，常用分箱、迴歸、電腦檢查和人工方法檢查結合等方法「光滑」資料

#分箱方法透過對資料排序，利用資料「近鄰」來光滑有序數據的局部光滑方法

set.seed(1)
s1 <- .Random.seed
x <- rnorm(12)
x <- sort(x)
dim(x) <- c(3, 4)
x[1, ] <- apply(x, 1, mean)[1]
x[2, ] <- apply(x, 1, mean)[2]
x[3, ] <- apply(x, 1, mean)[3]
x


###資料不一致的處理

x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE, FALSE, FALSE, TRUE))
x

probs <- c(1:3/4)
rt.value <- c(0, 0, 0) #設定傳回值為3個數字
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs)

probs <- c(1:4/4)
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs) #ERROR，傳回值與要求格式不一

rt.value <- c(0, 0, 0, 0) #設定傳回值為4個數字
vapply(x, quantile, FUN.VALUE = rt.value, probs = probs)
rt.value <- c(0, 0, 0, "") #設定傳回值為3個數字、1個字串


##資料整合

x <- cbind(sample(c(1:50), 10), sample(c(1:50), 10))
chisq.test(x) #對矩陣x進行卡方檢定，檢查兩列是否相關（Reject H0，變數不相關）

x <- cbind(rnorm(10), rnorm(10))
cor(x) #相關係數，相關係數大的話，其中一個可被視為容錯而刪除
cov(x) #協方差，正值表示兩個屬性趨向一起改變，負值則相反方向

x <- cbind(sample(c(1:10), 10, replace = TRUE), rnorm(10), rnorm(10))
head(x)
y <- unique(x[, 1])
sub <- rep(0, length(y))
for(i in 1:length(y)) {
        sub[i] <- which(x[, 1] == y[i])[1]
} #去掉重複觀測值
x <- x[sub, ]
head(x)


##資料轉換

#光滑：分箱、迴歸、分群等
#屬性建構：由指定的屬性建構出新屬性
#聚集：整理資料
#規範化：把資料按比例縮放，實質落入一個特定的小區間

set.seed(1)
s1 <- .Random.seed #設定亂數種子
a <- rnorm(5) #產生亂數
b <- scale(a) #對亂數標準化
b
a


#離散化：數值屬性的原始值用區間標籤或概念化標籤取代，將定量資料轉成定性資料

a <- c(.7063422, .7533599, .6675749, .6100253, .9341495, .6069284, .3462011)
n <- length(a)
la <- rep(0, n)
la[which(a > .5)] = 1
la #0~1的連續性資料，轉為0或1的離散型資料

#由額定資料產生概念分層：屬性，資料泛化可以視為資料合併

city <- c(6, 7, 6, 2, 2, 6, 2, 1, 5, 7, 2, 1, 1, 6, 1, 3, 8, 8, 1, 1)
province <- rep(0, 20)
province[which(city > 4)] <- 1
province


##資料精簡

#AIC準則（赤池資訊準則）：透過最佳模型選擇變數（AIC越小越好）
#LASSO：透過限制條件選擇變數
#分類樹、隨機森林：透過對分類效果的影響大小篩選屬性
#小波轉換、主成份分析：把原資料轉換或投影到較小的空間

x <- matrix(rnorm(100 * 20), 100, 20) #產生20行常態亂數
y <- rnorm(100) #產生常態亂數作為依變數
library(glmnet)
fit1 <- glmnet(x, y) #廣義線性回歸，自變數為分組的，預設為LASSO
b <- coef(fit1, s = .01) #s代表lambda值，lambda越小，約束越寬，篩選出的變數越多
b #b為變數係數，有值的代表有選入模型

predict(fit1, newx = x[1:10, ], s = c(.01, .005)) #lambda為0.01和0.005情況下的預測值

#數值精簡指用較小的資料表示形式取代原資料。如參數方法中使用模型估算資料，就可以只儲存模型參數代替實際資料，如回歸模型和對數線性模型
#非參數方法可以使用長條圖、分群、抽樣和資料立方體聚集等方法


#06連結分析==================================================

install.packages("arules")
library(arules)

#Apriori效率較低
#Eclat執行效率有所提升
#FP-Growth高效最佳化演算法

#Apriori參數預設值
##support = .1
##confidence = .8
##maxlen = 10
##minlen = 1
##target = "rules"/"frequent itemsets"（輸出連結規則/頻繁項集）

##appearance：對先決條件X（lhs）和連結結果Y（rhs）實際包含哪些項進行限制
##control：控制函數效能，對項集進行升冪（sort = 1）或降冪（sort = -1）排序，是否向使用者報告處理程序（verbose = TREU/FALSE）


library(arules)
data("Groceries")
summary(Groceries)
inspect(Groceries[1:10])  #inspect()

rules0 <- apriori(Groceries, parameter = list(support = .001, confidence = .5))
rules0 #5668筆規則
inspect(rules0[1:10])

#先設定得很低，再加強support或confidence來調整，設定值較高容易遺失有用資訊

rules1 <- apriori(Groceries, parameter = list(support = .005, confidence = .5))
rules1 #120筆

rules2 <- apriori(Groceries, parameter = list(support = .005, confidence = .60))
rules2 #22筆

rules3 <- apriori(Groceries, parameter = list(support = .005, confidence = .64))
rules3 #4筆
inspect(rules3)


rules.sorted_sup <- sort(rules0, by = "support") #透過support控制
inspect(rules.sorted_sup[1:5])

rules.sorted_con <- sort(rules0, by = "confidence") #透過confidence控制
inspect(rules.sorted_con[1:5])

rules.sorted_lift <- sort(rules0, by = "lift") #透過lift控制
inspect(rules.sorted_lift[1:5])

#用lift來篩選連結規則是最可靠的指標，結論也常常是最有用的


#只想知道芥末（mustard）的強連結商品，且只要2個商品連結
rules4 <- apriori(Groceries, parameter = list(maxlen = 2, supp = 0.001, conf = .1),
                  appearance = list(rhs = "mustard", default = "lhs"))
inspect(rules4)

itemsets_apr <- apriori(Groceries, parameter = list(supp = .001, target = "frequent itemsets"), 
                        control = list(sort = -1))
itemsets_apr
inspect(itemsets_apr[1:5]) #觀察銷量最高的商品（supp = 0.001, sort = -1）


#用eclat()來取得最適合進行bundle銷售的商品（eclat()無法產生關聯規則）
itemsets_ecl <- eclat(Groceries, parameter = list(minlen = 1, maxlen = 3, supp = .001, target = "frequent itemsets"),
                      control = list(sort = -1))
itemsets_ecl
inspect(itemsets_ecl[1:5]) #觀察前5個頻繁項目集（eclat(target = "frequent itemsets)）
#頻繁項集只和support設定有關，confidence值不影響


library(arulesViz)
rules5 <- apriori(Groceries, parameter = list(support = .002, confidence = .5))
plot(rules5) #顏色深淺為lift值高低

plot(rules5, measure = c("support", "lift"), shading = "confidence") #改成由confidence決定顏色

plot(rules5, interactive = TRUE) #設定互動參數interactive

plot(rules5, shading = "order", control = list(main = "Two key plot")) #shading = "order"點的顏色深淺代表連結規則中有多少樣商品

plot(rules5, method = "grouped") #lift是顏色深淺，support是尺寸大小

plot(rules5[1:50], method = "matrix", measure = "lift")
plot(rules5[1:50], method = "matrix3D", measure = "lift")
plot(rules5[1:50], method ="paracoord")


#07分群分析==================================================

countries <- read.csv("countries.csv")
dim(countries)
head(countries)
names(countries) <- c("country", "birth", "death")
var <- countries$country
var <- as.character(var)
head(var)
for(i in 1:68) row.names(countries)[i] <- var[i]
head(countries)

plot(countries$birth, countries$death)
C <- which(countries$country == "CHINA")
T <- which(countries$country == "TAIWAN")
H <- which(countries$country == "HONG KONG")
I <- which(countries$country == "INDIA")
U <- which(countries$country == "UNITED STATES")
J <- which(countries$country == "JAPAN")
M <- which.max(countries$birth) #which.max()找最大值
points(countries[c(C, T, H, I, U, J, M), -1], pch = 16)
legend(countries$birth[C], countries$death[C], "CHINA", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[T], countries$death[T], "TAIWAN", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[H], countries$death[H], "HONG-KONG", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[I], countries$death[I], "INDIA", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[U], countries$death[U], "UNITED STATES", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[J], countries$death[J], "JAPAN", bty = "n", xjust = .5, cex = .8)
legend(countries$birth[M], countries$death[M], countries$country[M], bty = "n", xjust = 1, cex = .8)


fit_km1 <- kmeans(countries[, -1], center = 3)
print(fit_km1)
fit_km1$centers
fit_km1$totss; fit_km1$tot.withinss; fit_km1$betweenss #希望組間差距遠大於組內差距
fit_km1$betweenss + fit_km1$tot.withinss

plot(countries[, -1], pch = (fit_km1$cluster-1))
points(fit_km1$centers, pch = 8)
legend(fit_km1$centers[1, 1], fit_km1$centers[1, 2], "Center_1", bty = "n", xjust = 1, yjust = 0, cex = .8)
legend(fit_km1$centers[2, 1], fit_km1$centers[2, 2], "Center_2", bty = "n", xjust = 0, yjust = 0, cex = .8)
legend(fit_km1$centers[3, 1], fit_km1$centers[3, 2], "Center_3", bty = "n", xjust = .5, cex = .8)


result <- rep(0, 67)
for(k in 1:67) {
        fit_km <- kmeans(countries[, -1], center = k)
        result[k] <- fit_km$betweenss/fit_km$totss
}
round(result, 2)

plot(1:67, result, type = "b", main = "Choosing the Optimal Number of Cluster",
     xlab = "number of cluster: 1 to 67", ylab = "betweenss/totss")
points(10, result[10], pch = 16)
legend(10, result[10], 
       paste("(10,", sprintf("%.1f%%", result[10] * 100),")", sep = ""), 
       bty = "n", xjust = .3, cex = .8)

fit_km2 <- kmeans(countries[, -1], center = 10)
cluster_CHINA <- fit_km2$cluster[which(countries$country == "CHINA")]
which(fit_km2$cluster == cluster_CHINA) #選出與中國同類別的國家和地區


library(cluster)
fit_pam <- pam(countries[, -1], 3)
print(fit_pam) #Medoids可以看中心點是誰

head(fit_pam$data)
fit_pam$call
fit_pam1 <- pam(countries[, -1], 3, keep.data = FALSE)
fit_pam1$data
fit_pam2 <- pam(countries[, -1], 3, cluster.only = TRUE)
print(fit_pam2)

which(fit_km$cluster != fit_pam$cluster)


##系譜分群（階層式）

fit_hc <- hclust(dist(countries[, -1]))
print(fit_hc)
plot(fit_hc)


group_k3 <- cutree(fit_hc, k = 3)
group_k3
table(group_k3)
group_h18 <- cutree(fit_hc, h = 18)
group_h18

table(group_h18)
sapply(unique(group_k3), function(g) countries$country[group_k3 == g])

plot(fit_hc)
rect.hclust(fit_hc, k = 4, border = "light grey")
rect.hclust(fit_hc, k = 3, border = "dark grey")
rect.hclust(fit_hc, k = 7, which =  c(2, 6), border = "dark grey")


#密度分群（DBSCAN）

library(fpc)
ds1 <- dbscan(countries[, -1], eps = 1, MinPts = 5)
ds2 <- dbscan(countries[, -1], eps = 4, MinPts = 5)
ds3 <- dbscan(countries[, -1], eps = 4, MinPts = 2)
ds4 <- dbscan(countries[, -1], eps = 8, MinPts = 2)
ds1; ds2; ds3; ds4

par(mfcol = c(2, 2))
plot(ds1, countries[, -1], main = "1: MinPts = 5 eps = 1")
plot(ds3, countries[, -1], main = "3: MinPts = 2 eps = 4")
plot(ds2, countries[, -1], main = "2: MinPts = 5 eps = 4")
plot(ds4, countries[, -1], main = "4: MinPts = 2 eps = 8")

d <- dist(countries[, -1])
max(d); min(d)

library(ggplot2)
interval <- cut_interval(d, 30)
table(interval)
which.max(table(interval))

for(i in 3:5) {
        for(j in 1:10) {
                ds <- dbscan(countries[, -1], eps = i, MinPts = j)
                print(ds)
        }
}


#期望值最大分群（EM）

library(mclust)
fit_EM <- Mclust(countries[, -1])
summary(fit_EM, parameters = TRUE)
countries_BIC <- mclustBIC(countries[, -1])
countries_BICsum <- summary(countries_BIC, data = countries[, -1])
countries_BICsum
countries_BIC
plot(countries_BIC, G = 1:7, col = "black")
names(countries_BICsum)
mclust2Dplot(countries[, -1], classification = countries_BICsum$classification, 
             parameters = countries_BICsum$parameters, col = "black")
countries_Dens <- densityMclust(countries[, -1])
plot(countries_Dens, countries[, -1], col = "grey", nlevels = 55)
plot(countries_Dens, type = "persp", col = grey(.8))


#08判別分析==================================================

##費希爾判別
#選出適當的投影軸，使類別內離差盡可能小，不同類別間離差盡可能大

#線性判別分析LDA：將樣本點投影到一維空間，實際常常更複雜
#二次判別分析QDA：投影至二維空間，使用許多二次曲面，常用的非線性判別函數


##貝氏判別：利用貝氏選出最大後驗機率的類別
#不怕雜訊和無關變數，但假設各特徵屬性之間是無關的（現實常常並非如此）


#距離判別：根據待判斷樣本與已知類別樣本之間的距離遠近做判別
#對類別域的交換或重疊較多的待分樣本集來說，此方法較為適合
#常用的有k最近鄰法（k-Nearest Neighbor; kNN）、有權數的k最近鄰（Weighted k-Nearest Neighbor）


##資料準備

library(kknn)
data(miete)
head(miete)
dim(miete)
summary(miete)

library(sampling)
n <- round(2 / 3 * nrow(miete) / 5) #訓練集佔資料總量2/3，計算每一等級中應取出的樣本數
n
sub_train <- strata(miete, stratanames = "nmkat", size = rep(n, 5), method = "srswor") #分層抽樣
head(sub_train)

data_train <- getdata(miete[, c(-1, -3, -12)], sub_train$ID_unit)
data_test <- getdata(miete[, -1, -3, -12], -sub_train$ID_unit)
dim(data_train)
dim(data_test)
head(data_test)


##線性判別分析==================================================

install.packages("MASS")
library(MASS)
fit_lda1 <- lda(nmkat ~., data_train) #執行線性判別
names(fit_lda1) #輸出項目名稱

fit_lda1$prior #先驗機率
fit_lda1$counts #各種類別的樣本數

fit_lda1$means #各變數在每一種類中的平均值

fit_lda1


fit_lda2 <- lda(data_train[, -12], data_train[, 12])
fit_lda2

plot(fit_lda1)

plot(fit_lda1, dimen = 1) #輸出1個判別式的圖形
plot(fit_lda1, dimen = 2)


pre_lda1 <- predict(fit_lda1, data_test) #預測
pre_lda1$class #各樣本的預測結果
pre_lda1$posterior #每一種類別的後驗機率

table(data_test$nmkat, pre_lda1$class) #混淆矩陣

error_lda1 <- sum(as.numeric(as.numeric(pre_lda1$class) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_lda1


##單純貝氏分類==================================================

install.packages("klaR")
library(klaR)
fit_Bayes1 <- NaiveBayes(nmkat ~., data_train)
names(fit_Bayes1)
fit_Bayes1$apriori
fit_Bayes1$tables #用於建立判別規則的所有變數在各種類別下的條件機率
fit_Bayes1$levels
fit_Bayes1$call
fit_Bayes1$usekernel
fit_Bayes$varnames

plot(fit_Bayes1, vars = "wfl", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))
plot(fit_Bayes1, vars = "mvdauer", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))
plot(fit_Bayes1, vars = "nmqm", n = 50, col = c(1, "darkgrey", 1, "darkgrey", 1))

fit_Bayes1 <- NaiveBayes(data_train[, -12], data_train[, 12])
pre_Bayes1 <- predict(fit_Bayes1, data_test)
pre_Bayes1

table(data_test$nmkat, pre_Bayes1$class)
error_Bayes1 <- sum(as.numeric(as.numeric(pre_Bayes1$class) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_Bayes1


##K最近鄰==================================================

install.packages("class")
library(class)
fit_pre_knn <- knn(data_tain[, -12], data_test[, -12], cl = data_train[, 12])
fit_pre_knn
table(data_test$nmkat, fit_pre_knn)
error_knn <- sum(as.numeric(as.numeric(fit_pre_knn) != as.numeric(data_test$nmkat)))/nrow(data_test)
error_knn

error_knn <- rep(0, 20)
for(i in 1:20){
        fit_pre_knn <- knn(data_train[, -12], data_test[, -12], cl = data_train[, 12], k = i)
        error_knn[i] <- sum(as.numeric(as.numeric(fit_pre_knn) != as.numeric(data_test$nmkat)))/nrow(data_test)
}
error_knn

plot(error_knn, type = "l", xlab = "K")


##有權數的K最近鄰演算法==================================================

install.packages("kknn")
library(kknn)
fit_pre_kknn <- kknn(nmkat ~., data_train, data_test[, -12], k = 5)
summary(fit_pre_kknn)
fit <- fitted(fit_pre_kknn)
fit

table(data_test$nmkat, fit)
error_kknn <- sum(as.numeric(as.numeric(fit) != as.numeric(data_test$nmkat))) / nrow(data_test)
error_kknn


#09決策樹==================================================

rpart(formula, data, weights, subset, na.action = na.rpart, method, model = FALSE, x = FALSE, y = TRUE, parms, control, cost, ...) #method參數用於選擇決策樹的類型，包含anova、poisson、class、exp四種，不進行設定的預設情況下，R會自行猜測，例如當y為因數型變數時，預設為class型，其中，anova對應的是回歸樹，class則是分類樹
rpart.control(minsplit = 20, minbucket = round(minsplit/3), cp = .01, maxcompete = 4, maxsurrogate = 5, usesurrogate = 2, xval = 10, surrogatestyle = 0, maxdepth = 30) #minsplit為每個節點最少需要幾個樣本，minbucket為每個葉節點最少需要幾個樣本，cp為複雜度參數，maxdepth控制樹的高度


library(mvpart)
data(car.test.frame)
head(car.test.frame)
str(car.test.frame)
summary(car.test.frame)

Group_Mileage <- matrix(0, 60, 1)
Group_Mileage[which(car.test.frame$Mileage >= 26)] = "A"
Group_Mileage[which(car.test.frame$Mileage <= 22)] = "C"
Group_Mileage[which(Group_Mileage == 0)] = "B"

car.test.frame$Group_Mileage = Group_Mileage
car.test.frame[1:10, c(4, 9)]

a <- round(1/4 * sum(car.test.frame$Group_Mileage == "A"))
b <- round(1/4 * sum(car.test.frame$Group_Mileage == "B"))
c <- round(1/4 * sum(car.test.frame$Group_Mileage == "C"))
a; b; c

library(sampling)
sub <- strata(car.test.frame, stratanames = "Group_Mileage", size = c(c, b, a), method = "srswor") #用strata()進行分層抽樣
sub

Train_Car <- car.test.frame[-sub$ID_unit, ]
Test_Car <- car.test.frame[sub$ID_unit, ]
nrow(Train_Car)
nrow(Test_Car)


library(rpart)
formula_Car_Reg <- Mileage ~ Price + Country + Reliability + Type + Weight + Disp. + HP
rp_Car_Reg <- rpart(formula_Car_Reg, Train_Car, method = "anova")
print(rp_Car_Reg)

printcp(rp_Car_Reg) #匯出回歸樹的CP表格

summary(rp_Car_Reg)

rp_Car_Reg1 <- rpart(formula_Car_Reg, Train_Car, method = "anova", minsplit = 10)
print(rp_Car_Reg1)
printcp(rp_Car_Reg1)

rp_Car_Reg2 <- rpart(formula_Car_Reg, Train_Car, method = "anova", cp = .1)
print(rp_Car_Reg2)
printcp(rp_Car_Reg2)

rp_Car_Reg3 <- prune.rpart(rp_Car_Reg, cp = .1)
print(rp_Car_Reg3)
printcp(rp_Car_Reg3)

rp_Car_Reg4 <- rpart(formula_Car_Reg, Train_Car, method = "anova", maxdepth = 1)
print(rp_Car_Reg4)
printcp(rp_Car_Reg4)


install.packages("rpart.plot")
library(rpart.plot)
rp_Car_Plot <- rpart(formula_Car_Reg, Train_Car, method = "anova", minsplit = 10)
print(rp_Car_Plot)
rpart.plot(rp_Car_Plot) #繪製決策樹

rpart.plot(rp_Car_Plot, type = 4) #更改決策樹類型
rpart.plot(rp_Car_Plot, type = 4, branch = 1) #branch
rpart.plot(rp_Car_Plot, type = 4, fallen.leaves = TRUE) #fallen.leaves

install.packages("maptree")
library(maptree)
draw.tree(rp_Car_Plot, col = rep(1, 7), nodeinfo = TRUE) #另一種決策樹製圖工具

plot(rp_Car_Plot, uniform = TRUE, main = "plot: Regression Tree") #plot()也可以畫決策樹
text(rp_Car_Plot, use.n = TRUE, all = TRUE)

post(rp_Car_Plot, file = "") #post()也可以畫決策樹


formula_Car_Cla <- Group_Mileage ~ Price + Country + Reliability + Type + Weight + Disp. + HP
rp_Car_Cla <- rpart(formula_Car_Cla, Train_Car, method = "class", minsplit = 5)
print(rp_Car_Cla)

rpart.plot(rp_Car_Cla, type = 4, fallen.leaves = TRUE)

pre_Car_Cla <- predict(rp_Car_Cla, Test_Car, type = "class")
pre_Car_Cla
p <- sum(as.numeric(pre_Car_Cla != Test_Car$Group_Mileage)) / nrow(Test_Car)
p
table(Test_Car$Group_Mileage, pre_Car_Cla)


install.packages("RWeka")
library(RWeka)
names(Train_Car) <- c("Price", "Country", "Reliability", "Mileage", "Type", "Weight", "Disp.", "HP", "Oil_Consumption")
formula <- Oil_Consumption ~ Price + Country + Reliability + Type + Weight + Disp. + HP
Train_Car$Oil_Consumption <- as.factor(Train_Car$Oil_Consumption) #C5.0的輸出變數資料型態必須是因素，若不是可用factor()函數來轉換
C45_0 <- J48(formula, Train_Car)
C45_0

summary(C45_0)

C45_1 <- J48(formula, Train_Car, control = Weka_control(M = 3))
C45_1
summary(C45_1)

plot(C45_1) #可能要先裝partykit套件


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


#12支援向量機==================================================

install.packages("e1071")
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


#13神經網路==================================================

install.packages("nnet")
library(nnet)

vector1 <- c("a", "b", "a", "c")
vector2 <- c(1, 2, 1, 3)
class.ind(vector1) #對模型中的y進行前置處理
class.ind(vector2)

url <- "http://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-white.csv"
wine <- read.csv(url, header = TRUE, sep = ";")

cha = 0
for(i in 1:4898) {
        if(wine[i, 12] > 6) cha[i] = "good"
        else if(wine[i, 12] > 5) cha[i] = "mid"
        else cha[i] = "bad"
}
wine[, 12] = factor(cha)
summary(wine$quality)

scale01 <- function(x) { #資料歸一化，NNet常見的前置處理，將所有資料轉化為[0, 1]之間的數，取消個維度數據間數量級的差別
        ncol = dim(x)[2] - 1
        nrow = dim(x)[1]
        new = matrix(0, nrow, ncol)
        for(i in 1:ncol) {
                max = max(x[, i])
                min = min(x[, i])
                for(j in 1:nrow) {
                        new[j, i] = (x[j, i] - min) / (max - min)
                }
        }
        new
}

names(wine)
set.seed(71)
samp <- sample(1:4898, 3000)
wine[samp, 1:11] <- scale01(wine[samp, ])
r = 1/max(abs(wine[samp, 1:11])) #參數rang的變化範圍
set.seed(101)
model1 <- nnet(quality ~., data = wine, subset = samp, #rang為初始隨機權重的範圍
               size = 4, rang = r, decay = 5e-4, maxit = 200)

x <- subset(wine, select = -quality)
y <- wine[, 12]
y <- class.ind(y)
set.seed(101)
model2 <- nnet(x, y, decay = 5e-4, maxit = 200, size = 4, rang = r) #decay為權重衰減精度

summary(model1)
summary(model2)

x <- wine[, 1:11]
pred <- predict(model1, x, type = "class") #用於預測的變數個數要跟建立模型時一樣！
set.seed(110)
pred[sample(1:4898, 8)]


xt <- wine[, 1:11]
pred <- predict(model2, xt) #第二種模型
dim(pred)
pred[sample(1:4898, 4), ] #呈現3個輸出結果的值，表示樣本為某種類別的機率

name <- c("bad", "good", "mid")
prednew <- max.col(pred) #max.col()確定每列中最大值的那行
prednewn <- name[prednew]
set.seed(201)
prednewn[sample(1:4898, 8)]


true <- max.col(y)
table(true, prednewn) #檢查預測精度


data(iris)

x <- iris[, -5]
y <- iris[, 5]

x <- subset(iris, select = -Species)
y <- class.ind(y)

model1 <- nnet(x, y, rang = 1/max(abs(x)), size = 4, maxit = 500, decay = 5e-4)
model2 <- nnet(x, y, rang = 1/max(abs(x)), size = 4, maxit = 500, decay = 5e-4)

#每次建模型使用的反覆運算初值都是不同的

model1$convergence #結果為0表示反覆運算會停止「並非」因為達到最大反覆運算次數
model2$convergence #所以最大反覆運算次數並非造成兩模型不同的主因

model1$value #最後值為模型擬合標準同模型權數衰減值的和，越小表示擬合效果越好
model2$value

name <- c("setosa", "versicolor", "virginica")
pred1 <- name[max.col(predict(model1, x))]
pred2 <- name[max.col(predict(model2, x))]
table(iris$Species, pred1) #展示模型精度
table(iris$Species, pred2)

#要得到最佳模型可以多嘗試不同的模型、測試每一節點數目下模型的誤判率


wine <- read.table("wine.txt", sep = ";")
names(wine) <- c("fixed", "volatile", "citric", "residual", "chlorides", 
                 "free", "total", "density", "PH", "sulphates", "alcohol", "quality")
set.seed(71)
wine <- wine[sample(1:4898, 3000), ]
nrow.wine <- dim(wine)[1]
scale01 <- function(x) { #歸一化程式
        ncol = dim(x)[2] - 1
        nrow = dim(x)[1]
        new <- matrix(0, nrow, ncol)
        for(i in 1:ncol) {
                max = max(x[, i])
                min = min(x[, i])
                for(j in 1:nrow) new[j, i] = (x[j, i] - min)/(max - min)
        }
        new
}

cha = 0
for(i in 1:nrow.wine) {
        if(wine[i, 12] > 6) cha[i] = "good"
        else if(wine[i, 12] > 5) cha[i] = "mid"
        else cha[i] = "bad"
}
wine[, 12] <- factor(cha)

set.seed(444)
samp <- sample(1:nrow.wine, nrow.wine * .7)
wine[samp, 1:11] <- scale01(wine[samp, ])
wine[-samp, 1:11] <- scale01(wine[-samp, ])
r = 1/max(abs(wine[samp, 1:11]))
n <- length(samp)

#嘗試不同隱藏層節點個數

err1 = 0
err2 = 0
for(i in 1:17) {
        set.seed(111)
        model <- nnet(quality ~., data = wine, maxit = 400, rang = r, size = i, subset = samp, decay = 5e-4)
        err1[i] <- sum(predict(model, wine[samp, 1:11], type = "class") != wine[samp, 12]) / n
        err2[i] <- sum(predict(model, wine[-samp, 1:11], type = "class") != wine[-samp, 12]) / (nrow.wine - n)
}

par(family = "Songti TC Light")
plot(1:17, err1, "l", col = 1, lty = 1, ylab = "模型誤判率", xlab = "隱藏層節點個數", 
     ylim = c(min(min(err1), min(err2)), max(max(err1), max(err2))))
lines(1:17, err2, col = 1, lty = 3)
points(1:17, err1, col = 1, pch = "+")
points(1:17, err2, col = 1, pch = "o")
legend(1, 0.53, "測試集誤判率", bty = "n", cex = 1.0) #後面有過度擬合的情形
legend(1, 0.35, "訓練集誤判率", bty = "n", cex = 1.0) #隱藏層節點3個就夠了

#測試不同訓練週期

errl1 <- 0
errl2 <- 0
for(i in 1:500) {
        set.seed(111)
        model <- nnet(quality ~., data = wine, maxit = i, rang = r, size = 3, subset = samp)
        errl1[i] <- sum(predict(model, wine[samp, 1:11], type = "class") != wine[samp, 12]) / n
        errl2[i] <- sum(predict(model, wine[-samp, 1:11], type = "class") != wine[-samp, 12]) / (nrow.wine - n)
}
plot(1:length(errl1), errl1, "l", ylab = "模型誤判率", xlab = "訓練週期", col = 1, 
     ylim = c(min(min(errl1), min(errl2)), max(max(errl1), max(errl2))))
lines(1:length(errl1), errl2, col = 1, lty = 3)
legend(250, .47, "測試集誤判率", bty = "n", cex = .8)
legend(250, .425, "訓練集誤判率", bty = "n", cex = .8) #模型的誤判率會趨於平穩

set.seed(111)
model <- nnet(quality ~., data = wine, maxit = 300, rang = r, size = 3, subset = samp)
x <- wine[-samp, 1:11]
pred <- predict(model, x, type = "class")
table(wine[-samp, 12], pred)


#14模型評估與選擇==================================================

install.packages("rattle", dependencies = TRUE)
install.packages("rattle", repos = "http://rattle.togaware.com", type = "source")
library(rattle) #要先安裝GTK+和GGobi
rattle()


###資料挖礦與大數據分析==================================================

#關聯規則==================================================

library(arules)
library(arulesViz)

data("IncomeESL")
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ] #刪除遺漏值
dim(IncomeESL)

Income <- as(IncomeESL, "transactions") #換成可以進行關聯分析的transactions物件，每個屬性值轉化為單一item
sort(itemFrequency(Income), decreasing = TRUE)
itemFrequencyPlot(Income, support = .2, cex.names = .8)

rules <- apriori(Income, parameter = list(support = .1, confidence = .6)) #支持度門檻.1，信賴度門檻.6
summary(rules)
plot(rules, measure = c("confidence", "lift"), shading = "support")
plot(rules, method = "grouped") #圓圈大小：支持度，顏色深淺：增益
rulesOwn <- subset(rules, subset = rhs %in% "householder status=own" & lift > 1) #想觀察特定族群
inspect(head(sort(rulesOwn, by = "support"), n = 5)) #inspect()


data("IncomeESL")
IncomeESL <- IncomeESL[complete.cases(IncomeESL), ]
IncomeESL[["income"]] <- factor((as.numeric(IncomeESL[["income"]]) > 6) + 1, 
                                levels = 1:2, labels = c("$40-", "$40+")) #對資料重新編碼

Income <- as(IncomeESL, "transactions")

rules <- apriori(Income, parameter = list(support = .2, confidence = .6))
rulesIncome <- subset(rules, subset = rhs %in% "income=$40+" & lift > 1)

inspect(sort(rulesIncome, by = "confidence"))


#決策樹分析==================================================

library(MASS)
library(rpart) #CART決策樹
data("Pima.tr")
summary(Pima.tr)
set.seed(1111)
cart <- rpart(type ~., Pima.tr, control = rpart.control(cp = 0)) #cp是複雜係數alpha
summary(cart)
par(xpd = TRUE)
plot(cart)
text(cart)

cart_prune <- prune(cart, cp = .03)
par(xpd = TRUE)
plot(cart_prune)
text(cart)

pre <- predict(cart, Pima.te, type = "class")
confusion_matrix <- table(Type = Pima.te$type, Predict = pre)
confusion_matrix
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
accuracy


install.packages("C50")
library(C50) #C5.0決策樹
library(MASS)
data("Pima.tr")
C50_tree <- C5.0(type ~., Pima.tr, control = C5.0Control(noGlobalPruning = T)) #F則會修剪樹
summary(C50_tree)
plot(C50_tree)

install.packages("CHAID", repos="http://R-Forge.R-project.org", type="source")
library(CHAID)
data("Pima.tr")
data("Pima.te")
Pima <- rbind(Pima.tr, Pima.te)
level_name <- {} #將資料離散化，CHAID只能做離散
for(i in 1:7){
        Pima[, i] = cut(Pima[, i], breaks = 3, ordered_result = T, include.lowest = T)
        level_name <- rbind(level_name, levels(Pima[, i]))
}
level_name <- data.frame(level_name)
rownames(level_name) <- colnames(Pima[1:7])
colnames(level_name) <- paste("L", 1:3, sep = "")
level_name

Pima.tr <- Pima[1:200, ]
Pima.te <- Pima[201:nrow(Pima), ]
set.seed(1111)
CHAID_tree <- chaid(type ~., Pima.tr)
CHAID_tree
plot(CHAID_tree)


#類神經網路==================================================

##倒傳遞類神經網路

library(MASS)
install.packages("RSNNS")
library(RSNNS)
data(Pima.tr)
set.seed(1111)
Pima.tr <- Pima.tr[sample(1:nrow(Pima.tr), length(1:nrow(Pima.tr))), ] #將資料順序重新排列
PimaValues <- Pima.tr[, 1:7]
PimaTargets <- decodeClassLabels(Pima.tr[, 8]) #目標屬性重新編碼，因為目標屬性為二分類變數
Pima.tr <- splitForTrainingAndTest(PimaValues, PimaTargets, ratio = .1) #資料切割
Pima.tr <- normTrainingAndTestSet(Pima.tr) #避免7個屬性不同尺度影響分析結果，進行標準化

model <- mlp(Pima.tr$inputsTrain, Pima.tr$targetsTrain, size = 14, learnFuncParams = .01, maxit = 100, inputsTest = Pima.tr$imputsTest, targetsTest = Pima.tr$targetsTest)
plotIterativeError(model)
weightMatrix(model)


p_table <- expand.grid(size = c(12, 13, 14, 15, 16), learning.rate = c(0.001, 0.01, 0.1))
for(i in 1:nrow(p_table)) {
        model <- mlp(Pima.tr$inputsTrain, Pima.tr$targetsTrain, size = p_table[i, 1], 
                     learnFuncParams = p_table[i, 2], 
                     maxit = 100, inputsTest = Pima.tr$inputsTest, targetsTest = Pima.tr$targetsTest)
        p_table$TestError[i] = model$IterativeTestError[100]
}
p_table

Pima.te[, 1:7] <- normalizeData(Pima.te[, 1:7])
predictions <- predict(model, Pima.te[, 1:7])
table <- confusionMatrix(Pima.te[, 8], predictions)
accuracy <- sum(diag(table)) / sum(table)
accuracy


##自我組織映射網路

library(MASS)
install.packages("kohonen")
library(kohonen)
data("Pima.tr")
Pima_class <- rbind(Pima.tr, Pima.te)[, 8]
Pima <- scale(rbind(Pima.tr, Pima.te)[, -8]) #將屬性標準化以避免不同尺度影響分群結果

set.seed(1111)
Pima_som <- som(data = Pima, grid = somgrid(4, 4, "hexagonal"), rlen = 1000, alpha = c(.05, .01))
plot(Pima_som, type = "changes")
plot(Pima_som, type = "dist.neighbours")
plot(Pima_som, type = "codes")
plot(Pima_som, type = "counts")


##適應性共振理論類神經網路

library(RSNNS)
data(snnsData)
patterns <- snnsData$art1_letters.pat
inputMaps <- matrixToActMapList(patterns, nrow = 7)
par(mfrow = c(3, 3))
for(i in 1:9) plotActMap(inputMaps[[i]])

model <- art1(patterns, dimX = 7, dimY = 5, learnFuncParams = c(.5, 0, 0), maxit = 100)
table(encodeClassLabels(model$fitted.values))


#群集分析==================================================

data(USArrests)
distance <- dist(USArrests, method = "euclidean") #還有manhattan, minkowski
hc <- hclust(distance, method = "single") #單一連結法用最小距離，還有complete, average, centroid, ward
plot(hc, hang = -1)

method = c("single", "complete", "average", "centroid", "ward.D")
overall_mean <- apply(USArrests, 2, mean)
tot_ss <- 0
for(i in 1:nrow(USArrests)) tot_ss <- tot_ss + sum((USArrests[i, ] - overall_mean)^2)
g_plot <- data.frame(group = 1:10, withinss = rep(0, 10), betweenss = rep(0, 10))
for(m in method) {
        hc <- hclust(dist(USArrests, method = "euclidean"), method = m)
        for(g in 2:nrow(g_plot)) {
                ##within SS
                withinss = 0
                group_mean = sapply(USArrests, tapply, cutree(hc, k = g), mean)
                for(i in 1:g) {
                        group = USArrests[cutree(hc, k = g) == i, ]
                        for(j in 1:nrow(group)) {
                                withinss = withinss + sum((group[j, ] - group_mean[i, ])^2)
                        }
                }
                g_plot$withinss[g] = withinss
                
                ##between SS
                g_plot$betweenss[g] = tot_ss - withinss
        }
        per = round(g_plot$betweenss/tot_ss, 3)
        plot(g_plot[2:10, "withinss"], type = "n", ylab = "Sum of squares", xlab = "Number of cluster", xlim = c(2, 10), ylim = c(2000, tot_ss), main = m)
        points(g_plot$group[2:10], g_plot[2:10, "withinss"], col = 1, lty = 1, pch = 1)
        lines(g_plot$group[2:10], g_plot[2:10, "withinss"], col = 1, lty = 1, pch = 1)
        points(g_plot$group[2:10], g_plot[2:10, "betweenss"], col = 2, lty = 2, pch = 2)
        lines(g_plot$group[2:10], g_plot[2:10, "betweenss"], col = 2, lty = 2, pch = 2)
        text(g_plot$group[2:10], g_plot[2:10, "betweenss"] - 10000, per[2:10], col = 4)
        legend("topleft", legend = c("wintin SS", "between SS"), col = 1:2, lty = 1:2, pch = 1:2)
}


distance <- dist(USArrests, method = "euclidean")
hc <- hclust(dist, method = "complete")
plot(hc, hang = -1)
rect.hclust(hc, 3)
cutree(hc, 3)


set.seed(1111)
g_plot <- data.frame(group = 1:10, withinss = rep(0, 10), betweenss = rep(0, 10))
for(g in 2:nrow(g_plot)) {
        for(g in 2:nrow(g_plot)) {
                km = kmeans(USArrests, centers = g)
                g_plot$withinss[g] = km$tot.withinss
                g_plot$betweenss[g] = km$betweenss
                
        }
        per = round(g_plot$betweenss/km$totss, 3)
        plot(g_plot[2:10, "withinss"], type = "n", ylab = "Sum of squares", xlab = "Number of cluster", xlim = c(2, 10), ylim = c(12000, 350000), main = "k-means")
        points(g_plot$group[2:10], g_plot[2:10, "withinss"], col = 1, lty = 1, pch = 1)
        lines(g_plot$group[2:10], g_plot[2:10, "withinss"], col = 1, lty = 1, pch = 1)
        points(g_plot$group[2:10], g_plot[2:10, "betweenss"], col = 2, lty = 2, pch = 2)
        lines(g_plot$group[2:10], g_plot[2:10, "betweenss"], col = 2, lty = 2, pch = 2)
        text(g_plot$group[2:10], g_plot[2:10, "betweens"] - 10000, per[2:10], col = 4)
        legend("topleft", legend = c("within SS", "between SS"), col = 1:2, lty = 1:2, pch = 1:2)
}

kmeans(USArrests, centers = 3)$cluster


#簡單貝氏分類法與貝氏網路==================================================

library(MASS)
library(RSNNS)
data("Pima.tr")
data("Pima.tr")
set.seed(1111)

Pima <- rbind(Pima.tr, Pima.te)
level_name <- {}

for(i in 1:7) { #貝氏僅支援類別變數，要將連續型屬性離散化
        Pima[, i] = cut(Pima[, i], breaks = 2, ordered_result = TRUE, include.lowest = TRUE)
        level_name <- rbind(level_name, levels(Pima[, i]))
}

level_name <- data.frame(level_name)
row.names(level_name) <- colnames(Pima)[1:7]
colnames(level_name) <- paste("L", 1:2, sep = "")
level_name
Pima.tr <- Pima[1:200, ]
Pima.te <- Pima[201:nrow(Pima), ]


install.packages("bnlearn")
library(bnlearn)

bn <- naive.bayes(Pima.tr, "type") #naive.bayes建立簡單貝氏分類法
plot(bn)
bn
pred <- predict(bn, Pima.te)
tab <- table(pred, Pima.te[, "type"])
tab
acc <- sum(diag(tab)) / sum(tab)
acc


tan <- tree.bayes(Pima.tr, "type") #tree.bayes建構貝氏網路
plot(tan)
tan
fitted <- bn.fit(tan, Pima.tr, method = "bayes")
pred <- predict(fitted, Pima.te)
tab <- table(pred, Pima.te[, "type"])
tab
acc <- sum(diag(tab)) / sum(tab)
acc


#多變量分析==================================================

library(car)
library(MASS)
install.packages("lmtest")

data(UScrime)
UScrime$So = factor(UScrime$So)
summary(UScrime)
scatterplotMatrix(UScrime)

reg <- lm(y~., data = UScrime)
summary(reg)

residualPlots(reg)
raintest(y~., data = UScrime) #檢定回歸模型線性假設（殘差檢定1）

qqnorm(residuals(reg))
qqline(residuals(reg))

(residuals(reg)) #檢定殘差是否服從常態性（殘差檢定2）

durbinWatsonTest(reg) #殘差是否服從獨立性（殘差檢定3）

bptest(reg) #檢定模式中是否具殘差均值性（殘差檢定4）

outlierTest(reg) #檢定模式中是否有離群值（殘差檢定5）

vif(reg) #檢定模式中的共線性


LG <- glm(type~., data = Pima.tr, family = binomial)
summary(LG)
pre_LG <- predict(LG, Pima.tr, type = "response")
acc = {}
for(t in seq(.1, .9, .1)){
        pre = ifelse(pre_LG > t, 1, 0)
        tab = table(Pima.tr$type, pre)
        acc = c(acc, sum(diag(tab))/nrow(Pima.tr))
}

t <- seq(.1, .9, .1)[which.max(acc)]
pre_LG <- predict(LG, Pima.te, type = "response")
pre <- ifelse(pre_LG >= t, 1, 0)
tab <- table(Pima.te$type, pre)
acc <- sum(diag(tab))/nrow(Pima.te)
acc

#時間資料分析==================================================

install.packages("forecast")
library(forecast)
install.packages("TSA")
library(TSA)
data(co2, package = "datasets")
tsdisplay(co2)
acf(co2)
pacf(co2)

train <- ts(co2[seq(1, length(co2)-12)], frequency = 12, start = c(1959, 1))
test <- ts(co2[seq(length(co2)-11, length(co2))], frequency = 12, start = c(1997, 1))

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


###中華R軟體協會

#中華R軟體協會 4-2 R軟體假設檢定==================================================

t.test.data <- rnorm(50, mean = 100, sd = 15)
t.test(t.test.data, mu = 95)

wilcox.test.data <- rnorm(50, mean = 100, sd = 15)
wilcox.test(wilcox.test.data, conf.int = TRUE, conf.level = .99)

prop.test(39, 215, .15)

shapiro.test(rnorm(100, mean = 5, sd = 3))
shapiro.test(runif(100, min = 2, max = 4))

data(Cars93, package = "MASS")
qqnorm(Cars93$Price, main = "Q-Q Plot: Price")
qqline(Cars93$Price)

qqnorm(log(Cars93$Price), main = "Q-Q Plot: log(Price")
qqline(log(Cars93$Price), col = 2)

#中華R軟體協會 4-3 R軟體變異數分析==================================================

drink.sales <- read.table("drink.csv", header = TRUE, sep = ",")
head(drink.sales)
drink.type <- gl(4, 5, label = c(letters[1:4]))
drink.type
drink <- data.frame(drink.type = drink.type, drink.sales)
head(drink)
class(drink)

drink.oneway <- oneway.test(drink$sales ~ drink$drink.type, var.equal = TRUE)
drink.oneway

drink.anova <- aov(drink$sales ~ drink$drink.type)
summary(drink.anova)

drink.lm <- lm(drink$sales ~ drink$drink.type)
anova(drink.lm)

#中華R軟體協會 4-4 R軟體線性回歸==================================================

head(cars)
dim(cars)
cars.lm <- lm(dist~speed, data = cars)
summary(cars.lm)


###應用R語言於資料分析==================================================

#簡介==================================================

#資料的讀取與寫入==================================================

#流程控制及自訂函數==================================================

#繪圖功能==================================================

#相關套件介紹==================================================

#監督式學習 6.1 決策樹==================================================

library(rpart) #CART
data(iris)
np <- ceiling(.1 * nrow(iris)) #抽10%當作測試資料
np
test.index <- sample(1:nrow(iris), np)
iris.testdata <- iris[test.index, ]
iris.traindata <- iris[-test.index, ]

iris.tree <- rpart(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, method = "class", data = iris.traindata)
iris.tree
summary(iris.tree)

plot(iris.tree)
text(iris.tree)

species.traindata <- iris$Species[-test.index]
train.predict <- factor(predict(iris.tree, iris.traindata, type = 'class'), levels = levels(species.traindata))
table.traindata <- table(species.traindata, train.predict)
table.traindata
correct.traindata <- sum(diag(table.traindata)) / sum(table.traindata) * 100
correct.traindata #訓練資料正確率

species.testdata <- iris$Species[test.index]
test.predict <- factor(predict(iris.tree, iris.testdata, type = 'class'), levels = levels(species.testdata))
table.testdata <- table(species.testdata, test.predict)
table.testdata
correct.testdata <- sum(diag(table.testdata)) / sum(table.testdata) * 100
correct.testdata


library(rpart)
data(iris)
np <- ceiling(.1 * nrow(iris))
np
test.index <- sample(1:nrow(iris), np)
iris.testdata <- iris[test.index, ]
iris.traindata <- iris[-test.index, ]

iris.tree <- rpart(Species ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, method = "class", data = iris.traindata, control = rpart.control(minsplit = 5, cp = .0001, maxdepth = 30))
species.traindata <- iris$Species[-test.index]
train.predict <- factor(predict(iris.tree, iris.traindata, type = "class"), levels = levels(species.traindata))

table.traindata <- table(species.traindata, train.predict)
table.traindata
correct.traindata <- sum(diag(table.traindata)) / sum(table.traindata) * 100
correct.traindata

species.testdata <- iris$Species[test.index]
test.predict <- factor(predict(iris.tree, iris.testdata, type = "class"), levels = levels(species.testdata))
table.testdata <- table(species.testdata, test.predict)
table.testdata
correct.testdata <- sum(diag(table.testdata)) / sum(table.testdata) * 100
correct.testdata


library(C50) #載入C50套件
data(iris)
np <- ceiling(.1 * nrow(iris))
np

test.index <- sample(1:nrow(iris), np)
iris.test <- iris[test.index, ]
iris.train <- iris[-test.index, ]


c <- C5.0Control(subset = FALSE, bands = 0, winnow = FALSE, noGlobalPruning = FALSE, CF = .25, minCases = 2, fuzzyThreshold = FALSE, sample = 0, seed = sample.int(4096, size = 1) - 1L, earlyStopping = TRUE)
iris_treeModel <- C5.0(x = iris.train[, -5], y = iris.train$Species, control = c)
summary(iris_treeModel)

test.output <- predict(iris_treeModel, iris.test[, -5], type = "class")
n <- length(test.output)
number = 0
for(i in 1:n) {
        if(test.output[i] == iris.test[i, 5]) {
                number = number + 1
        }
}
test.accuracy = number / n * 100
test.accuracy

iris.train$Species = factor(iris.train$Species) #C5.0的輸出變數資料型態必須是因素，若不是可用factor()函數來轉換


library(C50)
library(stringr)
data(iris)
c <- C5.0Control(subset = FALSE, bands = 0, winnow = FALSE, noGlobalPruning = FALSE, CF = .25, minCases = 2, fuzzyThreshold = FALSE, sample = .9, seed = sample.int(4096, size = 1) - 1L, earlyStopping = TRUE, label = "Species")
iris_treeModel <- C5.0(x = iris[, -5], y = iris$Species, control = c)
summary(iris_treeModel)

x <- str_locate_all(iris_treeModel$output, "%)") #用str_locate_all()取得iris_treeModel這個array的output元素中%的位置
y <- substr(iris_treeModel$output, x[[1]][2]-4, x[[1]][2]-1) #用substr()取得%前1~前4位置的文字
test.error <- as.numeric(y) #用as.numeric()轉換成數字
test.correct <- 100 - test.error #扣掉錯誤率得到正確率
test.correct


#監督式學習 6.2 支持向量機器==================================================

library(e1071)
data(iris)
index <- 1:nrow(iris)
np <- ceiling(.1 * nrow(iris)) #建立測試資料
np
test.index <- sample(1:nrow(iris), np)
iris.test <- iris[test.index, ]
iris.train <- iris[-test.index, ]

svm.model <- svm(Species~., data = iris.train, type = "C-classification", cost = 10, gamma = 10)
svm.pred <- predict(svm.model, iris.test[, -5])
table.svm.test <- table(pred = svm.pred, true = iris.test[, 5])
table.svm.test

correct.svm <- sum(diag(table.svm.test)) / sum(table.svm.test)
correct.svm <- correct.svm * 100
correct.svm

tuned <- tune.svm(Species~., data = iris.train, gamma = 10^(-3:-1), cost = 10^(-1:1)) #用tune.svm()搜尋最佳的cost與gamma組合
summary(tuned)

model <- svm(Species~., data = iris.train, kernel = "radial", gamma = .1, cost = 10)
summary(model)
svm.pred <- predict(model, iris.test[, -5])
table.svm.best.test <- table(pred = svm.pred, true = iris.test[, 5])
table.svm.best.test
correct.svm.best <- sum(diag(table.svm.best.test)) / sum(table.svm.best.test) * 100
correct.svm.best


#監督式學習 6.3 人工神經網路==================================================

install.packages("neuralnet")
library(neuralnet)

traininginput <- as.data.frame(runif(100, min = 0, max = 100))
trainingoutput <- sqrt(traininginput)

trainingdata <- cbind(traininginput, trainingoutput)
colnames(trainingdata) <- c("Input", "Output")

net.sqrt <- neuralnet(Output ~ Input, trainingdata, algorithm = "backprop", hidden = 10, threshold = .01, learningrate = .01, )

print(net.sqrt)
plot(net.sqrt)

testdata <- as.data.frame((1:10)^2) #產生測試資料
net.results <- compute(net.sqrt, testdata) #用compute()來預測模型結果

cleanoutput <- cbind(testdata, sqrt(testdata), as.data.frame(net.results$net.result))
colnames(cleanoutput) <- c("Input", "Expected Output", "Neural Net Output")
print(cleanoutput)


install.packages("DMwR")
library(DMwR)
regr.eval(cleanoutput[, 'Expected Output'],  #求平均絕對誤差MAE、均方根誤差RMSE
          cleanoutput[, 'Neural Net Output'], stats = c("mae", "rmse"))


#監督式學習 6.4 組合方法==================================================


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


#非監督式學習 7.1 階層式分群法==================================================

#非監督式學習 7.2 K平均算法==================================================

#非監督式學習 7.3 模糊C平均算法==================================================

#非監督式學習 7.4 分群指標==================================================


#演化式學習 8.1 基因演算法==================================================


#演化式學習 8.2 人工蜂群演算法==================================================

install.packages("ABCoptim")
library(ABCoptim)

fun <- function(x) {
        -cos(x[1]) * cos(x[2]) * exp(-((x[1] - pi) ^ 2 + (x[2] - pi) ^ 2))
}

abc_optim(rep(0, 2), fun, lb = -20, ub = 20, criter = 1000)


library(ABCoptim)
fw <- function(x) {
        10 * sin(.3 * x) * sin(1.3 * x ^ 2) + .00001 * x ^ 4 + .2 * x + 80
}
abc_optim(50, fw, lb = -100, ub = 100, criter = 1000)


#混合式學習==================================================

#關聯性規則==================================================

#社群網路分析==================================================

#圖形化資料分析工具==================================================

#巨量資料分析（R + Hadoop）==================================================


#Getting and Cleaning Data==================================================

#Reading local flat files

if(!file.exists("data")) {dir.create("data")}
fileUrl <- "http://data.baltimorecity.gov/api/views/dz54-2aru/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "cameras.csv")
dateDownloaded <- date()

cameraData <- read.table("./data/cameras.csv")
head(cameraData)

cameraData <- read.table("./data/cameras.csv", sep = ",", header = TRUE)
head(cameraData)

cameraData <- read.csv("./data/cameras.csv")
head(cameraData)

#Reading Excel files

fileUrl <- "http://data.baltimorecity.gov/api/views/dz54-2aru/rows.xlsx?accessType=DOWNLOAD"
download.file(fileUrl, destfile = "./data/cameras.xlsx")
dateDownloaded <- date()

library(xlsx)
cameraData <- read.xlsx("./data/cameras.xlsx", sheetIndex = 1, header = TRUE)
head(cameraData)

colIndex <- 2:3
rowIndex <- 1:4
cameraDataSubset <- read.xlsx("./data/cameras.xlsx", sheetIndex = 1, colIndex = colIndex, rowIndex = rowIndex)
cameraDataSubset

#Reading XML

library(XML)
fileUrl <- "http://www.w3schools.com/xml/simple.xml"
doc <- xmlTreeParse(fileUrl, useInternal = TRUE)
rootNode <- xmlRoot(doc)
xmlName(rootNode)

rootNode[[1]]
rootNode[[1]][[1]]

xmlSApply(rootNode, xmlValue)

xpathSApply(rootNode, "//name", xmlValue)

fileUrl <- "http://espn.go.com/nfl/team/_/name/bal/baltimore-ravens"
doc <- htmlTreeParse(fileUrl, useInternal = TRUE) #取得所有不同的node
scores <- xpathSApply(doc, "//li[@class='score']", xmlValue)
teams <- xpathSApply(doc, "//li[@class='team-name']", xmlValue)
scores
teams

#Reading JSON

library(jsonlite)
jsonData <- fromJSON("https://api.github.com/users/jtleek/repos")
names(jsonData)
names(jsonData$owner)
jsonData$owner$login

myjson <- toJSON(iris, pretty = TRUE) #把資料轉成JSON，讓某些api用，pretty會編排精美
cat(myjson) #cat就是print()

iris2 <- fromJSON(myjson) #從JSON變data.frame
head(iris2)

#Reading mySQL

install.packages("RMySQL")
library(RMySQL)
ucscDb <-dbConnect(MySQL(), user = "genome", host = "genome-mysql.cse.ucsc.edu")
result <- dbGetQuery(ucscDb, "show databases;")
dbDisconnect(ucscDb) #每次擷取完資料要記得輸入這個！！！
result

hg19 <- dbConnect(MySQL(), user = "genome", db = "hg19", host = "genome-mysql.cse.ucsc.edu")
allTables <- dbListTables(hg19) #dbListTables()
length(allTables)
allTables[1:5]

dbListFields(hg19, "affyU133Plus2") #從hg19資料庫，連上affyU133Plus2表格
dbGetQuery(hg19, "select count(*) from affyU133Plus2") #dbGetQuery()

affyData <- dbReadTable(hg19, "affyU133Plus2") #dbReadTable()
head(affyData)

query <- dbSendQuery(hg19, "select * from affyU133Plus2 where misMatches between 1 and 3") #dbSendQuery
affyMis <- fetch(query) #還沒結束query
quantile(affyMis$misMatches)
affyMisSmall <- fetch(query, n = 10)
dbClearResult(query) #停止query
dim(affyMisSmall)

dbDisconnect(hg19)

#Reading HDF5（groups包含多個data sets及metadata）

source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
library(rhdf5)
created <- h5createFile("example.h5")
created

created <- h5createGroup("example.h5", "foo")
created <- h5createGroup("example.h5", "baa")
created <- h5createGroup("example.h5", "foo/foobaa")
h5ls("example.h5")

A <- matrix(1:10, nr = 5, nc = 2)
h5write(A, "example.h5", "foo/A")
B <- array(seq(0.1, 2.0, by = .1), dim = c(5, 2, 2))
attr(B, "scale") <- "liter"
h5write(B, "example.h5", "foo/foobaa/B")
h5ls("example.h5")

df <- data.frame(1L:5L, seq(0, 1, length.out = 5), 
                 c("ab", "cde", "fghi", "a", "s"), stringsAsFactors = FALSE)
h5write(df, "example.h5", "df")
h5ls("example.h5")
readA <- h5read("example.h5", "foo/A")
readB <- h5read("example.h5", "foo/foobaa/B")
readdf <- h5read("example.h5", "df")
readA

h5write(c(12, 13, 14), "example.h5", "foo/A", index = list(1:3, 1))
h5read("example.h5", "foo/A")

#Reading data from the web

con <- url("http://scholar.google.com/citations?user=HI-I6C0AAAAJ&hl=en")
htmlCode <- readLines(con)
close(con)
htmlCode

library(XML)
url <- "http://scholar.google.com/citations?user=HI-I6C0AAAAJ&hl=en"
html <- htmlTreeParse(url, useInternalNodes = TRUE)
xpathSApply(html, "//title", xmlValue)
xpathSApply(html, "//td[@id='col-citedby']", xmlValue)

library(httr)
html2 <- GET(url)
content2 <- content(html2, as = "text")
parsedHtml <- htmlParse(content2, asText = TRUE)
xpathSApply(parsedHtml, "//title", xmlValue)

pg1 <- GET("http://httpbin.org/basic-auth/user/passwd")
pg1

pg2 <- GET("http://httpbin.org/basic-auth/user/passwd", authenticate("user", "passwd"))
pg2

google <- handle("http://google.com")
pg1 <- GET(handle = google, path = "/")
pg2 <- GET(handle = google, path = "search")

#Reading data from APIs

myapp <- oauth_app("twitter", key = "yourConsumerKeyHere", secret = "yourConsumerSecretHere")
sig <- sign_oauth1.0(myapp, token = "yourTokenHere", token_secret = "yourTokenSecretHere")
homeTL <- GET("http://api.twitter.com/1.1/statuses/home_timeline.json", sig)

json1 <- content(homeTL)
json2 <- jsonlite::fromJSON(toJSON(json1))
json2

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
