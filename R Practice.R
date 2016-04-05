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



#13字串處理==================================================

#14機率分佈==================================================

#15基本統計==================================================

#16線性模型==================================================

#17廣義線性模型==================================================

#18模型診斷==================================================

#19正規化和壓縮方法==================================================

#20非線性模型==================================================

#21時間序列與自相關性==================================================

##自迴歸移動平均模型（Autoregressive Moving Average）

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


##VAR向量自我迴歸
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


##GARCH

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



#08判別分析==================================================



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

#14模型評估與選擇==================================================

install.packages("rattle", dependencies = TRUE)
install.packages("rattle", repos = "http://rattle.togaware.com", type = "source")
library(rattle) #要先安裝GTK+和GGobi
rattle()


###資料挖礦與大數據分析==================================================

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

