#資料挖礦與大數據分析==========

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


#利用R語言打通大數據的經脈==================================================

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


##K-平均值分群==================================================

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


##系譜分群（階層式）==================================================

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


#密度分群==================================================

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


#期望值最大分群==================================================

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


#R軟體資料分析基礎與應用==================================================

#22資料分群==================================================




#Practical Machine Learning==================================================

#Unsupervised prediction==================================================

data(iris)
library(ggplot2)
inTrain <- createDataPartition(y = iris$Species, p = .7, list = FALSE)
training <- iris[inTrain, ]
testing <- iris[-inTrain, ]
dim(training); dim(testing)

kMeans1 <- kmeans(subset(training, select = -c(Species)), centers = 3)
training$clusters <- as.factor(kMeans1$cluster)
qplot(Petal.Width, Petal.Length, colour = clusters, data = training)

table(kMeans1$cluster, training$Species)
modFit <- train(clusters ~., data = subset(training, select = -c(Species)), method = "rpart")
table(predict(modFit, training), training$Species)

testClusterPred <- predict(modFit, testing)
table(testClusterPred, testing$Species)