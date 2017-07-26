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


## K-平均值分群 ==================================================

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
countries
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