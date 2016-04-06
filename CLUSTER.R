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


