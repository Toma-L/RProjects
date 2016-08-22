# 2. 快速探索數據 =====


## 2.1 散佈圖 =====

library(ggplot2)
qplot(mtcars$wt, mtcars$mpg)
qplot(wt, mpg, data = mtcars)
ggplot(mtcars, aes(x = wt, y = mpg)) + geom_point()


## 2.2 折線圖 =====

plot(pressure$temperature, pressure$pressure, type = "l")

plot(pressure$temperature, pressure$pressure, type = "l")
points(pressure$temperature, pressure$pressure) # 記得是point"s"
lines(pressure$temperature, pressure$pressure/2, col = "red")
points(pressure$temperature, pressure$pressure/2, col = "red")

library(ggplot2)
qplot(pressure$temperature, pressure$pressure, geom = "line") # geom = 
qplot(temperature, pressure, data = pressure, geom = "line") # 方法2
ggplot(pressure, aes(x = temperature, y = pressure)) + geom_line() #方法3

qplot(temperature, pressure, data = pressure, geom = c("line", "point"))


# 2.3 長條圖 =====

barplot(BOD$demand, names.arg = BOD$Time)

barplot(table(mtcars$cyl))

library(ggplot2)
# qplot(BOD$Time, BOD$demand, geom = "bar", stat = "identity") # Error
# qplot(factor(BOD$Time), BOD$demand, geom = "bar", stat = "identity") # Error

# qplot(Time, demand, data = BOD, geom = "bar", stat = "identity") # Error
ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")

qplot(factor(cyl), data = mtcars)
ggplot(mtcars, aes(x = factor(cyl))) + geom_bar()


## 2.4 直方圖 =====

hist(mtcars$mpg)
qplot(mtcars$mpg)

qplot(mpg, data = mtcars, binwidth = 4) # binwidth =


## 2.5 箱型圖 =====

plot(ToothGrowth$supp, ToothGrowth$len) # X為factor，所以自動畫出boxplot
qplot(ToothGrowth$supp, ToothGrowth$len, geom = "boxplot") # another way
ggplot(ToothGrowth, aes(x = supp, y = len)) + geom_boxplot() # the same

# 使用 interaction() 創造更多分組
qplot(interaction(ToothGrowth$supp, ToothGrowth$dose), ToothGrowth$len, geom = "boxplot")
qplot(interaction(supp, dose), len, data = ToothGrowth, geom = "boxplot") # another way
ggplot(ToothGrowth, aes(x = interaction(supp, dose), y = len)) + geom_boxplot() # the same


## 2.6 函數圖形 =====

curve(x ^ 3 - 5 * x, from = -4, to = 4)

myfun <- function(xvar) {
        1 / (1 + exp(-xvar + 10))
}
curve(myfun(x), from = 0, to = 20)
curve(1 - myfun(x), add = TRUE, col = "red")

# qplot(c(0, 20), fun = myfun, stat = "function", geom = "line") # Error
ggplot(data.frame(x = c(0, 20)), aes(x = x)) + stat_function(fun = myfun, geom = "line")


# 3. 長條圖 =====

## 3.1 =====



## 3.2 =====

## 3.3 =====

## 3.4 =====

## 3.5 =====

## 3.6 =====

## 3.7 =====

## 3.8 =====

## 3.9 =====

## 3.10 Cleveland點圖 =====

# Check the name!!!
install.packages("gcookbook")
library(gcookbook)

tophit <- tophitters2001[1:25, ]
ggplot(tophit, aes(x = avg, y = name)) + geom_point()


# 4. 折線圖 =====

# 5. 散佈圖 =====

# 6. 敘述統計 =====

# 6.1 簡單直方圖 =====

ggplot(faithful, aes(x = waiting)) + geom_histogram()

w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_histogram()

ggplot(faithful, aes(x = waiting)) + 
        geom_histogram(binwidth = 5, 
                       fill = "white",  # 填滿
                       colour = "black") # 邊框顏色

binsize <- diff(range(faithful$waiting)) / 15 # 計算binwidth
ggplot(faithful, aes(x = waiting)) + geom_histogram(binwidth = binsize, fill = "white", colour = "black")


h <- ggplot(faithful, aes(x = waiting)) # 儲存成變量以便重複利用

h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 31) # 設定組邊界
h + geom_histogram(binwidth = 8, fill = "white", colour = "black", origin = 35)

# geom_bar(stat = "bin") 可得相同結果


## 6.2 多組資料直方圖 =====

library(MASS) # 取用數據
ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + 
        facet_grid(smoke ~ .)

birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke)
levels(birthwt1$smoke) # 標籤是0, 1不方便

library(plyr) # 為了使用 revalue()
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))

ggplot(birthwt1, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") +
        facet_grid(smoke ~ .)


ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + facet_grid(race ~ .)
ggplot(birthwt, aes(x = bwt)) + geom_histogram(fill = "white", colour = "black") + 
        facet_grid(race ~ ., scales = "free") # scales = "free" 可以單獨設定各y軸尺度

birthwt1$smoke <- factor(birthwt1$smoke)
ggplot(birthwt1, aes(x = bwt, fill = smoke)) + 
        geom_histogram(position = "identity",  # 使直方圖重疊
                       alpha = .4) # 設定透明度


## 6.3 密度曲線 =====

ggplot(faithful, aes(x = waiting)) + geom_density()

ggplot(faithful, aes(x = waiting)) + 
        geom_line(stat = "density") +  # 不想要邊線和底線可以用 geom_line()
        expand_limits(y = 0) # 擴大y軸範圍，包含0


w <- faithful$waiting
ggplot(NULL, aes(x = w)) + geom_density() # 可以傳遞向量作為參數


ggplot(faithful, aes(x = waiting)) +  # 帶寬不同，曲線的光滑程度不同
        geom_line(stat = "density", adjust = .25, colour = "red") +
        geom_line(stat = "density") + 
        geom_line(stat = "density", adjust = 2, colour = "blue")


ggplot(faithful, aes(x = waiting)) + geom_density(fill = "blue", alpha = .2) + 
        xlim(35, 105) # 限制x軸範圍

ggplot(faithful, aes(x = waiting)) + geom_density(fill = "blue", colour = NA, alpha = .2) +
        geom_line(stat = "density") + xlim(35, 105)


ggplot(faithful, aes(x = waiting, y = ..density..)) +  # y = ..density.. 可以縮小直方圖的尺度以配合密度曲線
        geom_histogram(fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        xlim(35, 105)


## 6.4 多組資料密度曲線 =====

library(MASS)
birthwt1 <- birthwt
birthwt1$smoke <- factor(birthwt1$smoke) # 必須先轉化為factor
ggplot(birthwt1, aes(x = bwt, colour = smoke)) + geom_density()

ggplot(birthwt1, aes(x = bwt, fill = smoke)) + geom_density(alpha = .3)


ggplot(birthwt1, aes(x = bwt)) + geom_density() + facet_grid(smoke ~ .)
levels(birthwt1$smoke)
library(plyr)
birthwt1$smoke <- revalue(birthwt1$smoke, c("0" = "No Smoke", "1" = "Smoke"))
ggplot(birthwt1, aes(x = bwt)) + geom_density() + facet_grid(smoke ~ .)

ggplot(birthwt1, aes(x = bwt, y = ..density..)) + 
        geom_histogram(binwidth = 200, fill = "cornsilk", colour = "grey60", size = .2) +
        geom_density() + 
        facet_grid(smoke ~ .)


## 6.5 頻次多邊形 =====

# freqpoly 跟 histogram 非常類似，而核密度曲線只是一個估計

ggplot(faithful, aes(x = waiting)) + geom_freqpoly()

ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = 4) # 控制組距

binsize <- diff(range(faithful$waiting)) / 15
ggplot(faithful, aes(x = waiting)) + geom_freqpoly(binwidth = binsize)


## 6.6 基本箱型圖（盒鬚圖） =====

library(MASS)
ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot()

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot(width = .5) # 調整箱型圖寬度

ggplot(birthwt, aes(x = factor(race), y = bwt)) + 
        geom_boxplot(outlier.size = 1.5, outlier.shape = 21) # 修改outlier的外觀


ggplot(birthwt, aes(x = 1, y = bwt)) + geom_boxplot() + 
        scale_x_continuous(breaks = NULL) +  #只有1組，x軸不標記
        theme(axis.title.x = element_blank())


## 6.7 缺口箱型圖 =====

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot(notch = TRUE)

# 若各組缺口不重疊，表示各組中位數有差異


## 6.8 箱型圖加平均值 =====

ggplot(birthwt, aes(x = factor(race), y = bwt)) + geom_boxplot() + 
        stat_summary(fun.y = "mean", geom = "point", shape = 23, size = 3, fill = "white")


## 6.9 小提琴圖 =====

# install.packages("gcookbook")
library(gcookbook)
p <- ggplot(heightweight, aes(x = sex, y = heightIn))
p + geom_violin()

p + geom_violin() +  # 同時畫小提琴圖跟箱型圖
        geom_boxplot(width = .1, fill = "black", outlier.colour = NA) +
        stat_summary(fun.y = median, geom = "point", fill = "white", shape = 21, size = 2.5)


p + geom_violin(trim = FALSE)

# default情況是不同組的小提琴圖面積會相同，可以用scale = "count"調整

p + geom_violin(scale = "count")

p + geom_violin(adjust = 2) # adjust修改平滑程度
p + geom_violin(adjust = .5)


## 6.10 Wilkinson點圖 =====

library(gcookbook)
countries2009 <- subset(countries, Year == 2009 & healthexp > 2000)

p <- ggplot(countries2009, aes(x = infmortality))
p + geom_dotplot()

p + geom_dotplot(binwidth = .25) + geom_rug() +  # geom_rug() 標示資料點的具體位置
        scale_y_continuous(breaks = NULL) +  # 移除y軸尺度
        theme(axis.title.y = element_blank()) # 移除座標軸標籤
# 數據在水平方向上並非均勻分佈，而是堆在他表示的資料點的中心位置

p + geom_dotplot(method = "histodot", binwidth = .25) + geom_rug() +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())
# 用像直方圖的固定間距分組算法


# 奇數數量與偶數數量保持一致的中心堆疊方式，設定stackdir = "center"或stackdir = "centerwhole"
p + geom_dotplot(binwidth = .25, stackdir = "center") +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())
p + geom_dotplot(binwidth = .25, stackdir = "centerwhole") +
        scale_y_continuous(breaks = NULL) + theme(axis.title.y = element_blank())


## 6.11 多組資料點圖 =====

library(gcookbook)
ggplot(heightweight, aes(x = sex, y = heightIn)) + geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center")

# 同時畫箱型圖跟點圖
ggplot(heightweight, aes(x = sex, y = heightIn)) + geom_boxplot(outlier.colour = NA, width = .4) + 
        geom_dotplot(binaxis = "y", binwidth = .5, stackdir = "center", fill = NA)


# 把箱型圖畫在點圖旁邊
ggplot(heightweight, aes(x = sex, y = heightIn)) + 
        geom_boxplot(aes(x = as.numeric(sex) + .2, group = sex), width = .25) +
        geom_dotplot(aes(x = as.numeric(sex) - .2, group = sex), binaxis = "y", binwidth = .5, stackdir = "center") +
        scale_x_continuous(breaks = 1:nlevels(heightweight$sex), labels = levels(heightweight$sex))


## 6.12 二維資料密度圖 =====

p <- ggplot(faithful, aes(x = eruptions, y = waiting))
p + geom_point() + stat_density2d()


p + stat_density2d(aes(colour = ..level..)) # 把密度曲面的高度映射給等高線的顏色


# 將密度估計映射給填充色或瓦片圖的透明度
# 柵格: raster
# 瓦片: tile
p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
p + geom_point() + stat_density2d(aes(alpha = ..density..), geom = "tile", contour = FALSE)

p + stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE, h = c(.5, 5))


# 7. 註解 =====

# 8. 座標軸 =====

# 9. 整體外觀 =====

# 10. 圖例 =====

# 11. 分面 =====

# 12. 配色 =====

# 13. 其他圖形 =====

## 13.1 相關係數矩陣圖 =====

mcor <- cor(mtcars)
round(mcor, digits = 2)
library(corrplot)
corrplot(mcor)

corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
# corrplot(mcor, method = "shade", shade.col = NA, tl.col = "black", 
#         tl.srt = 45, col = col(200), addCoef.col = "black", 
#         cl.pos = "no", order = "AOE") # Error


# 14. 輸出圖形 =====

# 15. reshape =====