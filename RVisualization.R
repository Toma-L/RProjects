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

## 3.1 簡單長條圖 =====

library(gcookbook)
ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")

ggplot(BOD, aes(x = Time, y = demand)) + geom_bar(stat = "identity")
ggplot(BOD, aes(x = factor(Time), y = demand)) + geom_bar(stat = "identity")

ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", fill = "lightblue", colour = "black")


## 3.2 簇狀條形圖 =====

cabbage_exp
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) +  # fill要用factor變數
        geom_bar(position = "dodge",  # 水平排列
                 stat = "identity")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1") # scale_fill_brewer(palette = "")


ce <- cabbage_exp[1:5, ]
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(position = "dodge", stat = "identity", colour = "black") +
        scale_fill_brewer(palette = "Pastel1")


## 3.3 頻度長條圖 =====

ggplot(diamonds, aes(x = cut)) + geom_bar() # 離散型會得到bar chart

ggplot(diamonds, aes(x = carat)) + geom_bar() # 連續型會得到histogram


## 3.4 長條圖著色 =====

upc <- subset(uspopchange, rank(Change) > 40)
upc

ggplot(upc, aes(x = Abb, y = Change, fill = Region)) + geom_bar(stat = "identity")

ggplot(upc, aes(x = reorder(Abb, Change), y = Change, fill = Region)) + 
        geom_bar(stat = "identity", colour = "black") +
        scale_fill_manual(values = c("#669933", "#FFCC66")) +  # scale_fill_manual()
        xlab("State")


## 3.5 正負長條圖 =====

csub <- subset(climate, Source == "Berkeley" & Year >= 1900)
csub$pos <- csub$Anomaly10y >= 0
csub

ggplot(csub, aes(x = Year, y = Anomaly10y, fill = pos)) + 
        geom_bar(stat = "identity", position = "identity")
# position = "identity" 可以避免系統畫負值而產生的警告訊息

# 但是暖色通常用於正值才對

ggplot(csub, aes(x = Year, y = Anomaly10y, fill = pos)) + 
        geom_bar(stat = "identity", position = "identity", colour = "black", size = .25) + 
        scale_fill_manual(values = c("#CCEEFF", "#FFDDDD"), guide = FALSE)
# size 是調整邊線寬度用的，不是長條寬度！
# guide = FALSE 可以刪除圖例


## 3.6 調整長條寬度和間距 =====

ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity")


ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", width = .5)
ggplot(pg_mean, aes(x = group, y = weight)) + geom_bar(stat = "identity", width = 1)

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = "dodge")

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", width = .5, position = position_dodge(.7))
# position = position_dodge()


## 3.7 堆積長條圖 =====

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity")
cabbage_exp

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") + 
        guides(fill = guide_legend(reverse = TRUE))
# guides() 調整圖例，調整圖例的填充色順序用 fill = guide_legend(reverse = TRUE)

library(plyr) # 為了使用 order = desc()
ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar, order = desc(Cultivar))) + 
        geom_bar(stat = "identity")
# 調整堆疊順序用 order = desc()

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") + 
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")


## 3.8 百分比堆積長條圖 =====

library(plyr)
ce <- ddply(cabbage_exp, "Date", transform,  # 以Date為切割變數進行transform()
            percent_weight = Weight / sum(Weight) * 100)
ce
ggplot(ce, aes(x = Date, y = percent_weight, fill = Cultivar)) + 
        geom_bar(stat = "identity")

cabbage_exp
ce <- ddply(cabbage_exp, "Date", transform,
            percent_weight = Weight / sum(Weight) * 100)

ggplot(ce, aes(x = Date, y = percent_weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")


## 3.9 加上資料標籤 =====

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
        geom_bar(stat = "identity") + 
        geom_text(aes(label = Weight), vjust = 1.5, colour = "white") # 頂端下方

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Weight), vjust = -.2) # 頂端上方
# vjust 調整標籤在長條圖頂端的上或下

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) + 
        geom_bar(stat = "identity") +
        geom_text(aes(label = Weight), vjust = -.2) + 
        ylim(0, max(cabbage_exp$Weight) * 1.05) # 提高y軸上限

ggplot(cabbage_exp, aes(x = interaction(Date, Cultivar), y = Weight)) + 
        geom_bar(stat = "identity") + 
        geom_text(aes(y = Weight + .1, label = Weight)) # 設定標籤的y軸位置

ggplot(cabbage_exp, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = Weight), vjust = 1.5, colour = "white", 
                  position = position_dodge(.9), size = 3) # 分類間距的default為0.9


# 直向堆積長條圖要加入標籤之前，要先對各組資料求出累積和！
library(plyr)
ce <- arrange(cabbage_exp, Date, Cultivar)
ce <- ddply(ce, "Date", transform, label_y = cumsum(Weight))
ce

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity") +
        geom_text(aes(y = label_y, label = Weight), vjust = 1.5, colour = "white")


ce <- arrange(cabbage_exp, Date, Cultivar)
ce <- ddply(ce, "Date", transform, label_y = cumsum(Weight) - .5 * Weight) # 把標籤放在長條圖正中間
ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) +
        geom_bar(stat = "identity") + 
        geom_text(aes(y = label_y, label = Weight), colour = "White")

ggplot(ce, aes(x = Date, y = Weight, fill = Cultivar)) + 
        geom_bar(stat = "identity", colour = "black") +
        geom_text(aes(y = label_y, label = paste(format(Weight, nsmall = 2), "kg")), size = 4) +
        guides(fill = guide_legend(reverse = TRUE)) +
        scale_fill_brewer(palette = "Pastel1")


## 3.10 Cleveland點圖 =====

# Check the name!!!
library(gcookbook)

tophit <- tophitters2001[1:25, ]
ggplot(tophit, aes(x = avg, y = name)) + geom_point()


# 4. 折線圖 =====

## 4.1 簡單折線圖 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line()
BOD

BOD1 <- BOD
BOD1$Time <- factor(BOD1$Time) # 連續型變量轉factor
ggplot(BOD1, aes(x = Time, y = demand, group = 1)) + geom_line() # 對於factor變量，要加入 group = 1 確認為同一組資料


# 拓寬y軸範圍
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + ylim(0, max(BOD$demand))
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + expand_limits(y = 0)


# 4.2 加上資料標記 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point()


library(gcookbook)
ggplot(worldpop, aes(x = Year, y = Population)) + geom_line() + geom_point()
ggplot(worldpop, aes(x = Year, y = Population)) + geom_line() + geom_point() + scale_y_log10()


# 4.3 多重折線圖 =====

library(plyr)
tg <- ddply(ToothGrowth, c("supp", "dose"), summarise, length = mean(len))
ggplot(tg, aes(x = dose, y = length, colour = supp)) + geom_line()
ggplot(tg, aes(x = dose, y = length, linetype = supp)) + geom_line() # 改變線的形狀


# 一定一定要加 group = supp 讓電腦知道資料是一組的
ggplot(tg, aes(x = factor(dose), y = length, colour = supp, group = supp)) + geom_line()

# 錯誤示範
ggplot(tg, aes(x = dose, y = length)) + geom_line() # 沒有正確分組造成一個x對應不只一個點


# 缺頁 =====


# 4.5 修改資料標記 =====

ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point(size = 4, shape = 22, colour = "darkred", fill = "pink")
ggplot(BOD, aes(x = Time, y = demand)) + geom_line() + geom_point(size = 4, shape = 21, fill = "white")


# 缺頁 =====


# 4.7 堆積面積圖 =====

ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + geom_area()

# default 的堆積順序與圖標有時是相反的，可以用 breaks 參數調整
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup)) + geom_area(colour = "black", size = .2, alpha = .4) +
        scale_fill_brewer(palette = "Blues", breaks = rev(levels(uspopage$AgeGroup)))

library(plyr)
# 可以用 order = desc() 反轉堆積順序
ggplot(uspopage, aes(x = Year, y = Thousands, fill = AgeGroup, order = desc(AgeGroup))) +
        geom_area(colour = "black", size = .2, alpha = .4) + 
        scale_fill_brewer(palette = "Blues")


# 缺頁 =====


# 4.9 增加信賴區間 =====

clim <- subset(climate, Source == "Berkeley", select = c("Year", "Anomaly10y", "Unc10y"))
clim
# geom_ribbon() 要先畫，才不會讓 geom_line() 糊掉！
ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_ribbon(aes(ymin = Anomaly10y - Unc10y, ymax = Anomaly10y + Unc10y), alpha = .2) +
        geom_line()

ggplot(clim, aes(x = Year, y = Anomaly10y)) + 
        geom_line(aes(y = Anomaly10y - Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line(aes(y = Anomaly10y + Unc10y), colour = "grey50", linetype = "dotted") +
        geom_line()


# 5. 散佈圖 =====


# 5.1 散佈圖 =====


# 5.2 修改點的樣式 =====

heightweight[, c("sex", "ageYear", "heightIn")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = sex)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + geom_point()

ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex, colour = sex)) + geom_point() + 
        scale_shape_manual(values = c(1, 2)) +  # 人工選擇點的形狀
        scale_colour_brewer(palette = "Set1") # 人工選擇點的顏色


# 5.3 使用非內建的點形 =====

ggplot(heightweight, aes(x = ageYear, y = heightIn)) + geom_point(shape = 3)
ggplot(heightweight, aes(x = ageYear, y = heightIn, shape = sex)) + geom_point(size = 3) + scale_shape_manual(values = c(1, 4))


hw <- heightweight
hw$weightGroup <- cut(hw$weightLb, breaks = c(-Inf, 100, Inf), labels = c("< 100", ">= 100"))
ggplot(hw, aes(x = ageYear, y = heightIn, shape = sex, fill = weightGroup)) + 
        geom_point(size = 2.5) +
        scale_shape_manual(values = c(21, 24)) +
        scale_fill_manual(values = c(NA, "black"), 
                          guide = guide_legend(override.aes = list(shape = 21)))


# 5.4 將連續型變數映射到點的顏色或大小 =====

heightweight[, c("sex", "ageYear", "heightIn", "weightLb")]

ggplot(heightweight, aes(x = ageYear, y = heightIn, colour = weightLb)) + geom_point()
ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb)) + geom_point()

# 人天生對顏色和大小的變化不太敏銳，因此當一個變數不需要太精密的解釋時才適合用

ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) + 
        scale_fill_gradient(low = "black", high = "white") # 黑白漸層

ggplot(heightweight, aes(x = ageYear, y = heightIn, fill = weightLb)) + 
        geom_point(shape = 21, size = 2.5) +
        scale_fill_gradient(low = "black", high = "white", breaks = seq(70, 170, by = 20), guide = guide_legend())
# 黑白漸層搭配離散圖標

ggplot(heightweight, aes(x = ageYear, y = heightIn, size = weightLb, colour = sex)) + 
        geom_point(alpha = .5) +
        scale_size_area() +  # 使面積與資料值成正比
        scale_colour_brewer(palette = "Set1")

# 不同形狀的點很難比較面積大小，因此不要同時操作


# 5.5 處理資料點重疊問題 =====

# 圖形重疊（overplotting）的解決方案：
## 半透明的點
## 矩形資料分箱
## 六邊形資料分箱
## 箱型圖

sp <- ggplot(diamonds, aes(x = carat, y = price))
sp + geom_point()

sp + geom_point(alpha = .1) # 90%透明度
sp + geom_point(alpha = .01) # 99%透明度


sp + stat_bin2d() # 分別在x軸y軸分割30組，共900個箱子
sp + stat_bin2d(bins = 50) + scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 6000)) # 2500箱


library(hexbin) # 使用六邊形箱子
sp + stat_binhex() + scale_fill_gradient(low = "lightblue", high = "red", limits = c(0, 8000))
sp + stat_binhex() + scale_fill_gradient(low = "lightblue", high = "red", breaks = c(0, 250, 500, 1000, 2000, 4000, 6000), limits = c(0, 6000))
# 範圍外會變成灰色箱子


# 當其中一軸或兩軸為離散型變數時，也會出現overplotting，可用 position_jitter() 增加隨機擾動

sp1 <- ggplot(ChickWeight, aes(x = Time, y = weight))
sp1 + geom_point()
sp1 + geom_point(position = "jitter") # 跟 position_jitter() 意思一樣

sp1 + geom_point(position = position_jitter(width = .5, height = 0)) # width 和 height 調整擾動值的精度




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