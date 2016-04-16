#Factor Analysis

r.no <- factanal(sca[-c(3, 10, 13)], factors = 2, data = sca, rotation = "none")
r.no$loadings

r.va <- factanal(sca[-c(3, 10, 13)], factors = 2, data = sca, rotation = "varimax") #最大變異法進行正交轉軸
r.va$loadings

plot(r.no$loadings[, 1], r.no$loadings[, 2], pch = 16, cex = 1.6, col = 2)
grid(nx = 10, ny = 10)

plot(r.va$loadings[, 1], r.va$loadings[, 2], pch = 16, cex = 1.6, col = 4)
grid(nx = 10, ny = 10)


##因素分析函數語法==============================

###factanal()函數

#萃取共同因素的方法：主成份分析法/未加權最小平方法/一般化最小平方法/最大概似估計法/主軸因子法/最小化樣本加權卡方估計法/加權最小平方法
#正交轉軸法：最大變異法/四方最大法/bentlerT/均等最大法
#斜交轉軸法：最優斜交法/直接斜交法/簡單最大法
#因素分數型態：none/regression(Thompson分數)/Barlett加權最小平方分數
#要輸出何種因素矩陣：結構矩陣/樣式矩陣（組型矩陣）


##基本套件factanal()函數的應用==============================

temp <- read.csv("efa_0.csv", header = TRUE)
head(temp)
fac.m <- factanal(temp, factors = 2) #因速抽取方法為「最大概似估計法」
print(fac.m, digits = 2)
#Uniquenesses：變數無法被共同因素解釋的誤差變異
#Loadings：轉軸後的成分矩陣
#SS loadings：共同因素的特徵值
#Proportion Var：因素解釋提項變數的變異量
#Cumulative Var：共同因素可解釋指標變數的總變異量，高於60%為佳

#檢定共同因素是否足夠，p值0.162 > .05，接受虛無假設，抽取隻因素個數足夠反應題項面數的所有構面

round(fac.m$loadings[, 1], 2)
round(fac.m$loadings[, 2], 2)

sum(fac.m$loadings[, 1]^2) #共同因素的特徵值
sum(fac.m$loadings[, 2]^2)

sum(fac.m$loadings[1, ]^2) #第一個題項變數的共同性，指變數能被共同因素解釋的變異部分，越大表示變數越能有效反應潛在構念
for(i in 1:7) {
        h2 <- round(sum(fac.m$loadings[i, ]^2), 2)
        u2 = 1 - h2
        print(paste("h2=", h2, " u2=", u2))
}

fac.m$factors #因素個數
fac.m$method #最大概似估計法
round(fac.m$correlation, 2) #題項變數間的相關矩陣


##principal()函數的應用==============================

install.packages("psych")
install.packages("GPArotation")
library(psych)
library(GPArotation)


###直交轉軸

efa.v <- principal(temp, nfactors = 2, rotate = "varimax") #主成份分析法抽取，最大變異法轉軸
print(efa.v)

efa.v$loadings[, 1]
sum((efa.v$loadings[, 1])^2)
efa.v$loadings[, 2]
sum((efa.v$loadings[, 2])^2)

print(efa.v, cut = .45, sort = TRUE) #因素負荷量<.45的值不會輸出，sort = T表示共同因素依照因素負荷量絕對值大小排列

names(efa.v)

round(efa.v$values, 3) #共同因素的特徵值
sum(efa.v$values) #七個成分的特徵值總和等於題項個數

sum(efa.v$factors)
round(efa.v$scores, 3) #樣本觀察值在兩個共同因素的分數

efa.v$loading #因素負荷量

(efa.v$loading[1, ])^2 #因素負荷量平方
round(sum((efa.v$loading[1, ])^2), 3) #共同性

efa.v$rotation #使用何種轉軸法

efa.v <- principal(temp, nfactors = 2, rotate = "varimax")
efa.n <- principal(temp, nfactors = 2, rotate = "none")

n.varx <- round(efa.n$loadings[, 1], 2)
n.vary <- round(efa.n$loadings[, 2], 2)
plot(n.varx, n.vary, cex = 2, pch = 15, col = 4) #未轉軸前因素結構圖
grid(nx = 10, ny = 10)

#可看出七個指標變數很難分成兩群

n.varx <- round(efa.v$loadings[, 1], 2)
n.vary <- round(efa.v$loadings[, 2], 2)
plot(n.varx, n.vary, cex = 2, pch = 15, col = 4) #轉軸後因素結構圖
grid(nx = 9, ny = 9)

#七個指標變數容易分成兩群

itn <- names(n.varx)
itn
itn <- substring(itn, 1, 1) #簡化題項變數名稱
itn

v.varx <- round(efa.v$loadings[, 1], 2)
v.vary <- round(efa.v$loadings[, 2], 2)
plot(v.varx, v.vary, cex = 2, pch = 16, col = 4)
text(v.varx, v.vary, itn, pos = 2, offset = 1.5, font = 2, cex = )
grid(nx = 10, ny = 10)


###斜交轉軸

efa.op <- principal(temp, nfactor = 2, 
                    rotate = "oblimin", oblique.scores = TRUE) #主成份法+直接斜交
print(efa.op) #會呈現共同因素間的相關矩陣

round(efa.v$communality, 3) #共同性估計值


ml.p <- fa(temp, nfactors = 2, rotate = "promax", #最大概似估計法+最優斜交
           oblique.scores = FALSE, fm = "ml")
print(ml.p)


##人格特質量表的EFA==============================

temp <- read.csv("efa_1.csv", header = TRUE)
efa <- temp[, 2:18]
names(efa)

library(psych)
library(GPArotation)

###量表共同因素的檢核

KMO(efa) #KMO看量表題項變數是否適合進行因素分析，需達0.60以上，越接近1越好，又稱「抽樣適切性量測值」
#MSA是取樣適切性量數，代表該指標變數與其他指標變數間的相關程度，越小越無法反應共同的潛在因素

cortest.bartlett(efa) #Bartlett球形檢定為變數間淨相關係數的檢核，顯著表示指標變數間「有共同因素存在」

###因素個數的判別

pc.no <- principal(efa, nfactors = 3, rotate = "none")
print(pc.no, digits = 2) #拒絕虛無假設，萃取的因素個數不足夠

round(pc.no$values, 2) #因素分析物件的特徵值，保留大於1的（有4個）
eig <- round(pc.no$values, 2) 
plot(eig, pch = 16, cex = 1.5) #繪製陡坡圖
lines(eig, col = 4, lwd = 2)
grid(nx = 17, ny = 10)

pr4.no <- principal(efa, 4, rotate = "none") #主成份法+不轉軸
print(pr4.no, digits = 2)


###因素萃取方法為主成份，轉軸方法為最大變異法

principal(efa, 4, rotate = "varimax")

pc.m <- principal(efa, nfactors = 4, rotate = "varimax")
print(pc.m, cut = .40, sort = TRUE) #因素負荷量小於.40的不會呈現


###刪除第二題（DA2）的EFA

efa.1 <- efa[, -2]
names(efa.1)