library(plyr)
library(maptools)
library(sp)
library(dismo)
library(raster)
library(sp)
library(rgdal)
library(gstat)
library(raster)
library(maptools)
library(RColorBrewer)
library(grid)
library(lattice)
library(maptools)
library(dplyr)
library(reshape2)
library(Hmisc)
library(data.table)



tif_path="D:/OneDrive/桌面/搞科研/拟南芥or柳树/工作数据/2022.10.12_15/wc2.1_10m_bio_19"
name <- strsplit(tif_path,split="/")
number <- length(name[[1]]) #判断目录级数
lst <- list.files(path=tif_path,pattern='tif$',full.names = T) # 读取目录下所有指定类型的文件，并返回列表。这里用的是tiff格式图层，可根据自己的需要修改，例如# bil、asc 等。数据来源WorldClid网站。
n <- length(lst) #检测数据文件个数
file<-lst[13]
file
BIO <- raster(file) #读取第i个气候图层
BIO_centroid <- extract(BIO, data, method='simple', fun=mean, df=TRUE) #提取每个点的气候信息，第i个气候图层，其中用的是均值(座标处不一定有数值)，同时也可选择其他计算方法。
result<-BIO_centroid[,2]

location_all=acc[ft1016$id,c("long","lat")]
pdf("all figs/fig7/7loaction.pdf",width = 8,height = 6)
plot(BIO,xlim=c(0,20),ylim=c(55,65),
     col=colorRampPalette(c("green","lightgreen","yellow","red"))(100))
points(location_all[,1],location_all[,2],pch=17,col="gray",cex=1)
points(location_all[,1],location_all[,2],pch=2,col="black",cex=1)
dev.off()