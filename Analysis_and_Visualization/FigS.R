
# S1 ----------------------------------------------------------------------
library(plyr)
library(maptools)
library(sp)

map_info=read.csv("D:\\OneDrive\\桌面\\搞科研\\拟南芥or柳树\\工作文件\\1001_genome_accessions.txt",sep=",",header=F)
load("D:\\OneDrive\\桌面\\搞科研\\拟南芥or柳树\\分析代码\\2022.10.02\\worldmap.RData")
names(map_info)[6:7]=c("lat2","long2")
location_all=map_info[map_info$V1%in%row.names(p10),c(6,7)]
plot(worldMap,xlim=c(0,20),ylim=c(55,65))
points(location_all$long2,location_all$lat2,pch=20,col="red",cex=1)


# S2 ----------------------------------------------------------------------
c_all=c()
for(i in 1:length(g1016_qc)){
  a=cor.test(p10_z[,g1016_qc[i]],p16_z[,g1016_qc[i]],method = "spearman")
  c_all=c(c_all,a$estimate)
}
plot(density(c_all),frame.plot=F,main="",xlab="Spearman coefficient",xlim=c(-1,1))
polygon(density(c_all),col="gray")
abline(v=0)
abline(v=median(c_all),lty="dashed")
text(0.34,3,"Median = 0.16")

summary(c_all)

# S3 ----------------------------------------------------------------------
id=pve_single_FT10/pve_single_FT16<2&pve_single_FT10/pve_single_FT16>0.5
cols=ifelse(id,"gray","skyblue")
plot(pve_single_FT10,pve_single_FT16,pch=20,frame.plot=F,col=cols,
     xlab="Single gene PVE at 16°C",ylab="Single gene PVE at 16°C")
abline(a=0,b=1,lty="dashed",col="black")
abline(a=0,b=2,lty="dashed",col="blue")
abline(a=0,b=0.5,lty="dashed",col="blue")
table(cols)
15140/19007

# S4 ----------------------------------------------------------------------
library(clusterProfiler)
library(org.At.tair.db)
library(ggplot2)
go=enrichGO(gene_list, OrgDb = "org.At.tair.db", keyType="TAIR",ont="ALL")
dotplot(go,split="ONTOLOGY",showCategory=20,font.size=8)+facet_grid(ONTOLOGY~., scale="free")

