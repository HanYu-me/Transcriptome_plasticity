
# 1A -----------------------------------------------------------------------
pdf(file="all figs/fig1/1A.pdf",width=8,height = 6)
hist(ft1016$FT10_mean,col="skyblue",main="",xlab = "Flowering time (days)",
     xlim=c(40,150))
axis(1,seq(40,150,10),c(seq(40,150,10)),las=1,cex.axis=1)
hist(ft1016$FT16_mean,col=rgb(1,0,0,1/4),add=T)
legend("topright", c("10C", "16C"), col=c("skyblue", "pink"), pch=15)
dev.off()


# 1B ----------------------------------------------------------------------
a=fts[ft1016[,1][fts[ft1016[,1],"lat"]>60],1]
b=fts[ft1016[,1][fts[ft1016[,1],"lat"]<60],1]
mean(ft1016[a,"FT10_mean"]-ft1016[a,"FT16_mean"],na.rm=T)
mean(ft1016[b,"FT10_mean"]-ft1016[b,"FT16_mean"],na.rm=T)
pdf(file="all figs/fig1/1B.pdf",width=8,height = 6)
boxplot(ft1016[a,"FT10_mean"],ft1016[a,"FT16_mean"],ft1016[b,"FT10_mean"],ft1016[b,"FT16_mean"],col=c("skyblue","pink","skyblue","pink"),
        main="",names=c("10C north","16C north","10C south","16C south"),ylab="Flowering time (days)",frame.plot=F)

points(x=rep(1,length(a)),y=ft1016[a,"FT10_mean"],pch=20,cex=1,col="black")
points(x=rep(2,length(a)),y=ft1016[a,"FT16_mean"],pch=20,cex=1,col="black")
for(i in 1:length(a)){
  segments(1,ft1016[a,"FT10_mean"][i], 2, ft1016[a,"FT16_mean"][i], col = "tomato", lwd = 1,lty="dashed")
}

points(x=rep(3,length(b)),y=ft1016[b,"FT10_mean"],pch=20,cex=1,col="black")
points(x=rep(4,length(b)),y=ft1016[b,"FT16_mean"],pch=20,cex=1,col="black")
for(i in 1:length(b)){
  segments(3,ft1016[b,"FT10_mean"][i], 4, ft1016[b,"FT16_mean"][i], col = "#5A8FE4", lwd = 1,lty="dashed")
}
dev.off()


# 1C ----------------------------------------------------------------------
p10_z=cbind.data.frame(p10[,1:2],apply(p10[,g1016_qc],2,zscore))
p16_z=cbind.data.frame(p16[,1:2],apply(p16[,g1016_qc],2,zscore))
H <- rep(0,length(g1016_qc))
head(H)
summary(H)
lid <- as.character(c(idnames(mag10),idnames(mag16)))
for(i in 1:length(g1016_qc)){
  y10 <- p10_z[,g1016_qc[i]]
  y16 <- p16_z[,g1016_qc[i]]
  y <- c(y10,y16)
  E <- c(rep(1,nrow(p10)),rep(2,nrow(p16)))
  lm <- lm(y ~ E + lid)
  av <- anova(lm) 
  H[i] <- av$`Sum Sq`[2]/sum(av$`Sum Sq`)
  cat(i,"\n")
}

rep.ibs=ibs(mag10,weight="freq")
rep.ibs[upper.tri(rep.ibs)] <- t(rep.ibs)[upper.tri(rep.ibs)]
#save(rep.ibs,mag10,p10,p16,g1016_qc,file="interdata.RData")
#load("interdata.RData")
hs <- data.frame("p10"=rep(0,length(g1016_qc)),"p16"=rep(0,length(g1016_qc)))
rownames(hs) <- g1016_qc
head(hs)
zscore <- function(x) qnorm((rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))
for(i in 1:length(g1016_qc)){
  cat(i,"\n")
  h10 <-  try(polygenic(zscore(p10[,g1016_qc[i]]),kinship.matrix = rep.ibs,data = mag10),silent = T)
  h16 <-  try(polygenic(zscore(p16[,g1016_qc[i]]),kinship.matrix = rep.ibs,data = mag10),silent = T)
  if(!inherits(h10,"try-error"))
    hs[i,1] <- round(h10$esth2,digits = 2)
  if(!inherits(h16,"try-error"))
    hs[i,2] <- round(h16$esth2,digits = 2)
}
save(hs,file="interdata.RData")
pdf(file="all figs/fig1/1C.pdf",width=8,height = 6)
plot(density(H),xlim=c(-0.1,1),ylim=c(0,7),frame.plot=F,main="",xlab="Heritability")
polygon(density(H),col="#9EEA9E")
lines(density(hs$p10))
polygon(density(hs$p10),col="skyblue")
lines(density(hs$p16))
polygon(density(hs$p16),col=rgb(1,0,0,1/4))
legend("topright", c("H2","h2 10C","h216C"), col=c("lightgreen","skyblue", "pink"), pch=15)
abline(v=median(H),lty="dashed",col="black")
abline(v=median(hs$p10),lty="dashed",col="black")
abline(v=median(hs$p16),lty="dashed",col="black")
dev.off()

median(H)
median(hs$p10)
median(hs$p16)
# 1D ----------------------------------------------------------------------
library(hglm)

a=onestep_pve(idname=ft1016[,1],phe=ft1016[,3],exp=p10_z,genes=g1016_qc,name="aa",filepath = "data/PVE/",nrandom=1)
b=onestep_pve(idname=ft1016[,1],phe=ft1016[,4],exp=p16_z,genes=g1016_qc,name="aaa",filepath = "data\\PVE\\")
d=onestep_pve(idname=ft1016[,1],phe=ft1016[,3]-ft1016[,4],exp=p10_z,genes=g1016_qc,name="aaaa",filepath = "data\\PVE\\")
e=onestep_pve(idname=ft1016[,1],phe=ft1016[,3]-ft1016[,4],exp=p16_z,genes=g1016_qc,name="aaaaa",filepath = "data\\PVE\\")

data=data.frame(rbind(a,b,d,e))
names(data)="X1"
data$group=c("FT 10C","FT 16C","FTP 10C","FTP 16C")
pdf(file="all figs/fig1/1D.pdf",width=8,height = 6)
barplot(data$X1,col="skyblue",names=data$group,ylab="PVE")
dev.off()


