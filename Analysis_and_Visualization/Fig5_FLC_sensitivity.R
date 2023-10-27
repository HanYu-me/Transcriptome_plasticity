

# 5A ----------------------------------------------------------------------
pdf(file="all figs/fig5/5A.pdf",width=5,height=5)
plot(p10[,flc],ft1016[,3],pch=20,main="",col="skyblue",xlab="FLC expression",ylab="Flowering time (days)",xlim=c(0,250),ylim=c(40,150),frame.plot=F,xaxt="n")
axis(1,seq(0,250,50),seq(0,250,50),las=1,cex.axis=1)
abline(lm(ft1016[,3]~p10[,flc]),col="skyblue")
points(p16[,flc],ft1016[,4],pch=20,col="tomato")
abline(lm(ft1016[,4]~p16[,flc]),col="tomato")
legend("topright",fill=c("skyblue","tomato"),legend = c("10C","16C"),bty="n")
dev.off()


# 5B ----------------------------------------------------------------------
pdf(file="all figs/fig5/5B.pdf",width=5,height=5)
plot(p10[,flc],ft1016[,3]-ft1016[,4],pch=20,main="",col="skyblue",xlab="FLC expression",ylab="Flowering time plasticity (days)",xlim=c(0,250),ylim=c(-40,40),frame.plot=F,xaxt="n")
axis(1,seq(0,250,50),seq(0,250,50),las=1,cex.axis=1)
abline(lm(ft1016[,3]-ft1016[,4]~p10[,flc]),col="skyblue")
points(p16[,flc],ft1016[,3]-ft1016[,4],pch=20,col="tomato")
abline(lm(ft1016[,3]-ft1016[,4]~p16[,flc]),col="tomato")
legend("topright",fill=c("skyblue","tomato"),legend = c("10C","16C"),bty="n")
dev.off()


# 5C ----------------------------------------------------------------------
data=data.frame(fts[ft1016$id,],flc10=p10[ft1016$id,flc],flc16=p16[ft1016$id,flc],FTP=ft1016$FT10_mean-ft1016$FT16_mean,ft10=ft1016$FT10_mean,ft16=ft1016$FT16_mean)
pdf(file="all figs/fig5/5C.pdf",width=5,height=5)
boxplot(data$flc10-data$flc16~data$lat<60,names=c("Northern","Southern"),
        col=c("#00BFC4","#F8766D"),xlab="",frame.plot=F,ylab="FLC sensitivity")
t.test(data$flc10-data$flc16~data$lat<60)
dev.off()


# 5D ----------------------------------------------------------------------
pdf(file="all figs/fig5/5D.pdf",width=5,height=5)
plot(data$bio_3,data$flc10-data$flc16,col="#00BFC4",cex=1,pch=20,frame.plot=F,xlim=c(25,35),xaxt="n",xlab="Bio3",ylab="FLC sensitivity")
axis(1,seq(25,35,2),seq(25,35,2),las=1,cex.axis=1)
abline(lm(data$flc10-data$flc16~data$bio_3),col="#00BFC4")
dev.off()
cor.test(data$bio_3,data$flc10-data$flc16,method="spearman")


# 5E ----------------------------------------------------------------------
pdf(file="all figs/fig5/5E.pdf",width=5,height=5)
plot(data$bio_3,data$FTP,col="#00BFC4",cex=1,pch=20,frame.plot=F,xlim=c(25,35),xaxt="n",xlab="Bio3",ylab="Flowering time plasticity")
axis(1,seq(25,35,2),seq(25,35,2),las=1,cex.axis=1)
abline(lm(data$FTP~data$bio_3),col="#00BFC4")
dev.off()
cor.test(data$bio_3,data$FTP,method="pearson")
cor.test(data$bio_3,data$FTP,method="spearman")

# 5F ----------------------------------------------------------------------
pdf(file="all figs/fig5/5F.pdf",width=5,height=5)
boxplot(data$bio_3~data$lat<60,names=c("Northern","Southern"),
        col=c("#00BFC4","#F8766D"),xlab="",frame.plot=F,ylab="Bio3")
dev.off()       
t.test(data$bio_3~data$lat<60)$p.value



# data_analysis -----------------------------------------------------------
moa10[flc,8]
moa16[flc,8]
moa_FTP_10[flc,8]
moa_FTP_16[flc,8]


t.test(data$flc10-data$flc16~data$lat>60)
53.99390/32.06891-1

for(i in 9:19){
  print(paste(round(cor.test(data[,i],data$flc10-data$flc16,method="spearman")$estimate,4),round(cor.test(data[,i],data$flc10-data$flc16,method="spearman")$p.value,4)))
}
