# 2ABC selsetion single gene pve ------------------------------------------------------------
##estimate selestion
s <- data.frame("r10"=rep(0,length(g1016_qc)),"p10"=rep(0,length(g1016_qc)),
                "r16"=rep(0,length(g1016_qc)),"p16"=rep(0,length(g1016_qc)))
for(i in 1:length(g1016_qc)){
  est_S <- function(x,y){
    y<- y/mean(y,na.rm=T)
    ct10 <- lm(y~scale(x))
    #plot(x,y)
    #abline(ct10)
    sct10 <- summary(ct10)
    return(c(sct10$coefficients[2,1],sct10$coefficients[2,4]))
  }
  s[i,c(1,2)] <- est_S(x =p10[idnames(mag10),g1016_qc[i]] ,y = ft1016[idnames(mag10),"FT10_mean"])
  s[i,c(3,4)] <- est_S(x =p16[idnames(mag10),g1016_qc[i]] ,y = ft1016[idnames(mag10),"FT16_mean"])
  cat(i,"\n")
}
##single gene pve
pve_single_FT10=c()
pve_single_FT16=c()
pve_single_FTP10=c()
pve_single_FTP16=c()
for(i in 1:length(g1016_qc)){
  setTxtProgressBar(pb, i/length(g1016_qc))
  pve_single_FT10=c(pve_single_FT10,lm_fixed(p10_z,g1016_qc[i],ft1016[,3],idname=""))
  pve_single_FT16=c(pve_single_FT16,lm_fixed(p16_z,g1016_qc[i],ft1016[,4],idname=""))
  pve_single_FTP10=c(pve_single_FTP10,lm_fixed(p10_z,g1016_qc[i],ft1016[,3]-ft1016[,4],idname=""))
  pve_single_FTP16=c(pve_single_FTP16,lm_fixed(p16_z,g1016_qc[i],ft1016[,3]-ft1016[,4],idname=""))
}

##plot
png(file="all figs/fig2/2ABC.png",width=4,height = 9.6,units="in",res=3000)
par(mfrow=c(3,1))
par(mar = c(2, 5, 2, 5))
s3 <- s[,3]
s1 <- s[,1]
plot(density(s1),lwd=2,xlim=c(-0.3,0.3),col = rgb(135/255,206/255,235/255,1),xlab="",ylab="",frame.plot=F,main="")
lines(density(s3),col = rgb(255/255,99/255,71/255,1),lwd=2)
abline(v=0,col="black",lty="dashed")
legend("topright",col=c("skyblue","tomato"),legend = c("",""),lty=1,lwd=2,bty="n",inset = 0.07)


plot(s[,1],s[,3],pch=19,cex=0.5,col="gray",xlab="",ylab="",frame.plot=F,
     xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
points(s[s$p10<0.01,1],s[s$p10<0.01,3],pch=19,cex=0.5,col="skyblue",xlab=expression("S"["10C"]),ylab=expression("S"["16C"]),frame.plot=F,
       xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
points(s[s$p16<0.01,1],s[s$p16<0.01,3],pch=19,cex=0.5,col="tomato",xlab=expression("S"["10C"]),ylab=expression("S"["16C"]),frame.plot=F,
       xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
points(s[s$p10<0.01&s$p16<0.01,1],s[s$p16<0.01&s$p10<0.01,3],pch=19,cex=0.5,col="lightgreen",xlab=expression("S"["10C"]),ylab=expression("S"["16C"]),frame.plot=F,
       xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
#legend("topright",fill=c("skyblue","tomato","lightgreen"),legend = c("10C","16C","both"),bty="n",inset = 0.07)
abline(lm(s[,3]~s[,1]),lty="dashed",col="tomato")
abline(v=0,lty="dashed")
abline(h=0,lty="dashed")



plot(density(pve_single_FT10),col="skyblue",frame.plot=F,lwd=2,xlim=c(-0.01,0.3),
     main="",xlab="",ylab="")
lines(density(pve_single_FT16),col="tomato",lwd=2)
#lines(density(pve_single_FTP10),col="purple",lwd=2)
#lines(density(pve_single_FTP16),col="darkorange",lwd=2)
legend("topright",col=c("skyblue","tomato"),
       legend = c("",""),lty=1,lwd=2,bty="n",inset = 0.07)
# legend("topright",fill=c("skyblue","tomato","purple","darkorange"),
#        legend = c("FT 10C","FT 16C","FTP 10C","FTP 16C"),bty="n",inset = 0.07)
dev.off()

median(abs(s1))
max(abs(s1))
median(abs(s3))
a=g1016_qc[s$p10<0.01]
b=g1016_qc[s$p16<0.01]
c=g1016_qc[s$p10<0.01&s$p16<0.01]
length(setdiff(unique(c(a,b)),intersect(a,b)))
length(intersect(a,b))

median(pve_single_FT10)
median(pve_single_FT16)
# 2D manhattan ------------------------------------------------------------
##write exp data
exp=p10_z
exp[,1:2]=exp$id
names(exp)[1:2]=c("FID","IID")
write.table(exp,file="final_data/twa/exp_10C.txt",row.names = F,col.names = T,quote=F)
exp=p16_z
exp[,1:2]=exp$id
names(exp)[1:2]=c("FID","IID")
write.table(exp,file="final_data/twa/exp_16C.txt",row.names = F,col.names = T,quote=F)
##write phe data
names(ft1016)
write.table(ft1016[,c(1,1,3)],file="final_data/twa/FT_10C.txt",row.names = F,col.names = T,quote=F)
write.table(ft1016[,c(1,1,4)],file="final_data/twa/FT_16C.txt",row.names = F,col.names = T,quote=F)
write.table(cbind.data.frame(ft1016[,c(1,1)],ft1016$FT10_mean-ft1016$FT16_mean),file="final_data/twa/FTP.txt",row.names = F,col.names = T,quote=F)

###TWA use my shiny app
moa10=data.frame(fread("final_data/twa/moa_FT10.moa"))
moa16=data.frame(fread("final_data/twa/moa_FT16.moa"))
moa_FTP_10=data.frame(fread("final_data/twa/moa_FTP10.moa"))
moa_FTP_16=data.frame(fread("final_data/twa/moa_FTP16.moa"))
row.names(moa10)=moa10$Probe
row.names(moa16)=moa16$Probe
row.names(moa_FTP_10)=moa_FTP_10$Probe
row.names(moa_FTP_16)=moa_FTP_16$Probe

gene=data.frame(sub("gene_id (\\S+) .+", "\\1", gsub(";", "", gtf$V9[which(gtf$V3=="gene")])))
gene$chr=gtf$V1[which(gtf$V3=="gene")]
gene$BP=(gtf$V4[which(gtf$V3=="gene")]+gtf$V5[which(gtf$V3=="gene")])/2
names(gene)[1]="Probe"
row.names(gene)=gene$Probe
chrs=c(0,36513169,62296673,91841993,116512553)
for(i in 1:nrow(gene)){
  setTxtProgressBar(pb, i/nrow(gene))
  gene[i,3]=gene[i,3]+chrs[gene[i,2]]
}
##plot
d=0.05
FT_10=row.names(moa10[p.adjust(moa10$p,method="fdr")<d,])
FT_16=row.names(moa16[p.adjust(moa16$p,method="fdr")<d,])
FTP_10=row.names(moa_FTP_10[p.adjust(moa_FTP_10$p,method="fdr")<d,])
FTP_16=row.names(moa_FTP_16[p.adjust(moa_FTP_16$p,method="fdr")<d,])



cols=ifelse(gene[row.names(moa10),2]%in%c(2,4),"tomato","skyblue")

png(file="all figs/fig2/2D.png",width=8,height = 6.4,units="in",res=3000)
plot(c(min(gene$BP,na.rm=T),max(gene$BP,na.rm=T)),c(0,max(-log10(moa10[,8]))),
     xaxt="n",yaxt="n",frame.plot=F,ylab="",xlab="",type = "n")
chr=c(chrs[1]+15e6,chrs[2]+10e6,chrs[3]+13e6,chrs[4]+10e6,chrs[5]+15e6)
axis(1,chr,c("1","2","3","4","5"),las=1,cex.axis=1)
axis(2,c(0,2,4,6,8,10,12,14,16,18,20),c(0,2,4,6,8,10,12,14,16,18,20),las=1,cex.axis=1)
points(gene[row.names(moa10),3],-log10(moa10[,8]),col=cols,pch=20,cex=0.6)
points(gene[row.names(moa16),3],-log10(moa16[,8]),col=cols,pch=20,cex=0.6)
points(gene[row.names(moa_FTP_10),3],-log10(moa_FTP_10[,8]),col=cols,pch=20,cex=0.6)
points(gene[row.names(moa_FTP_16),3],-log10(moa_FTP_16[,8]),col=cols,pch=20,cex=0.6)
#abline(h=-log10(0.05/length(g1016_qc)),col="red")
abline(h=-log10(max(moa10[p.adjust(moa10$p,method="fdr")<d,8])),lty="dashed",col="tomato")
abline(h=-log10(max(moa16[p.adjust(moa16$p,method="fdr")<d,8])),lty="dashed",col="skyblue")
abline(h=-log10(max(moa_FTP_10[p.adjust(moa_FTP_10$p,method="fdr")<d,8])),lty="dashed",col="purple")
abline(h=-log10(max(moa_FTP_16[p.adjust(moa_FTP_16$p,method="fdr")<d,8])),lty="dashed",col="darkorange")
points(gene[FT_10,3],-log10(moa10[FT_10,8]),pch=17,cex=0.8,col="tomato")
points(gene[FT_16,3],-log10(moa16[FT_16,8]),pch=17,cex=0.8,col="skyblue")
points(gene[FTP_10,3],-log10(moa_FTP_10[FTP_10,8]),pch=17,cex=0.8,col="purple")
points(gene[FTP_16,3],-log10(moa_FTP_16[FTP_16,8]),pch=17,cex=0.8,col="darkorange")
legend("topright",col=c("skyblue","tomato","purple","darkorange"),
       legend = c("FT10","FT16","FTP10","FTP16"),bty="n",pch=17)
abline(v=gene[FTgenes[FTgenes%in%gene_list],3],col="black",lty="dashed")
dev.off()

a=gene_list
gene_list=unique(c(FT_10,FT_16,FTP_10,FTP_16))

FT_gene=unique(c(FT_10,FT_16))
FTP_gene=unique(c(FTP_10,FTP_16))

# 2E ----------------------------------------------------------------------
library(pheatmap)
data=data.frame(moa10[gene_list,"p"],moa16[gene_list,"p"],moa_FTP_10[gene_list,"p"],moa_FTP_16[gene_list,"p"])
row.names(data)=gene_list
names(data)=c("FT 10C","FT 16C","FTP 10C","FTP 16C")
data=-log10(data)
a=setdiff(unique(c(FT_10,FT_16)),intersect(unique(c(FT_10,FT_16)),unique(c(FTP_10,FTP_16))))
b=intersect(unique(c(FT_10,FT_16)),unique(c(FTP_10,FTP_16)))
c=setdiff(unique(c(FTP_10,FTP_16)),intersect(unique(c(FT_10,FT_16)),unique(c(FTP_10,FTP_16))))
data=rbind(rbind(data[a,],data[b,]),data[c,])
pdf(file="all figs/fig2/2E.pdf",width=8,height = 3)
pheatmap(t(data),cluster_cols = F,cluster_rows = F,col=colorRampPalette(c("white","orange","red"))(100))
dev.off()




# 2F ----------------------------------------------------------------------
library(VennDiagram)
ven_p=venn.diagram(
  x=list(
    "10C FT"=FT_10,"16C FT"=FT_16,
    "10C FTP"=FTP_10,"16C FTP"=FTP_16,"FTgene"=FTgenes
  ),
  filename = NULL,    #保存路径
  col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明          
  fill = c("blue", "green", "tomato","pink","orange"),  #填充颜色
  alpha = 0.50,     #透明度
  cex = 1,    #每个区域label名称的大小
  cat.col = c("blue", "green", "tomato","pink","orange"),  #分类颜色
  reverse=TRUE,
  main.just = c(2, 1)
)
grid.draw(ven_p)

ven_p=venn.diagram(
  x=list(
    "10C FT"=FT_10,"16C FT"=FT_16,
    "10C FTP"=FTP_10,"16C FTP"=FTP_16
  ),
  filename = NULL,    #保存路径
  col = "transparent",      #指定图形的圆周边缘颜色  transparent 透明          
  fill = c("skyblue", "tomato", "#E2C7E7","#EEC296"),  #填充颜色
  alpha = 0.50,     #透明度
  cex = 1,    #每个区域label名称的大小
  cat.col = c("skyblue", "tomato", "#E2C7E7","#FFA143"),  #分类颜色
  reverse=TRUE,
  main.just = c(2, 1)
)
grid.draw(ven_p)


#  2G ---------------------------------------------------------------------
a=onestep_pve(phe=ft1016[,3],exp=p10_z,genes=setdiff(gene_list,FT_10),name="aa",filepath = "data/PVE/",nrandom=1)
b=onestep_pve(phe=ft1016[,4],exp=p16_z,genes=setdiff(gene_list,FT_16),name="aaa",filepath = "data\\PVE\\")
d=onestep_pve(phe=ft1016[,3]-ft1016[,4],exp=p10_z,genes=setdiff(gene_list,FTP_10),name="aaaa",filepath = "data\\PVE\\")
e=onestep_pve(phe=ft1016[,3]-ft1016[,4],exp=p16_z,genes=setdiff(gene_list,FTP_16),name="aaaaa",filepath = "data\\PVE\\")

data=data.frame(rbind(a,b,d,e))
names(data)="X1"
data$group=c("FT10","FT16","FTP10","FTP16")
pdf(file="all figs/fig2/2G.pdf",width=4,height = 5)
barplot(data$X1,col=c("skyblue","tomato","#E2C7E7","#FFC285"),names=data$group,ylab="PVE")
axis(1,c(1,2,3,4),c("","","",""),las=1,cex.axis=1)
dev.off()


data=data.frame(gene_list,moa10[gene_list,7:8],moa16[gene_list,7:8],moa_FTP_10[gene_list,7:8],moa_FTP_16[gene_list,7:8])
write.csv(data,file="test.csv",row.names = F,quote=F)
