
# Get methylation of each geen --------------------------------------------
# get gene position and strand --------------------------------------------
gene=data.frame(sub("gene_id (\\S+) .+", "\\1", gsub(";", "", gtf$V9[which(gtf$V3=="gene")])))
gene$chr=gtf$V1[which(gtf$V3=="gene")]
gene$start=gtf$V4[which(gtf$V3=="gene")]
gene$end=gtf$V5[which(gtf$V3=="gene")]
gene$strand=gtf$V7[which(gtf$V3=="gene")]
names(gene)[1]="Probe"
row.names(gene)=gene$Probe
gene=na.omit(gene[g1016_qc,])
write.table(gene,file="D:\\Methylation\\gene_position.txt",col.names = T,row.names = F,quote = F)
library(stringr)
gene=read.table(file="D:\\Methylation\\gene_position.txt",header = T,stringsAsFactors = F)
row.names(gene)=gene$Probe
gene$chr=paste0("chr",gene$chr)
gene$start=gene$start-3000
gene$end=gene$end+3000


#BiocManager::install("methylKit")

library(methylKit)
region=GRanges(seqnames=gene$chr,ranges=IRanges(start=gene$start,end=gene$end),strand=gene$strand)

all_gene=paste0(gene$chr,"_",gene$start,"_",gene$strand)

# CG ----------------------------------------------------------------------
files= list.files("D:\\Methylation\\Meth_types\\CG\\")
path="D:\\Methylation\\Meth_types\\CG\\"
gene_CG_10=gene
gene_CG_16=gene
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  #10s a data
  name=str_split(files[i],"_")[[1]]
  if(length(name)==3){
    if(substring(name[3],1,3)=="10C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpG", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CG=regionCounts(object=obj,regions=region)
      subname=paste0(CG$chr,"_",CG$start,"_",CG$strand)
      row.names(CG)=subname
      CG=getData(CG)
      gene_CG_10[,name]=CG[all_gene,"numCs"]/CG[all_gene,"coverage"]
    }else if(substring(name[3],1,3)=="16C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpG", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CG=regionCounts(object=obj,regions=region)
      subname=paste0(CG$chr,"_",CG$start,"_",CG$strand)
      row.names(CG)=subname
      CG=getData(CG)
      gene_CG_16[,name]=CG[all_gene,"numCs"]/CG[all_gene,"coverage"]
    }
  }
}
save(gene_CG_10,gene_CG_16,file="data/Methylome/CG.RData")
# CHG ---------------------------------------------------------------------
files= list.files("E:\\Methylation\\Meth_types\\CHG\\")
path="D:\\Methylation\\Meth_types\\CHG\\"
gene_CHG_10=gene
gene_CHG_16=gene
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  #10s a data
  name=str_split(files[i],"_")[[1]]
  if(length(name)==3){
    if(substring(name[3],1,3)=="10C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHG=regionCounts(object=obj,regions=region)
      subname=paste0(CHG$chr,"_",CHG$start,"_",CHG$strand)
      row.names(CHG)=subname
      CHG=getData(CHG)
      gene_CHG_10[,name]=CHG[all_gene,"numCs"]/CHG[all_gene,"coverage"]
    }else if(substring(name[3],1,3)=="16C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHG=regionCounts(object=obj,regions=region)
      subname=paste0(CHG$chr,"_",CHG$start,"_",CHG$strand)
      row.names(CHG)=subname
      CHG=getData(CHG)
      gene_CHG_16[,name]=CHG[all_gene,"numCs"]/CHG[all_gene,"coverage"]
    }
  }
}
save(gene_CHG_10,gene_CHG_16,file="data/Methylome/CHG.RData")
# CHH ---------------------------------------------------------------------
files= list.files("D:\\Methylation\\Meth_types\\CHH\\")
path="D:\\Methylation\\Meth_types\\CHH\\"
gene_CHH_10=gene
gene_CHH_16=gene
for(i in 208:length(files)){
  if(i%%20==0){
    gc()
  }
  setTxtProgressBar(pb, i/length(files))
  #90s a data
  name=str_split(files[i],"_")[[1]]
  if(length(name)==3){
    if(substring(name[3],1,3)=="10C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CHH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHH=regionCounts(object=obj,regions=region)
      subname=paste0(CHH$chr,"_",CHH$start,"_",CHH$strand)
      row.names(CHH)=subname
      CHH=getData(CHH)
      gene_CHH_10[,name]=CHH[all_gene,"numCs"]/CHH[all_gene,"coverage"]
    }else if(substring(name[3],1,3)=="16C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CHH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHH=regionCounts(object=obj,regions=region)
      subname=paste0(CHH$chr,"_",CHH$start,"_",CHH$strand)
      row.names(CHH)=subname
      CHH=getData(CHH)
      gene_CHH_16[,name]=CHH[all_gene,"numCs"]/CHH[all_gene,"coverage"]
    }
  }
}
save(gene_CHH_10,gene_CHH_16,file="data/Methylome/CHH.RData")

names=data.frame(str_split_fixed(files,"_",3))
names=names[setdiff(1:nrow(names),grep("rep",names$X3)),]
names(table(names$X2)[table(names$X2)==2])

load("data/Methylome/CG.RData")
load("data/Methylome/CHG.RData")
load("data/Methylome/CHH.RData")

data_10=data.frame(t(gene_CG_10),stringsAsFactors = F)
data_16=data.frame(t(gene_CG_16),stringsAsFactors = F)
data_10=data_10[6:nrow(data_10),]
data_16=data_16[6:nrow(data_16),]
all_id=intersect(row.names(data_10),row.names(data_16))
all_id=all_id[!is.na(acc[all_id,"lat"])]
all_id=intersect(all_id,ft1016$id)
plot(fts_sum[all_id,"lat"])
data_10=data_10[all_id,]
data_16=data_16[all_id,]
data_10=data.frame(apply(data_10,2,function(x){as.numeric(x)}),stringsAsFactors = F)
data_16=data.frame(apply(data_16,2,function(x){as.numeric(x)}),stringsAsFactors = F)
row.names(data_10)=all_id
row.names(data_16)=all_id
a1=data_10
a2=data_16

data_10=data.frame(t(gene_CHG_10),stringsAsFactors = F)
data_16=data.frame(t(gene_CHG_16),stringsAsFactors = F)
data_10=data_10[6:nrow(data_10),]
data_16=data_16[6:nrow(data_16),]
all_id=intersect(row.names(data_10),row.names(data_16))
all_id=all_id[!is.na(acc[all_id,"lat"])]
all_id=intersect(all_id,ft1016$id)
plot(fts_sum[all_id,"lat"])
data_10=data_10[all_id,]
data_16=data_16[all_id,]
data_10=data.frame(apply(data_10,2,function(x){as.numeric(x)}),stringsAsFactors = F)
data_16=data.frame(apply(data_16,2,function(x){as.numeric(x)}),stringsAsFactors = F)
row.names(data_10)=all_id
row.names(data_16)=all_id
b1=data_10
b2=data_16

cols=ifelse(acc[all_id,"lat"]>60,"blue","red")
t.test((b1[,flc]-b2[,flc])~cols)
plot(acc[all_id,"lat"],b1[,flc]-b2[,flc])
abline(v=60)
ids=intersect(all_id,fts_sum$id)
points(fts_sum[ids,"lat"],b1[ids,flc]-b2[ids,flc],col="red")
t.test(b1[ids,flc]-b2[ids,flc]~fts_sum[ids,"lat"]>60)


data_10=data.frame(t(gene_CHH_10),stringsAsFactors = F)
data_16=data.frame(t(gene_CHH_16),stringsAsFactors = F)
data_10=data_10[6:nrow(data_10),]
data_16=data_16[6:nrow(data_16),]
all_id=intersect(row.names(data_10),row.names(data_16))
all_id=all_id[!is.na(acc[all_id,"lat"])]
all_id=intersect(all_id,ft1016$id)
plot(fts_sum[all_id,"lat"])
data_10=data_10[all_id,]
data_16=data_16[all_id,]
data_10=data.frame(apply(data_10,2,function(x){as.numeric(x)}),stringsAsFactors = F)
data_16=data.frame(apply(data_16,2,function(x){as.numeric(x)}),stringsAsFactors = F)
row.names(data_10)=all_id
row.names(data_16)=all_id
c1=data_10
c2=data_16

nrow(a1)
nrow(b1)
nrow(c1)


# 6A ----------------------------------------------------------------------
pdf("all figs/fig6/6A.pdf",width = 10,height=6)
par(mar=c(3,5,3,1))

cols=ifelse(acc[all_id,"lat"]>60,"blue","red")
boxplot(a1[,flc],a2[,flc],b1[,flc],b2[,flc],c1[,flc],c2[,flc],col=c("skyblue","pink","skyblue","pink","skyblue","pink"),
        main="",names=c("10C CG","16C CG","10C CHG","16C CHG","10C CHH","16C CHH"),ylab="Average methylation level",frame.plot=F,
        ylim=c(0,0.08))
for(i in 1:length(cols)){
  segments(1,a1[,flc][i], 2, a2[,flc][i], col = "gray", lwd = 1,lty="dashed")
  segments(3,b1[,flc][i], 4, b2[,flc][i], col = "gray", lwd = 1,lty="dashed")
  segments(5,c1[,flc][i], 6, c2[,flc][i], col = "gray", lwd = 1,lty="dashed")
}
points(x=rep(1,length(cols)),y=a1[,flc],pch=20,cex=1,col=cols)
points(x=rep(2,length(cols)),y=a2[,flc],pch=20,cex=1,col=cols)
points(x=rep(3,length(cols)),y=b1[,flc],pch=20,cex=1,col=cols)
points(x=rep(4,length(cols)),y=b2[,flc],pch=20,cex=1,col=cols)
points(x=rep(5,length(cols)),y=c1[,flc],pch=20,cex=1,col=cols)
points(x=rep(6,length(cols)),y=c2[,flc],pch=20,cex=1,col=cols)
dev.off()
t.test(a1[,flc]~cols)
t.test(a2[,flc]~cols)
t.test(b1[,flc]~cols)
t.test(b2[,flc]~cols)
t.test(c1[,flc]~cols)
t.test(c2[,flc]~cols)

# no significant ----------------------------------------------------------------------


t.test(a1[,flc]-a2[,flc]~cols)
t.test(b1[,flc]-b2[,flc]~cols)
t.test(c1[,flc]-c2[,flc]~cols)

wilcox.test(a1[,flc]-a2[,flc]~cols)
wilcox.test(b1[,flc]-b2[,flc]~cols)
wilcox.test(c1[,flc]-c2[,flc]~cols)

t.test(abs(a1[,flc]-a2[,flc])~cols)
t.test(abs(b1[,flc]-b2[,flc])~cols)
t.test(abs(c1[,flc]-c2[,flc])~cols)

pdf("all figs/fig6/6B.pdf",width = 5,height=5)
par(mar=c(3,5,3,1))
boxplot(b1[,flc]-b2[,flc]~cols,names=c("Northern","Southern"),
        col=c("#00BFC4","tomato"),xlab="",frame.plot=F,ylab="CHG change")
dev.off()


# 6B ----------------------------------------------------------------------
pdf("all figs/fig6/6B.pdf",width = 10,height=5)
par(mar=c(3,5,3,1))
plot(b1[,flc]-b2[,flc],p10[all_id,flc]-p16[all_id,flc],frame.plot=F,
     col="#00BFC4",xlab="FLC sensitivity",ylab="CHG change",pch=20,cex=1)
abline(lm(p10[all_id,flc]-p16[all_id,flc]~b1[,flc]-b2[,flc]),col="#00BFC4")
dev.off()
cor.test(a1[,flc]-a2[,flc],p10[all_id,flc]-p16[all_id,flc],method="pearson")
cor.test(b1[,flc]-b2[,flc],p10[all_id,flc]-p16[all_id,flc],method="pearson")
cor.test(b1[,flc]-b2[,flc],p10[all_id,flc]-p16[all_id,flc],method="spearman")
cor.test(c1[,flc]-c2[,flc],p10[all_id,flc]-p16[all_id,flc],method="pearson")


# 6C ----------------------------------------------------------------------


pdf("all figs/fig6/6C.pdf",width = 8,height=8)
nor=na.omit(all_id[acc[all_id,"lat"]>60])
sou=na.omit(all_id[acc[all_id,"lat"]<60])
par(mfrow=c(3,1))
par(mar=c(1,5,1,1))
type="CHH"
{
  FLC_M_10_n=FLC_M_10[,str_split_fixed(names(FLC_M_10),"_",2)[,1]%in%nor]
  FLC_M_10_s=FLC_M_10[,str_split_fixed(names(FLC_M_10),"_",2)[,1]%in%sou]
  FLC_M_16_n=FLC_M_16[,str_split_fixed(names(FLC_M_16),"_",2)[,1]%in%nor]
  FLC_M_16_s=FLC_M_16[,str_split_fixed(names(FLC_M_16),"_",2)[,1]%in%sou]
  mean_M_10_n=data.frame(FLC_M_10_n[,1],FLC_M_10_n[,2],apply(FLC_M_10_n[,seq(3,ncol(FLC_M_10_n),3)],1,function(x){mean(x,na.rm=T)}))
  mean_M_16_n=data.frame(FLC_M_16_n[,1],FLC_M_16_n[,2],apply(FLC_M_16_n[,seq(3,ncol(FLC_M_16_n),3)],1,function(x){mean(x,na.rm=T)}))
  mean_M_10_s=data.frame(FLC_M_10_s[,1],FLC_M_10_s[,2],apply(FLC_M_10_s[,seq(3,ncol(FLC_M_10_s),3)],1,function(x){mean(x,na.rm=T)}))
  mean_M_16_s=data.frame(FLC_M_16_s[,1],FLC_M_16_s[,2],apply(FLC_M_16_s[,seq(3,ncol(FLC_M_16_s),3)],1,function(x){mean(x,na.rm=T)}))
  names(mean_M_10_n)=c("strand","type","ratio")
  names(mean_M_16_n)=c("strand","type","ratio")
  names(mean_M_10_s)=c("strand","type","ratio")
  names(mean_M_16_s)=c("strand","type","ratio")
  info1=!is.nan(mean_M_10_n$ratio)&!is.na(mean_M_10_n$strand)&!is.na(mean_M_10_n$type)
  info2=!is.nan(mean_M_16_n$ratio)&!is.na(mean_M_16_n$strand)&!is.na(mean_M_16_n$type)
  table(info1)
  table(info2)
  info=(info1&info2)
  mean_M_10_n=mean_M_10_n[info,]
  mean_M_16_n=mean_M_16_n[info,]
  mean_M_10_s=mean_M_10_s[info,]
  mean_M_16_s=mean_M_16_s[info,]
  if(type=="CG"){
    lim=c(0,0.7*4)
    mean_M_10_n=mean_M_10_n[mean_M_10_n$type=="CG",]
    mean_M_16_n=mean_M_16_n[mean_M_16_n$type=="CG",]
    mean_M_10_s=mean_M_10_s[mean_M_10_s$type=="CG",]
    mean_M_16_s=mean_M_16_s[mean_M_16_s$type=="CG",]
  }else if(type=="CHG"){
    lim=c(0,0.3*4)
    mean_M_10_n=mean_M_10_n[mean_M_10_n$type=="CHG",]
    mean_M_16_n=mean_M_16_n[mean_M_16_n$type=="CHG",]
    mean_M_10_s=mean_M_10_s[mean_M_10_s$type=="CHG",]
    mean_M_16_s=mean_M_16_s[mean_M_16_s$type=="CHG",]
  }else if(type=="CHH"){
    lim=c(0,0.3*4)
    mean_M_10_n=mean_M_10_n[mean_M_10_n$type=="CHH",]
    mean_M_16_n=mean_M_16_n[mean_M_16_n$type=="CHH",]
    mean_M_10_s=mean_M_10_s[mean_M_10_s$type=="CHH",]
    mean_M_16_s=mean_M_16_s[mean_M_16_s$type=="CHH",]
  }
  mean_M_10_n=mean_M_10_n[mean_M_10_n$strand=="-",]
  mean_M_16_n=mean_M_16_n[mean_M_16_n$strand=="-",]
  mean_M_10_s=mean_M_10_s[mean_M_10_s$strand=="-",]
  mean_M_16_s=mean_M_16_s[mean_M_16_s$strand=="-",]
  
  mean_M_10_n=na.omit(mean_M_10_n)
  mean_M_16_n=na.omit(mean_M_16_n)
  mean_M_10_s=na.omit(mean_M_10_s)
  mean_M_16_s=na.omit(mean_M_16_s)
  
  
  
  plot(x=3170328:3182448,y=rep(1,length(3170328:3182448)), type="n",ylim=lim,main="",
       frame.plot=F,xaxt="n",yaxt="n",ylab="",xlab="")
  data=smooth.spline(x=row.names(mean_M_10_s),y=mean_M_10_s$ratio)
  data$y[data$y<0]=0
  data$y[1]=(max(lim)/4)*0;data$y[length(data$y)]=(max(lim)/4)*0
  polygon(data,col="#00BFC4",border = F)
  data=smooth.spline(x=row.names(mean_M_16_s),y=mean_M_16_s$ratio)
  data$y[data$y<0]=0;data$y=data$y+max(lim)/4
  data$y[1]=(max(lim)/4)*1;data$y[length(data$y)]=(max(lim)/4)*1
  polygon(data,col="tomato",border = F)
  data=smooth.spline(x=row.names(mean_M_10_n),y=mean_M_10_n$ratio)
  data$y=abs(data$y);data$y=data$y+(max(lim)/4)*2
  data$y[1]=(max(lim)/4)*2;data$y[length(data$y)]=(max(lim)/4)*2
  polygon(data,col="#00BFC4",border = F)
  data=smooth.spline(x=row.names(mean_M_16_n),y=mean_M_16_n$ratio)
  data$y[data$y<0]=0;data$y=data$y+(max(lim)/4)*3
  data$y[1]=(max(lim)/4)*3;data$y[length(data$y)]=(max(lim)/4)*3
  polygon(data,col="tomato",border = F)

  axis(2,c(seq(0,max(lim),max(lim)/4)+(max(lim)/4)/2)[1:4],c("Southern 10C","Southern 16C","Northern 10C","Northern 16C"),las=1,cex.axis=1)
}
dev.off()



gene=read.table(file="E:elife_data\\Methylation\\gene_position.txt",header = T,stringsAsFactors = F)
row.names(gene)=gene$Probe
gene$chr=paste0("chr",gene$chr)
gene$start=gene$start-3000
gene$end=gene$end+3000
gene=gene[flc,]
gene$start=gene$end-3000
library(methylKit)
region=GRanges(seqnames=gene$chr,ranges=IRanges(start=gene$start,end=gene$end),strand=gene$strand)
all_gene=paste0(gene$chr,"_",gene$start,"_",gene$strand)

files= list.files("E:elife_data\\Methylation\\Meth_types\\CHG\\")
path="E:elife_data\\Methylation\\Meth_types\\CHG\\"
gene_CHG_10=gene
gene_CHG_16=gene
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  #10s a data
  name=str_split(files[i],"_")[[1]]
  if(length(name)==3){
    if(substring(name[3],1,3)=="10C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHG=regionCounts(object=obj,regions=region)
      subname=paste0(CHG$chr,"_",CHG$start,"_",CHG$strand)
      row.names(CHG)=subname
      CHG=getData(CHG)
      gene_CHG_10[,name]=CHG[all_gene,"numCs"]/CHG[all_gene,"coverage"]
    }else if(substring(name[3],1,3)=="16C"){
      name=name[2]
      obj=methRead(paste0(path,files[i]),sample.id="test",assembly="hg18",header=F, context="CpH", resolution="base", mincov = 1,
                   pipeline=list(fraction=TRUE,chr.col=1,start.col=2,end.col=2,
                                 coverage.col=6,strand.col=3,freqC.col=5 )
      )  
      CHG=regionCounts(object=obj,regions=region)
      subname=paste0(CHG$chr,"_",CHG$start,"_",CHG$strand)
      row.names(CHG)=subname
      CHG=getData(CHG)
      gene_CHG_16[,name]=CHG[all_gene,"numCs"]/CHG[all_gene,"coverage"]
    }
  }
}

gene_CHG_10_flc=gene_CHG_10
gene_CHG_16_flc=gene_CHG_16
data_10=data.frame(t.data.frame(gene_CHG_10_flc),stringsAsFactors = F)
data_16=data.frame(t.data.frame(gene_CHG_16_flc),stringsAsFactors = F)

data_16=data.frame(data_16[6:nrow(data_16),])
all_id=intersect(row.names(data_10),row.names(data_16))
all_id=all_id[!is.na(acc[all_id,"lat"])]
all_id=intersect(all_id,ft1016$id)
plot(fts_sum[all_id,"lat"])
data_10=data_10[all_id,]
data_16=data_16[all_id,]
data_10=data.frame(apply(data_10,2,function(x){as.numeric(x)}),stringsAsFactors = F)
data_16=data.frame(apply(data_16,2,function(x){as.numeric(x)}),stringsAsFactors = F)
row.names(data_10)=all_id
row.names(data_16)=all_id
b1=data_10
b2=data_16

t.test(as.numeric(data_10[all_id,])-as.numeric(data_16[all_id,])~acc[all_id,"lat"]>60)
plot(acc[all_id,"lat"],as.numeric(data_10[all_id,])-as.numeric(data_16[all_id,]))
cor.test(as.numeric(data_10[all_id,])-as.numeric(data_16[all_id,]),p10[all_id,flc]-p16[all_id,flc])
