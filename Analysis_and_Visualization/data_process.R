# get TPM and genename ----------------------------------------------------
###get gene length
gene=gtf[gtf$V3=="exon",]
gene$name=sub("gene_id (\\S+) .+", "\\1", gsub(";", "", gene$V9))
gene$len=gene$V5-gene$V4
gene_len=aggregate(gene$len~gene$name,FUN=sum)
row.names(gene_len)=gene_len[,1]
names(gene_len)=c("gene","len")
###get exp
ids=str_split_fixed(files,"_",4)[,2]
temps=str_split_fixed(files,"_",4)[,3]
ids=str_replace(ids,"s","")
id_10=ids[temps=="10C"]
id_16=ids[temps=="16C"]
length(intersect(id_10,id_16))
ids=intersect(id_10,id_16)
table(ids%in%fl.data@phdata$id)
aa=fl.data@phdata[ids,]
aa=aa[!is.na(aa$id),]
aa
aa=aa[!is.na(aa$FT10_mean)&!is.na(aa$FT16_mean),]
aa=aa[!is.na(aa$FT10_sd)&!is.na(aa$FT16_sd),]
all_id=aa$id
aa
for(i in 1:length(files)){
  setTxtProgressBar(pb, i/length(files))
  data=data.frame(fread(paste0("E:\\elife_data\\Transcriptome\\GSE54680_RAW (1)\\",files[i])))
  data=data[,c(1,2)]
  names(data)=c("gene","exp")
  if(i==1){
    all_exp=data
    names(all_exp)[2]=paste0(str_split_fixed(files[i],"_",4)[,2:3],collapse = "_")
  }else{
    all_exp=merge.data.frame(all_exp,data,all=T,by="gene")
    names(all_exp)[i+1]=paste0(str_split_fixed(files[i],"_",4)[,2:3],collapse = "_")
  }
}
ids=str_split_fixed(names(all_exp),"_",4)[,1][-1]
temps=str_split_fixed(names(all_exp),"_",4)[,2][-1]
ids=str_replace(ids,"s","")
exp_10=all_exp[,c("gene",paste0("s",ids[temps=="10C"],"_10C"))]
exp_16=all_exp[,c("gene",paste0("s",ids[temps=="16C"],"_16C"))]
names(exp_10)[2:ncol(exp_10)]=ids[temps=="10C"]
names(exp_16)[2:ncol(exp_16)]=ids[temps=="16C"]
gene=as.character(exp_16$gene)
row.names(exp_10)=gene
row.names(exp_16)=gene
exp_10=exp_10[,-1]
exp_16=exp_16[,-1]
all_id[acc[all_id,"long"]<0]
all_id=all_id[acc[all_id,"long"]>0]
exp_10=exp_10[gene_len$gene,all_id]
exp_16=exp_16[gene_len$gene,all_id]
exp_10=exp_10[row.names(exp_10)%in%gene_len$gene,]
exp_16=exp_16[row.names(exp_16)%in%gene_len$gene,]



dim(exp_10)

countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

countToFpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

exp_10=data.frame(apply(exp_10,2,function(x)countToTpm(x,gene_len[row.names(exp_10),2])))
exp_16=data.frame(apply(exp_16,2,function(x)countToTpm(x,gene_len[row.names(exp_16),2])))
dim(exp_10)

exp_10=t.data.frame(exp_10)
exp_16=t.data.frame(exp_16)
row.names(exp_10)=all_id  
row.names(exp_16)=all_id

mag10=fl.data[all_id,]
mag16=fl.data[all_id,]

p10=cbind.data.frame(mag10@phdata[,1:2],exp_10)
p16=cbind.data.frame(mag16@phdata[,1:2],exp_16)

mag10@phdata=p10
mag16@phdata=p16

g10 <- colnames(p10)[-c(1,2)]
g16 <- colnames(p16)[-c(1,2)]
g1016 <- intersect(g10,g16)

na10 <- apply(p10[,-c(1,2)],2,function(x) sum(x>0,na.rm = T) )
na16 <- apply(p16[,-c(1,2)],2,function(x) sum(x>0,na.rm = T) )
if10 <- apply(p10[,-c(1,2)],2,function(x) sum(is.infinite(x),na.rm = T) )
if16 <- apply(p16[,-c(1,2)],2,function(x) sum(is.infinite(x),na.rm = T) )
na_t <- na10 + na16 > 107
na_if <- if10 + if10 < 50
g1016_qc <- g10[na_t & na_if]

ft1016=aa[all_id,]
