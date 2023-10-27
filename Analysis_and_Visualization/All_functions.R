library(hglm)
library(data.table)

TWAs_fun=function(idname,phe,exp,name,filepath){
  exp_sub=exp
  orm_fun(IID=idname,filepath=filepath,name=name,exp=exp_sub,filetype="bin")
  phe_all=cbind.data.frame(idname,idname,phe)
  write.table(phe_all,file = paste(filepath,name,"_phe.txt",sep=""),sep = "\t",quote = F,col.names = F,row.names = F)
  cmd_orm=paste(osca,"--moment --befile",paste0(filepath,name,"_be"),
        "--reml-maxit 150 --pheno",paste(filepath,name,"_phe.txt",sep=""), 
        "--out",paste(filepath,name,"_moa",sep=""))
  system(cmd_orm)
  data=read.table(paste(filepath,name,"_moa.mlma",sep=""),stringsAsFactors = F,header = T)
  row.names(data)=data$Probe
  print("##########this method is based moment not moa#########")
  print("##########please use osca linux version to run moa#########")
  return(data)
}
  
onestep_osca_pve=function(idname,phe,exp,genes,name,filepath,nrandom=1,genes2=c(),method="merge"){
  project=paste0(getwd(),"/")
  if(nrandom==1){
    exp_sub=exp[,genes]
    orm_fun(IID=idname,filepath=filepath,name=name,exp=exp_sub,filetype="bin")
    phe_all=cbind.data.frame(idname,idname,phe)
    write.table(phe_all,file = paste(project,filepath,name,"_phe.txt",sep=""),sep = "\t",quote = F,col.names = F,row.names = F)
    cmd_orm=paste(osca,"--reml","--orm",paste0(project,filepath,name,"_orm"),
                  "--pheno",paste0(project,filepath,name,"_phe.txt"),
                  "--out",paste0(project,filepath,name,"_pve"))
    shell(cmd_orm)
    pve=as.numeric(fread(file=paste0(filepath,name,"_pve.hsq"))[4,2:3])
    return(pve)
  }else if(nrandom==2){
    exp_sub=exp[,genes]
    orm_fun(IID=idname,filepath=filepath,name=paste0(name,"_1"),exp=exp_sub,filetype="bin")
    exp_sub=exp[,genes2]
    orm_fun(IID=idname,filepath=filepath,name=paste0(name,"_2"),exp=exp_sub,filetype="bin")
    write.table(c(paste0(filepath,paste0(name,"_1_orm")),paste0(filepath,paste0(name,"_2_orm"))),
                file=paste0(filepath,name,"_files.txt"),row.names = F,col.names = F,quote=F)
    phe_all=cbind.data.frame(idname,idname,phe)
    
    write.table(phe_all,file = paste(filepath,name,"_phe.txt",sep=""),sep = "\t",quote = F,col.names = F,row.names = F)
    
    cmd_orm=paste(osca,"--reml","--merge-orm",paste0(filepath,name,"_files.txt"),
                  "--pheno",paste0(filepath,name,"_phe.txt"),
                  "--out",paste0(filepath,name,"_pve"))
    system(cmd_orm)
    pve=as.numeric(fread(file=paste0(filepath,name,"_pve.hsq"))[5,2:3])
    pve2=as.numeric(fread(file=paste0(filepath,name,"_pve.hsq"))[6,2:3])
    return(c(pve,pve2))
  }
}
onestep_pve=function(idname="",phe,exp,genes,name,filepath,nrandom=1,genes2=c(),method="interaction",outputdata="pve"){
  library(data.table)
  library(hglm)
  
  if(length(idname)<10){
    idname=ft1016[,1]
  }
  if(nrandom==1){
    Z=orm_fun(IID=idname,filepath=filepath,name=name,exp=exp[,genes])
    pve=hglm_exp_pve(phe,Z,ncol(Z),nr=nrandom,outputdata=outputdata)  
  }else{
    if(method=="interaction"){
      Z1=orm_fun(IID=idname,filepath,name,exp[,genes])
      Z2=orm_fun(IID=idname,filepath,name,exp[,genes2])
      Z3=Z1%*%Z2
      Z=cbind(Z1,Z2,Z3)
      randc=c(ncol(Z1),ncol(Z2),ncol(Z3))
      pve=hglm_exp_pve(phe,Z,randc,nr=3,outputdata=outputdata)
    }else if(method=="merge"){
      Z=orm_fun(IID=idname,filepath,name,exp[,unique(c(genes,genes2))])
      pve=hglm_exp_pve(phe,Z,ncol(Z),1,outputdata=outputdata)  
    }
    else if(method=="separate"){
      Z1=orm_fun(IID=idname,filepath,name,exp[,genes])
      Z2=orm_fun(IID=idname,filepath,name,exp[,genes2])
      Z=cbind(Z1,Z2)
      randc=c(ncol(Z1),ncol(Z2))
      pve=hglm_exp_pve(phe,Z,randc,nr=2,outputdata=outputdata) 
    }
  }
  return(pve)
}
c.z.hglm <- function(kin){
  relmat <- kin
  relmat[upper.tri(relmat)] <- t(relmat)[upper.tri(relmat)]#上三角变成下三角
  svd <- svd(relmat)
  Z <- svd$u %*% diag(sqrt(svd$d)) #左边的矩阵乘奇异值
  return(Z)
}
hglm_exp_pve=function(phe,Z,RandC,nr,outputdata="pve"){
  if(nr==1){
    na_phe=complete.cases(phe)
    na_z=complete.cases(Z)
    na_all=na_phe&na_z
    Z=c.z.hglm(Z[na_all,na_all])
    Z=as.matrix(Z)
    X=matrix(rep(1, length(phe[na_all])))
    hglm_exp=hglm(X=X,y=phe[na_all],Z=Z,calc.like=T)
    if(outputdata=="pve"){
      return(hglm_exp$varRanef/(hglm_exp$varRanef+hglm_exp$varFix))    
    }else if(outputdata=="randeff"){
      ranef=hglm_exp$ranef
      names(ranef)=ft1016[na_all,1]
      return(ranef)
    }else if(outputdata=="breedingvalue"){
      bred=c(Z%*%hglm_exp$ranef)
      names(bred)=ft1016[na_all,1]
      return(bred)
    }else if(outputdata=="likelihood"){
      return(hglm_exp$likelihood)
    }else if(outputdata=="model"){
      return(hglm_exp)
    }
  }else if(nr>=2){
    na_phe=complete.cases(phe)
    na_z=complete.cases(Z)
    na_all=na_phe&na_z
    Z_all=matrix()
    nums=c()
    for(i in 1:length(RandC)){
      Z_sub=Z[,(RandC[1]*(i-1)+1):(RandC[1]*i)]
      Z_sub=c.z.hglm(Z_sub[na_all,na_all])
      if(i==1){
        Z_all=Z_sub  
        nums=c(nums,ncol(Z_sub))
      }else{
        Z_all=cbind(Z_all,Z_sub)
        nums=c(nums,ncol(Z_sub))
      }
    }
    Z=as.matrix(Z_all)
    X=matrix(rep(1, length(phe[na_all])))
    hglm_exp=hglm(X=X,y=phe[na_all],Z=Z,RandC=nums,calc.like=T)
    if(outputdata=="pve"){
      return( hglm_exp$varRanef/(sum(hglm_exp$varRanef)+hglm_exp$varFix))
    }else if(outputdata=="randeff"){
      ranef=hglm_exp$ranef
      names(ranef)=rep(ft1016[na_all,1],nr)
      return(ranef)
    }else if(outputdata=="breedingvalue"){
      bred=c(Z%*%hglm_exp$ranef)
      names(bred)=ft1016[na_all,1]
      return(bred)
    }else if(outputdata=="likelihood"){
      return(hglm_exp$likelihood)
    }else if(outputdata=="model"){
      return(hglm_exp)
    }
  }
}
onestep_hglm_snp_pve=function(snp,gwa,phe,nrandom=1,snp2=0,ibs1=0,ibs2=0){
  if(nrandom==1){
    if(ibs1==0){
      ibs1=ibs(gwa[,snp],weight="freq")
    }
    na_phe <- complete.cases(phe)
    geno=as.double.gwaa.data(gwa[,snp])
    na_snp=complete.cases(geno)
    a=cbind.data.frame(na_phe,na_snp)
    na_all=a[,1]&a[,2]
    X=matrix(rep(1, length(phe[na_all])))
    Z_snp<- c.z.hglm(kin =ibs1[na_all,na_all])
    onlyIBS_snp_hglm <- hglm(X = X, y = phe[na_all], Z = Z_snp)
    h2_snp <- onlyIBS_snp_hglm$varRanef/(sum(onlyIBS_snp_hglm$varRanef)+ onlyIBS_snp_hglm$varFix)
  }
  if(nrandom==2){
    if(ibs1==0){
      ibs1=ibs(gwa[,snp],weight="freq")
    }
    if(ibs2==0){
      ibs2=ibs(gwa[,snp2],weight="freq")
    }
    na_phe <- complete.cases(phe)
    geno=as.double.gwaa.data(gwa[,snp])
    na_snp=complete.cases(geno)
    geno=as.double.gwaa.data(gwa[,snp2])
    na_snp2=complete.cases(geno)
    a=cbind.data.frame(na_phe,na_snp,na_snp2)
    na_all=a[,1]&a[,2]&a[,3]
    X=matrix(rep(1, length(phe[na_all])))
    Z_snp1<- c.z.hglm(kin =ibs1[na_all,na_all])
    Z_snp2<- c.z.hglm(kin =ibs2[na_all,na_all])
    onlyIBS_snp_hglm <- hglm(X = X, y = phe[na_all], Z = cbind(Z_snp1,Z_snp2),RandC = c(ncol(Z_snp1),ncol(Z_snp2)))
    h2_snp=onlyIBS_snp_hglm$varRanef/(sum(onlyIBS_snp_hglm$varRanef)+ onlyIBS_snp_hglm$varFix)
  }
  return(h2_snp)  
}
grm_fun=function(filepath,name,gwa,filetype=""){
  gcta <- "C:\\Users\\MSI\\OneDrive\\桌面\\搞科研\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\gcta_v1.94.0Beta_windows_x86_64"
  plink=plink <- "C:\\Users\\MSI\\OneDrive\\桌面\\搞科研\\拟南芥or柳树\\分析代码\\2022.07.28\\plink_win64_20220402\\plink"
  export.plink(gwa,filebasename = paste0(filepath,name))
  cmd_plink=paste(plink,"--tfile",paste0(filepath,name),"--make-bed --out",paste0(filepath,name))
  system(cmd_plink)
  if(filetype==""){
    cmd_k <- paste(gcta,"--bfile",paste0(filepath,name), "--make-grm-gz","--out",paste0(filepath,name))
    system(cmd_k)
    grmfile=data.frame(fread(paste0(filepath,name,".grm.gz")),stringsAsFactors = F)
    grm_data=matrix(ncol=max(grmfile[,1]),nrow=max(grmfile[,1]))
    for(i in 1:nrow(grmfile)){
      a=as.numeric(grmfile[i,1:4])
      grm_data[a[1],a[2]]=a[4] 
    }
    grm_data[upper.tri(grm_data)] <- t(grm_data)[upper.tri(grm_data)]
    return(grm_data)    
  }else if(filetype=="b"){
    cmd_k <- paste(gcta,"--bfile",paste0(filepath,name), "--make-grm","--out",paste0(filepath,name))
    system(cmd_k)
  }
}
orm_fun=function(IID,filepath="data/PVE/",name,exp,filetype=""){
  if(filetype==""){
    osca="D:\\OneDrive\\桌面\\搞科研\\软件\\osca_Win\\OSCA_x64.exe"
    exp=cbind.data.frame(IID,exp)
    write.table(exp,file=paste(filepath,name,"_exp.txt",sep=""),
                row.names = F,col.names = T,quote = F,sep="\t")
    cmd_befile=paste(osca,"--efile",paste(filepath,name,"_exp.txt",sep=""),
                     "--gene-expression","--make-bod","--no-fid",
                     "--out",paste0(filepath,name,"_be"))
    system(cmd_befile)
    cmd_orm=paste(osca,"--befile",paste0(filepath,name,"_be"),"--make-orm-gz",
                  "--out",paste0(filepath,name,"_orm"))
    system(cmd_orm)
    ormfile=data.frame(fread(paste0(filepath,name,"_orm",".orm.gz")),stringsAsFactors = F)
    orm_data=matrix(ncol=max(ormfile[,1]),nrow=max(ormfile[,1]))
    for(i in 1:nrow(ormfile)){
      a=as.numeric(ormfile[i,1:4])
      orm_data[a[1],a[2]]=a[4] 
    }
    orm_data[upper.tri(orm_data)] <- t(orm_data)[upper.tri(orm_data)]
    return(orm_data)  
  }else if(filetype=="bin"){
    osca="D:\\OneDrive\\桌面\\搞科研\\软件\\osca_Win\\OSCA_x64.exe"
    exp=cbind.data.frame(IID,exp)
    write.table(exp,file=paste(filepath,name,"_exp.txt",sep=""),
                row.names = F,col.names = T,quote = F,sep="\t")
    cmd_befile=paste(osca,"--efile",paste(filepath,name,"_exp.txt",sep=""),
                     "--gene-expression","--make-bod","--no-fid",
                     "--out",paste0(filepath,name,"_be"))
    system(cmd_befile)
    cmd_orm=paste(osca,"--befile",paste0(filepath,name,"_be"),"--make-orm",
                  "--out",paste0(filepath,name,"_orm"))
    system(cmd_orm)
  }
}

var.exp=function(ph,kin,exp,genes){
  na <- !is.na(ph)
  Z <- c.z.hglm(kin = kin[na,na])
  X = as.matrix(cbind(1, exp[na,genes]))
  withsnp.hglm <- hglm(X, y = ph[na], Z = Z)
  eff <- withsnp.hglm$fixef[-1]
  return(eff)
}

est_S_s<- function(x,y){
  y<- y/mean(y,na.rm=T)
  ct10 <- lm(y~scale(x))
  #plot(x,y)
  #abline(ct10)
  sct10 <- summary(ct10)
  return(c(sct10$coefficients[2,1],sct10$coefficients[2,4]))
}    
est_S_c <- function(x,y){
  y<- y/mean(y,na.rm=T)
  ct10 <- lm(y~ scale(x)+I(as.numeric(scale(x))^2))
  sct10 <- summary(ct10)
  if(nrow(sct10$coefficients)==3){
    return(c(sct10$coefficients[3,1],sct10$coefficients[3,4]))
  }else{
    return(c(NA,NA))
  }
}

hglm_fixed=function(exp,fixgene,randomgene,phe,idname="",filepath="data/PVE/",name="a"){
  if(idname==""){
    idname=ft1016[,1]
  }
  print(length(randomgene))
  X=cbind(matrix(rep(1, length(idname))),exp[,fixgene])
  Z=orm_fun(IID=idname,filepath,name,exp[,randomgene])
  na1=!is.na(phe)  
  na2=complete.cases(X)
  na_all=na1&na2
  
  Z=c.z.hglm(Z[na_all,na_all])
  Z=as.matrix(Z)
  X=as.matrix(X[na_all,])
  phe=phe[na_all]
  hglm=hglm(X=X,y=phe,Z=Z,calc.like=T)
  # pve_fixed= 1- var(phe- X%*%hglm$fixef)/var(phe)
  # pve_random= 1-var(phe-X%*%hglm$fixef-Z%*%hglm$ranef)/var(phe)-pve_fixed
  pve_fixed=var(X%*%hglm$fixef)/var(phe)
  pve_random=var(Z%*%hglm$ranef)/var(phe)
  bv_fixed=X%*%hglm$fixef
  bv_random=Z%*%hglm$ranef
  bv_all=X%*%hglm$fixef+Z%*%hglm$ranef
  names(bv_fixed)=idname[na_all]
  names(bv_random)=idname[na_all]
  names(bv_all)=idname[na_all]
  return(list("pves"=c(pve_fixed,pve_random),"bv_fixed"=bv_fixed,"bv_random"=bv_random,"bv_all"=bv_all,"AIC"=hglm$likelihood$cAIC,"hglm"=hglm))
}

lm_fixed=function(exp,fixgene,phe,idname=""){
  if(idname==""){
    idname=ft1016[,1]
  }
  X=cbind(matrix(rep(1, length(idname))),exp[,fixgene])
  na1=!is.na(phe)  
  na2=complete.cases(X)
  na_all=na1&na2
  X=as.matrix(X[na_all,])
  phe=phe[na_all]
  a=lm(phe~X)
  a=summary(a)
  var(X%*%a$coefficients[,1])/var(phe)
  pve=var(X%*%a$coefficients[,1])/var(phe)
  return(pve)
}

plot_pheat=function(data,thresh,datatype="p",datap){
  if(datatype=="b"){
    if(is.numeric(thresh)){
      data2=datap
      for(i in 1:nrow(data)){
        data[i,][data2[i,]>thresh]=0
      }
      #pheatmap(data,cluster_rows = F,cluster_cols = F,col=colorRampPalette(c("blue","white","red"))(100)) 
      pheatmap(data,col=colorRampPalette(c("blue","white","red"))(100)) 
    }else if(thresh=="fdr"){
      datap=apply(datap, 2, function(x)p.adjust(x,method = "fdr"))
      data2=datap
      for(i in 1:nrow(data)){
        data[i,][datap[i,]>0.05]=0
      }
      pheatmap(data,cluster_rows = F,cluster_cols = F,col=colorRampPalette(c("blue","white","red"))(100)) 
    }
  }else if(datatype=="p"){
    if(is.numeric(thresh)){
      data2=data
      for(i in 1:nrow(data)){
        data[i,][data2[i,]<=thresh]=-log10(data[i,][data2[i,]<=thresh])
        data[i,][data2[i,]>thresh]=0
      }
      pheatmap(data,cluster_rows = F,cluster_cols = F,col=colorRampPalette(c("white","orange","red"))(100)) 
    }else if(thresh=="fdr"){
      data=apply(data, 2, function(x)p.adjust(x,method = "fdr"))
      data2=data
      for(i in 1:nrow(data)){
        data[i,][data2[i,]<=0.05]=-log10(data[i,][data2[i,]<=0.05])
        data[i,][data2[i,]>0.05]=0
      }
      pheatmap(data,cluster_rows = F,cluster_cols = F,col=colorRampPalette(c("white","orange","red"))(100)) 
    }
  }
}

var.exp=function(ph,kin,exp,genes){
  na <- !is.na(ph)
  Z <- c.z.hglm(kin = kin[na,na])
  withsnp.hglm <- hglm(X = as.matrix(cbind(1, exp[na,genes])), y = ph[na], Z = Z)
  eff <- withsnp.hglm$fixef
  eff=eff[-1]
  return(eff)
}


add_row=function(data1,data2){
  if(class(data2)=="data.frame"){
    
  }else{
    data2=data.frame(t(data2))  
  }
  names(data2)=names(data1)
  a=rbind.data.frame(data1,data2)
  names(a)=names(data1)
  return(a)
}

snp_filter=function(snp,p,gwaa){
  snps=c()
  chr=as.numeric(gsub("chr","",str_split_fixed(snp,"_",2)[,1]) )
  window <- seq(0,max(GenABEL::map(gwaa)),5e04)
  pos_all=GenABEL::map(gwaa[,snp])
  for(i in unique(chr)){
    snp_in_chr=snp[chr==i]
    pos=pos_all[chr==i]
    inter_id=findInterval(pos,window)
    inter_num=unique(inter_id)
    for(j in inter_num){
      snp_in_win=snp_in_chr[inter_id==j][order(p[inter_id==j])==1]
      snps=c(snps,snp_in_win)
    }
  }
  return(snps)
}

gwa_fun=function(gwaa,trait,name){
  phe=trait
  hy <- polygenic(formula=phe,kinship.matrix =ibs,data = gwaa)
  mm=mmscore(hy,gwaa)
  par(mfrow=c(1,1))
  plot(mm,cex=0.5,main=paste(name,"  lambda =",round(mm@lambda$estimate,2)))
  abline(h=-log10(0.05/625070),col="red",lty="dashed")
}

genetic_cor=function(phe1,phe2,grm,outfile,name){
  gcta <- "C:\\Users\\MSI\\OneDrive\\桌面\\搞科研\\拟南芥or柳树\\分析代码\\2022.07.28\\gcta_v1.94.0Beta_windows_x86_64\\gcta_v1.94.0Beta_windows_x86_64\\bin\\gcta_v1.94.0Beta_windows_x86_64"
  data=data.frame(1:107,id=ft1016[,1],phe1=phe1,phe2=phe2)
  write.table(data,file = paste0(outfile,name,".txt"),row.names = F,col.names = F,quote = F)  
  cmd=paste(gcta,"--reml-bivar --grm",grm,"--pheno", paste0(outfile,name,".txt")," --reml-maxit 200 --out",paste0(outfile,name))
  system(cmd)
}

est_gcta_h2=function(idname,phe,gwadata,outfile,name){
  na=!is.na(phe)
  a=paste0(outfile,name,"_plink")
  export.plink(gwadata[na,],filebasename = a)
  cmd_plink=paste(plink,"--tfile",paste0(outfile,name,"_plink"),"--make-bed --out",paste0(outfile,name,"_plink"))
  shell(cmd_plink)
  cmd_k <- paste(gcta,"--bfile",paste0(outfile,name,"_plink"), "--make-grm" ,"--make-grm-alg 1","--out",paste0(outfile,name,"_grm"),sep=" ")
  shell(cmd_k)
  phe=data.frame("id"=idname,"family"=idname,"phe"=phe)
  phe=phe[na,]
  phe$id=1:nrow(phe)
  write.table(phe,file = paste0(outfile,name,"_phe.txt"),col.names = F,row.names = F,quote = F)
  cmd_biGREML <- paste(gcta,"--pheno",paste0(outfile,name,"_phe.txt"),"--grm",paste0(outfile,name,"_grm"),"--reml","--mpheno",1, "--out",paste0(outfile,name,"_h2"),sep=" ")
  shell(cmd_biGREML)
  data=fread(paste0(outfile,name,"_h2.hsq"))
  return(data[4,])
}
