
# 3A ----------------------------------------------------------------------

##spearman co-expression network -----------------------------------------
library(parallel)
gene_all_10=g1016_qc
get_p=function(i){
  data=apply(p10_z[,gene_all_10], 2, function(x) cor.test(as.numeric(x), as.numeric(p10_z[,gene_all_10[i]]),method="spearman")$p.value)
  write.table(gene_all_10[i],file="test.txt",append = T,col.names = F,row.names = F,quote=F)
  return(data)
}
clnum<-15
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(get_p)))
clusterExport(cl,c("p10_z","gene_all_10"))
data=parLapply(cl,1:length(gene_all_10),get_p)
stopCluster(cl)
gene_cor_10p=as.data.frame(data,stringsAsFactors = F)
names(gene_cor_10p)=row.names(gene_cor_10p)
save(gene_cor_10p,file="final_data/gene_cor_10p.RData")
rm(data,gene_cor_10p)
gc()

gene_all_16=g1016_qc
get_p=function(i){
  data=apply(p16_z[,gene_all_16], 2, function(x) cor.test(as.numeric(x), as.numeric(p16_z[,gene_all_16[i]]),method="spearman")$p.value)
  write.table(gene_all_16[i],file="test.txt",append = T,col.names = F,row.names = F,quote=F)
  return(data)
}
clnum<-15
cl <- makeCluster(getOption("cl.cores", clnum));
clusterExport(cl,deparse(substitute(get_p)))
clusterExport(cl,c("p16_z","gene_all_16"))
data=parLapply(cl,1:length(gene_all_16),get_p)
stopCluster(cl)
gene_cor_16p=as.data.frame(data,stringsAsFactors = F)
names(gene_cor_16p)=row.names(gene_cor_16p)
save(gene_cor_16p,file="final_data/gene_cor_16p.RData")
rm(data,gene_cor_16p)
gc()

gene_cor_10b=cor(p10_z[,g1016_qc],method="spearman")
dim(gene_cor_10b)
names(gene_cor_10b)=g1016_qc
row.names(gene_cor_10b)=g1016_qc
save(gene_cor_10b,file="final_data/gene_cor_10b.RData")
rm(gene_cor_10b)
gene_cor_16b=cor(p16_z[,g1016_qc],method="spearman")
dim(gene_cor_16b)
names(gene_cor_16b)=g1016_qc
row.names(gene_cor_16b)=g1016_qc
save(gene_cor_16b,file="final_data/gene_cor_16b.RData")
rm(gene_cor_16b)


## summary network ---------------------------------------------------------
##get connection, mean coefficient, intersection
load("final_data/gene_cor_10p.RData")
load("final_data/gene_cor_10b.RData")
dim(gene_cor_10p)
dim(gene_cor_10b)
a=gene_cor_10p[gene_list,]
edge_10=data.frame(from="",to="",value=0)[-1,]
for(i in 1:nrow(a)){
  setTxtProgressBar(pb, i/nrow(a))
  data=a[i,a[i,]<0.05/19007**2]
  if(is.data.frame(data)){
    for(j in 1:ncol(data)){
      row=row.names(data)
      col=colnames(data)[j]
      edge_10=rbind.data.frame(edge_10,data.frame(from=row,to=col,value=gene_cor_10b[row,col]))
    }
  }
}
length(unique(edge_10$to))


load("final_data/gene_cor_16p.RData")
load("final_data/gene_cor_16b.RData")
dim(gene_cor_16p)
dim(gene_cor_16b)
a=gene_cor_16p[gene_list,]
edge_16=data.frame(from="",to="",value=0)[-1,]
for(i in 1:nrow(a)){
  setTxtProgressBar(pb, i/nrow(a))
  data=a[i,a[i,]<0.05/19007**2]
  if(is.data.frame(data)){
    for(j in 1:ncol(data)){
      row=row.names(data)
      col=colnames(data)[j]
      edge_16=rbind.data.frame(edge_16,data.frame(from=row,to=col,value=gene_cor_16b[row,col]))
    }
  }
}
length(unique(edge_16$to))
table(unique(edge_16$to)%in%peripheral)
table(unique(edge_10$to)%in%peripheral)

#save(gene_list,peripheral,g1016_qc,file="interdata.RData")

a=row.names(gene_cor_10p)
b=row.names(gene_cor_16p)
a=intersect(a,b)
#gent all gene connection
data=lapply(1:length(a),function(x){
  setTxtProgressBar(pb, x/length(a))
  p_10=as.numeric(gene_cor_10p[a[x],])
  p_10=ifelse(p_10<0.05/18436**2,T,F)
  x2=sum(p_10,na.rm = T)
  p_16=as.numeric(gene_cor_16p[a[x],])
  p_16=ifelse(p_16<0.05/18343**2,T,F)
  x3=sum(p_16,na.rm = T)
  x4=length(intersect(names(gene_cor_10p)[p_10],names(gene_cor_16p)[p_16]))
  return(c(a[x],x2,x3,x4))
})
degree_all_final=as.data.frame(t(data.frame(data)),stringsAsFactors = F)
row.names(degree_all_final)=degree_all_final$V1
names(degree_all_final)=c("gene","n10","n16","inter")
degree_all_final[,2]=as.numeric(degree_all_final[,2])
degree_all_final[,3]=as.numeric(degree_all_final[,3])
degree_all_final[,4]=as.numeric(degree_all_final[,4])

## WGCNA -------------------------------------------------------------------

### 10C ---------------------------------------------------------------------
library(WGCNA)
power1=6
cutheight=0.5
deepsplit=2
type = "unsigned"
corType = "pearson"
robustY=T
maxPOutliers=1
#change 
p10_qc=p10_z[,g1016_qc]
p10_qc=data.frame(t(p10_qc))
m.mad <- apply(p10_qc,1,mad)
dataExprVar <- p10_qc[which(m.mad > quantile(m.mad, probs=seq(0, 1, 0.25))[mad1]),]
dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.80,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power=16
#run time 1h
net_10= blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = type, minModuleSize = 20,
                         reassignThreshold = 0, mergeCutHeight = 0.2,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=T,  saveTOMFileBase="final_data\\10C_wgcna",
                         corType = corType, deepSplit = 1,
                         maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                         verbose = 3)
net_10_2=recutBlockwiseTrees(datExpr=dataExpr,TOMFiles = net_10$TOMFiles, 
                             goodGenes = net_10$goodGenes,goodSamples = net_10$goodSamples,
                             blocks = net_10$blocks,dendrograms=net_10$dendrograms,
                             minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.2,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             corType = corType, deepSplit =2,
                             maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                             verbose = 3)
table(net_10$color)
table(net_10_2$color)

### 16C ---------------------------------------------------------------------
p10_qc=p16_z[,g1016_qc]
p10_qc=data.frame(t(p10_qc))
m.mad <- apply(p10_qc,1,mad)
dataExprVar <- p10_qc[which(m.mad > quantile(m.mad, probs=seq(0, 1, 0.25))[mad1]),]
dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, networkType=type, verbose=5)

par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.80,col="red")

# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
power=14
net_16 = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                          TOMType = type, minModuleSize = 20,
                          reassignThreshold = 0, mergeCutHeight = 0.2,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs=T,  saveTOMFileBase="final_data\\16C_wgcna",
                          corType = corType, deepSplit = 2,
                          maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                          verbose = 3)
net_16_2=recutBlockwiseTrees(datExpr=dataExpr,TOMFiles = net_16$TOMFiles, 
                             goodGenes = net_16$goodGenes,goodSamples = net_16$goodSamples,
                             blocks = net_16$blocks,dendrograms=net_16$dendrograms,
                             minModuleSize = 30,
                             reassignThreshold = 0, mergeCutHeight = 0.2,
                             numericLabels = TRUE, pamRespectsDendro = FALSE,
                             corType = corType, deepSplit =2,
                             maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                             verbose = 3)

table(net_16$colors)
table(net_16_2$colors)

table(net_10_2$colors)
table(net_16_2$colors)
## plot --------------------------------------------------------------------
library(ggplot2)
library(enrichplot)
library(topGO)
library(data.table)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(parallel)
library(RColorBrewer)


### 10C ---------------------------------------------------------------------
edges=edge_10[edge_10$from%in%gene_list,]
edges$main=T
nodes=data.frame(gene=g1016_qc,stringsAsFactors = F)
graph=graph_from_data_frame(edges, vertices = nodes)
deg=degree(graph)
nodes$degree=deg
nodes$group=""
for(i in 1:nrow(nodes)){
  if(i==nrow(nodes)){
    break
  }
  if(nodes[i,1]%in%gene_list){
    nodes[i,3]=nodes[i,1]
  }else{
    if(length(edges[edges$to==nodes[i,1],1]>1)){
      nodes[i,3]=edges[edges$to==nodes[i,1],1][1]
    }else if(length(edges[edges$to==nodes[i,1],1]==1)){
      nodes[i,3]=edges[edges$to==nodes[i,1],1]
    }else{
      nodes[i,3]=NA
    }
  }
}
graph <- graph_from_data_frame(edges, vertices = nodes)
layout <- create_layout(graph, layout = 'circle')
outer=layout[!layout$name%in%gene_list,]
outer_group=data.frame()
for(i in 0:max(net_10_2$colors)){
  genes=names(net_10_2$colors[net_10_2$colors==i])
  outer_group=rbind.data.frame(outer_group,outer[outer$name%in%genes,])
}
genes=outer$name[!outer$name%in%names(net_10_2$colors)]
outer_group=rbind.data.frame(outer[outer$name%in%genes,],outer_group)
outer=outer_group
row.names(outer)=1:nrow(outer)
outer$x=cos((as.numeric(row.names(outer)) - 1) / nrow(outer) * 2 * pi)
outer$y=sin((as.numeric(row.names(outer)) - 1) / nrow(outer) * 2 * pi)

inner=layout[layout$name%in%gene_list,]
inner=bind_rows(inner[inner$group%in%setdiff(FT_gene,FTP_gene),],inner[inner$group%in%intersect(FT_gene,FTP_gene),],inner[inner$group%in%setdiff(FTP_gene,FT_gene),])
angles <- seq(1,360,360/(nrow(inner)))
radii <- rep(0.4,nrow(inner))
inner$x= radii * cos(angles / 180 * pi)
inner$y = radii * sin(angles / 180 * pi)
layout <- bind_rows(outer, inner) %>%
  arrange(.ggraph.index)
networkcolor=c("gray",ifelse(inner$group%in%FT_gene,"tomato","skyblue"))
networkcolor=c(rep("gray",length(genes)+table(net_10_2$colors)[1]),
               rep(sample(colorRampPalette(brewer.pal(length(unique(net_10_2$colors)), 'Set3'))(length(unique(net_10_2$colors))-1)),
                   times=as.numeric(table(net_10_2$colors)[2:length(unique(net_10_2$colors))])))
aa=net_10_2$colors
name=c()
for(i in 0:length(table(aa))){
  if(i==0){
    name=c(genes,names(aa[aa==0]))
  }else{
    name=c(name,names(aa[aa==i]))
  }
}
names(networkcolor)=name
networkcolor_10=networkcolor
plot(1:length(networkcolor),col=networkcolor)
networkcolor[FT_gene]="tomato"
networkcolor[FTP_gene]="skyblue"
networkcolor[intersect(FT_gene,FTP_gene)]="lightgreen"

g=ggraph(layout) +
  geom_edge_diagonal(
    aes(edge_color = node1.group, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE,alpha=0.2
  ) +
  geom_node_point(
    #size = degree
    aes(size = 0.8, color = name,filter = name %in% gene_list),
    alpha = 0.8, show.legend = FALSE,
    
  ) +
  geom_node_point(
    #size = degree
    aes(size = 0.1, color = name,filter = ! name %in% gene_list),
    alpha = 0.8, show.legend = FALSE,
    
  )+
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = (!name %in% gene_list)&(degree>1000000)
    ),
    size = 3, hjust = 'outward', family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = name,
      filter = name %in% gene_list&(degree>1000000)
    ),
    size = 4, hjust = 0.5, family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y - 0.045,
      label = ifelse(
        degree > 1000,
        format(degree, big.mark = ','),
        degree
      ),
      filter = (name %in% gene_list)&(degree>1000000)
    ),
    size = 4, hjust = 0.5, family = 'Oswald'
  ) +
  scale_edge_color_manual(values = networkcolor) +
  scale_color_manual(values = networkcolor) +
  scale_size_area(max_size = 20) +
  scale_edge_alpha_manual(values = c(0.5, 1)) +
  coord_fixed() +
  labs(
    title = ''
  ) +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
ggsave(g,filename = "all figs/fig3/3A1.png",width =20,height = 20)


### 16C ---------------------------------------------------------------------
edges=edge_16[edge_16$from%in%gene_list,]
edges$main=T
nodes=data.frame(gene=g1016_qc,stringsAsFactors = F)
graph=graph_from_data_frame(edges, vertices = nodes)
deg=degree(graph)
nodes$degree=deg
nodes$group=""
for(i in 1:nrow(nodes)){
  if(i==nrow(nodes)){
    break
  }
  if(nodes[i,1]%in%gene_list){
    nodes[i,3]=nodes[i,1]
  }else{
    if(length(edges[edges$to==nodes[i,1],1]>1)){
      nodes[i,3]=edges[edges$to==nodes[i,1],1][1]
    }else if(length(edges[edges$to==nodes[i,1],1]==1)){
      nodes[i,3]=edges[edges$to==nodes[i,1],1]
    }else{
      nodes[i,3]=NA
    }
  }
}
graph <- graph_from_data_frame(edges, vertices = nodes)
layout <- create_layout(graph, layout = 'circle')
outer=layout[!layout$name%in%gene_list,]
outer_group=data.frame()
for(i in 0:max(net_16_2$colors)){
  genes=names(net_16_2$colors[net_16_2$colors==i])
  outer_group=rbind.data.frame(outer_group,outer[outer$name%in%genes,])
}
genes=outer$name[!outer$name%in%names(net_16_2$colors)]
outer_group=rbind.data.frame(outer[outer$name%in%genes,],outer_group)
outer=outer_group
row.names(outer)=1:nrow(outer)
outer$x=cos((as.numeric(row.names(outer)) - 1) / nrow(outer) * 2 * pi)
outer$y=sin((as.numeric(row.names(outer)) - 1) / nrow(outer) * 2 * pi)

inner=layout[layout$name%in%gene_list,]
inner=bind_rows(inner[inner$group%in%setdiff(FT_gene,FTP_gene),],inner[inner$group%in%intersect(FT_gene,FTP_gene),],inner[inner$group%in%setdiff(FTP_gene,FT_gene),])
angles <- seq(1,360,360/(nrow(inner)))
radii <- rep(0.4,nrow(inner))
inner$x= radii * cos(angles / 180 * pi)
inner$y = radii * sin(angles / 180 * pi)
layout <- bind_rows(outer, inner) %>%
  arrange(.ggraph.index)
networkcolor=c(rep("gray",length(genes)+table(net_16_2$colors)[1]),
               rep(sample(colorRampPalette(brewer.pal(length(unique(net_16_2$colors)), 'Set3'))(length(table(net_16_2$colors))-1)),
                   times=as.numeric(table(net_16_2$colors)[2:length(unique(net_16_2$colors))])))
aa=net_16_2$colors
name=c()
for(i in 0:length(table(aa))){
  if(i==0){
    name=c(genes,names(aa[aa==0]))
  }else{
    name=c(name,names(aa[aa==i]))
  }
}
names(networkcolor)=name
networkcolor_16=networkcolor
plot(1:length(networkcolor),col=networkcolor)
networkcolor[FT_gene]="tomato"
networkcolor[FTP_gene]="skyblue"
networkcolor[intersect(FT_gene,FTP_gene)]="lightgreen"

g=ggraph(layout) +
  geom_edge_diagonal(
    aes(edge_color = node1.group, edge_alpha = as.factor(main)),
    edge_width = 0.2, show.legend = FALSE,alpha=0.2
  ) +
  geom_node_point(
    #size = degree
    aes(size = 0.8, color = name,filter = name %in% gene_list),
    alpha = 0.8, show.legend = FALSE,
    
  ) +
  geom_node_point(
    #size = degree
    aes(size = 0.1, color = name,filter = ! name %in% gene_list),
    alpha = 0.8, show.legend = FALSE,
    
  )+
  geom_node_text(
    aes(
      x = 1.0175 * x,
      y = 1.0175 * y,
      label = name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90,
      filter = (!name %in% gene_list)&(degree>10000)
    ),
    size = 3, hjust = 'outward', family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y,
      label = name,
      filter = name %in% gene_list&(degree>10000)
    ),
    size = 4, hjust = 0.5, family = 'Oswald'
  ) +
  geom_node_text(
    aes(
      x = x,
      y = y - 0.045,
      label = ifelse(
        degree > 1000,
        format(degree, big.mark = ','),
        degree
      ),
      filter = name %in% gene_list&(degree>10000)
    ),
    size = 4, hjust = 0.5, family = 'Oswald'
  ) +
  scale_edge_color_manual(values = networkcolor) +
  scale_color_manual(values = networkcolor) +
  scale_size_area(max_size = 20) +
  scale_edge_alpha_manual(values = c(0.5, 1)) +
  coord_fixed() +
  labs(
    title = ''
  ) +
  theme_void() +
  theme(
    text = element_text(family = 'Oswald'),
    legend.position = c(0.645, 0.51),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )
ggsave(g,filename = "all figs/fig3/3A2.png",width =20,height = 20)

# 3B ----------------------------------------------------------------------
png("all figs/fig3/3B.png",width=4,height=5,units="in",res=1000)
plot(c(5, 25), c(0, length(networkcolor_10)*1.1),col="white",frame.plot=F)
rect(9.5,1:length(networkcolor_10) - 0.5, 
     11, 1:length(networkcolor_10)+ 0.5,
     col =networkcolor_10, border = NA)

rect(19,1:length(networkcolor_16) - 0.5, 
     20.5, 1:length(networkcolor_16)+ 0.5,
     col =networkcolor_16, border = NA)

length_left=c(0,table(networkcolor_10)[unique(networkcolor_10)])
length_right=c(0,table(networkcolor_16)[unique(networkcolor_16)])
for(i in 2:length(length_left)){
  length_left[i]=sum(length_left[(i-1):i])
}
for(i in 2:length(length_right)){
  length_right[i]=sum(length_right[(i-1):i])
}
for(i in 1:length(unique(networkcolor_10))){
  genes=names(networkcolor_10)[networkcolor_10==unique(networkcolor_10)[i]]
  for(j in 1:length(unique(networkcolor_16))){
    genes2=names(networkcolor_16)[networkcolor_16==unique(networkcolor_16)[j]]
    if(length(intersect(genes2,genes))>1){
      
      x=c(10*1.1, 10*1.3, 10*1.5, 10*1.7, 10*1.9, 
          10*1.9, 10*1.7, 10*1.5, 10*1.3, 10*1.1, 10*1.1) 
      left.start=length_left[i]-0.5
      left.end=length_left[i]+length(intersect(genes2,genes))+0.5
      right.start=length_right[j]-0.5
      right.end=length_right[j]+length(intersect(genes2,genes))+0.5
      y=c(left.start, left.start, (left.start + right.start)/2, right.start,
          right.start, right.end, right.end, (left.end + right.end)/2, left.end, 
          left.end, left.start)
      s=c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0)
      xspline(x=x,y=y,shape=s,open=F,border = NA,
              col=rgb(col2rgb(unique(networkcolor_10)[i])[1, 1], col2rgb(unique(networkcolor_10)[i])[2, 1], col2rgb(unique(networkcolor_10)[i])[3, 1], 75, max = 255))  
    }
    length_right[j]=length_right[j]+length(intersect(genes2,genes))
    length_left[i]=length_left[i]+length(intersect(genes2,genes))
  }
}
for(i in 1:length(gene_list)){
  left.start=which(names(networkcolor_10)==gene_list[i])-2
  right.start=which(names(networkcolor_16)==gene_list[i])-2
  left.end=left.start+2
  right.end=right.start+2
  x=c(10*1.1, 10*1.3, 10*1.5, 10*1.7, 10*1.9, 
      10*1.9, 10*1.7, 10*1.5, 10*1.3, 10*1.1, 10*1.1)
  y=c(left.start, left.start, (left.start + right.start)/2, right.start,
      right.start, right.end, right.end, (left.end + right.end)/2, left.end, 
      left.end, left.start)
  s=c(0, 1, 0, 1, 0, 0, 1, 0, 1, 0)
  # if(gene_list[i]%in%intersect(FTP_gene,FT_gene)){
  #   xspline(x=x,y=y,shape=s,open=F,border = "lightgreen",
  #           col="lightgreen")    
  # }else if(gene_list[i]%in%setdiff(FTP_gene,FT_gene)){
  #   xspline(x=x,y=y,shape=s,open=F,border = "skyblue",
  #           col="skyblue")    
  # }else if(gene_list[i]%in%setdiff(FT_gene,FTP_gene)){
  #   xspline(x=x,y=y,shape=s,open=F,border = "tomato",
  #           col="tomato")    
  # }
}
dev.off()


# 3C ----------------------------------------------------------------------
a=aggregate(value ~ from, data = edge_10, FUN = mean)
b=aggregate(value ~ from, data = edge_10, FUN = length)
data=cbind.data.frame(a,b[,2])
a=aggregate(value ~ from, data = edge_16, FUN = mean)
b=aggregate(value ~ from, data = edge_16, FUN = length)
data2=cbind.data.frame(a,b[,2])
data=merge.data.frame(data,data2,by="from",all=T)
names(data)=c("gene","m10","n10","m16","n16")
degree_twa_all=data
row.names(degree_twa_all)=degree_twa_all$gene
data=degree_twa_all
data$inter=0
for(i in 1:nrow(data)){
  gene=as.character(data[i,1])
  set_10=as.character(edge_10[edge_10$from==gene,2])
  set_16=as.character(edge_16[edge_16$from==gene,2])
  data$inter[i]=length(intersect(set_10,set_16))
}
data=data[order(data$inter),]
a=setdiff(gene_list,unique(c(as.character(edge_10$from),as.character(edge_16$from))))
data=rbind(data.frame(gene=a,m10=0,n10=0,m16=0,n16=0,inter=0),data)
  
pdf(file="all figs/fig3/3C.pdf",width=10.5,height=6)
barplot(data$n10,names.arg="",ylim=c(-3000,5000),border = "black",
        col="lightblue",yaxt="n",ylab="Peripheral genes num of core gene", xaxt = "n")
par(new = TRUE)
barplot(data$inter,names.arg="",ylim=c(-3000,5000),border = "black",
        col="lightgreen",yaxt="n")
par(new = TRUE)
barplot(-data$n16,names.arg="",ylim=c(-3000,5000),border = "black",
        col="pink",yaxt="n")
par(new = TRUE)
barplot(-data$inter,names.arg="",ylim=c(-3000,5000),border = "black",
        col="lightgreen",yaxt="n")
axis(side = 2, at = seq(-3000, 5000, by = 1000), 
     labels = c(seq(3000,0, by = -1000),seq(1000,5000, by = 1000)), las = 1)
legend("topleft", c("10C", "16C","intersect"), fill = c("lightblue", "pink","lightgreen"),bty="n")
dev.off()


# 3D ----------------------------------------------------------------------

pdf(file="all figs/fig3/3D.pdf",width=10.5,height=6)
barplot(degree_twa_all[as.character(data$gene),"m10"],ylim=c(-1,1),
        col="lightblue",ylab="Mean correlation coeffcient at 10C")
dev.off()

# 3E ----------------------------------------------------------------------


pdf(file="all figs/fig3/3E.pdf",width=10.5,height=6)
par(mar=c(7,5,1,1))
barplot(degree_twa_all[as.character(data$gene),"m16"],ylim=c(-1,1),names.arg=data$gene,las=2,
        col="pink",ylab="Mean correlation coeffcient at 16C")
dev.off()


# 3F ----------------------------------------------------------------------
##estimate PVE
data=degree_twa_all
data[is.na(data)]=0
##filter genes which have at least one gene in each env
genes=as.character(data$gene[data$n10+data$n16>=2&data$n16>0&data$n10>0])
cor_all=data.frame(matrix(nrow=length(genes),ncol=16))
edge_10$from=as.character(edge_10$from)
edge_10$to=as.character(edge_10$to)
edge_16$from=as.character(edge_16$from)
edge_16$to=as.character(edge_16$to)
for(i in 1:length(genes)){
  setTxtProgressBar(pb, i/length(genes))
  re10=edge_10[edge_10$from==genes[i],2]
  re16=edge_16[edge_16$from==genes[i],2]
  reall=unique(c(re10,re16))
  reinter=intersect(re10,re16)
  reouter=setdiff(reall,reinter)
  row.names(cor_all)[i]=genes[i]
  
  #cor_all[i,1]=hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=g1016_qc,phe=ft1016[,3],idname="")$pves[1]
  cor_all[i,1]=hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,3],idname="")$pves[1]
  cor_all[i,2]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,3],idname="")$pves)
  #cor_all[i,3]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reinter,phe=ft1016[,3],idname="")$pves)
  cor_all[i,4]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reouter,phe=ft1016[,3],idname="")$pves)
  
  #cor_all[i,5]=hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=g1016_qc,phe=ft1016[,4],idname="")$pves[1]
  cor_all[i,5]=hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,4],idname="")$pves[1]
  cor_all[i,6]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,4],idname="")$pves)
  #cor_all[i,7]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reinter,phe=ft1016[,4],idname="")$pves)
  cor_all[i,8]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reouter,phe=ft1016[,4],idname="")$pves)
  
  #cor_all[i,9]=hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=g1016_qc,phe=ft1016[,3]-ft1016[,4],idname="")$pves[1]
  cor_all[i,9]=hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,4],idname="")$pves[1]
  cor_all[i,10]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
  #cor_all[i,11]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reinter,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
  cor_all[i,12]=sum(hglm_fixed(exp=p10_z,fixgene=genes[i],randomgene=reouter,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
  
  #cor_all[i,13]=hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=g1016_qc,phe=ft1016[,3]-ft1016[,4],idname="")$pves[1]
  cor_all[i,13]=hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,4],idname="")$pves[1]
  cor_all[i,14]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reall,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
  #cor_all[i,15]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reinter,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
  cor_all[i,16]=sum(hglm_fixed(exp=p16_z,fixgene=genes[i],randomgene=reouter,phe=ft1016[,3]-ft1016[,4],idname="")$pves)
}

names(cor_all)=c(paste(rep(c("FT 10","FT 16","FTP 10","FTP 16"),each=4),rep(c("single gene","all gene","inner gene","differ gene"),time=4)))
data=cor_all
names(data)
data=data[,c(1,2,4,5,6,8,9,10,12,13,14,16)]
names(data)
(mean(data$`FT 10 all gene`-data$`FT 10 single gene`)+mean(data$`FT 16 all gene`-data$`FT 16 single gene`))/2
(mean(data$`FTP 16 all gene`-data$`FTP 10 single gene`)+mean(data$`FTP 16 all gene`-data$`FTP 16 single gene`))/2

##PVE plot
cor_gg=data.frame(matrix(nrow=12*nrow(data),ncol=3))
a=melt(data)
cor_gg[,1]=rep(c("FT10","FT16","FTP10","FTP16"),each=length(genes)*3)
cor_gg[,2]=as.numeric(a$value)
cor_gg[,3]=rep(rep(c("Central gene","Central+connecter genes","Central+rewired genes"),each=length(genes)),time=4)
names(cor_gg)=c("trait","value","class")
cor_gg$class=factor(cor_gg$class,levels=c("Central gene","Central+connecter genes","Central+rewired genes"))

comp=list(c("Central gene","Central+connecter genes"),c("Central+connecter genes","Central+rewired genes"),c("Central gene","Central+rewired genes"))
g=ggplot(cor_gg,aes(x=trait,y=value,fill=class))+
  geom_boxplot(alpha=0.2,width=0.6,position=position_dodge(width=0.8),
               size=0.75,outlier.color = "white",outlier.size=0.01)+
  theme_bw()+ylim(c(0,1))+ theme(axis.line = element_line(colour = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank(),
                                 panel.border = element_blank(),
                                 panel.background = element_blank()) 
ggsave(g,filename = "all figs/fig3/3F.pdf",width =15,height =5)


t.test(data[,1],data[,2])
t.test(data[,1],data[,3])
t.test(data[,2],data[,3])
t.test(data[,4],data[,5])
t.test(data[,4],data[,6])
t.test(data[,5],data[,6])
t.test(data[,7],data[,8])
t.test(data[,7],data[,9])
t.test(data[,8],data[,9])
t.test(data[,10],data[,11])
t.test(data[,10],data[,12])
t.test(data[,11],data[,12])


# data_analysis -----------------------------------------------------------


a=paste0(edge_10$from,edge_10$to)
b=paste0(edge_16$from,edge_16$to)
length(unique(c(a,b)))
length(unique(c(edge_10$to,edge_16$to)))

length(unique(c(a,b)))-length(intersect(a,b))
16391/19583

mean(median(data[,2])-median(data[,1]),
median(data[,5])-median(data[,4]))
mean(median(data[,8])-median(data[,7]),
median(data[,11])-median(data[,10]))

