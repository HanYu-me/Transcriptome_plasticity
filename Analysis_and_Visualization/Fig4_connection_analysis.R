peripheral=unique(c(as.character(unique(edge_16$to)),as.character(unique(edge_10$to))))
core=gene_list
unlink_gene=setdiff(g1016_qc,c(core,peripheral))

# 4A ----------------------------------------------------------------------

pdf("all figs/fig4/4A.pdf",width = 5,height=5)
par(mar=c(3,5,3,1))
boxplot(degree_all_gene[core,"n10"],degree_all_gene[peripheral,"n10"],degree_all_gene[unlink_gene,"n10"],main="",
        col=c("pink","skyblue","lightgreen"),frame.plot=F,names=c("Core gene","Peripheral gene","Unlinked gene"),
        ylab="Connectivity at 10C")
dev.off()
t.test(degree_all_gene[core,"n10"],degree_all_gene[peripheral,"n10"])$p.value
t.test(degree_all_gene[unlink_gene,"n10"],degree_all_gene[peripheral,"n10"])$p.value

# 4B ----------------------------------------------------------------------


pdf("all figs/fig4/4B.pdf",width = 5,height=5)
par(mar=c(3,5,3,1))
boxplot(degree_all_gene[core,"n16"],degree_all_gene[peripheral,"n16"],degree_all_gene[unlink_gene,"n16"],main="",
        col=c("pink","skyblue","lightgreen"),frame.plot=F,names=c("Core gene","Peripheral gene","Unlinked gene"),
        ylab="Connectivity at 16C")
dev.off()
t.test(degree_all_gene[core,"n16"],degree_all_gene[peripheral,"n16"])$p.value
t.test(degree_all_gene[unlink_gene,"n16"],degree_all_gene[peripheral,"n16"])


# 4C ----------------------------------------------------------------------


pdf("all figs/fig4/4C.pdf",width = 5,height=5)
par(mar=c(3,5,3,1))
degree_all_final$inter_rate=degree_all_final$inter/(degree_all_final$n10+degree_all_final$n16)
boxplot(degree_all_final[core,"inter_rate"],degree_all_final[peripheral,"inter_rate"],degree_all_final[unlink_gene,"inter_rate"],
        main="",col=c("pink","skyblue","lightgreen"),frame.plot=F,names=c("Core gene","Peripheral gene","Unlinked gene"),
        ylab="Connectivity conservation")
dev.off()


# 4D ----------------------------------------------------------------------
mean_10=apply(p10[,g1016_qc],2,function(x){mean(x,na.rm=T)})
mean_16=apply(p16[,g1016_qc],2,function(x){mean(x,na.rm=T)})
diff_1016=apply(p10[,g1016_qc]-p16[,g1016_qc],2,function(x){mean(abs(x),na.rm=T)})
names(mean_10)=g1016_qc
names(mean_16)=g1016_qc
names(diff_1016)=g1016_qc
rate_diff=diff_1016/(mean_10+mean_16)
pdf("all figs/fig4/4D.pdf",width = 5,height=5)
par(mar=c(3,5,3,1))
boxplot(rate_diff[core],rate_diff[peripheral],rate_diff[unlink_gene],main="",
        col=c("pink","skyblue","lightgreen"),frame.plot=F,names=c("Core gene","Peripheral gene","Unlinked gene"),
        ylab="Expression stability")
dev.off()


# 4E ----------------------------------------------------------------------
library(clusterProfiler)
library(org.At.tair.db)
go=enrichGO(peripheral, OrgDb = "org.At.tair.db", keyType="TAIR",ont="BP",qvalueCutoff = 0.05)
pdf("all figs/fig4/4E.pdf",width = 6,height=10)
par(mar=c(3,5,3,1))
barplot(go,showCategory=40,font.size=8)
dev.off()



# 4F ----------------------------------------------------------------------
go=enrichGO(unlink_gene, OrgDb = "org.At.tair.db", keyType="TAIR",ont="BP",qvalueCutoff = 0.05)
pdf("all figs/fig4/4F.pdf",width = 6,height=10)
par(mar=c(3,5,3,1))
barplot(go,showCategory=40,font.size=8)
dev.off()



# data analysis -----------------------------------------------------------
length(unlink_gene)
