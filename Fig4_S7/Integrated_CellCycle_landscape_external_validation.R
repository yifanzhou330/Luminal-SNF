rm(list = ls()) ; graphics.off()

library(RColorBrewer)
library(pheatmap)
library(export)
CDK = c("CCND1","CDK1","CDK2")
load("/TCGA_hg38_FPKM_Mut_hg19ascatCNV.Rdata")
SNF.file = "RawNum100_CorThresh0.75_RF_scale_fold0_TCGA"
SNFnew = openxlsx::read.xlsx(paste0(SNF.file,".xlsx"),1,rowNames = T)
SNF_Cluster = SNFnew$inferSNF ; names(SNF_Cluster) = rownames(SNFnew)

#### RNA -------- 
CDK.hp = exp_fpkm_TCGA[CDK,names(SNF_Cluster)]
bk = unique(c(seq(-1,1,length = 100)))
p = pheatmap::pheatmap(CDK.hp[,names(sort(SNF_Cluster))],cluster_rows = F,cluster_cols = F,
                       scale = "row", show_colnames = F,
                       #annotation_col =annot_col, annotation_colors = annotation_color,
                       border_color = NA,
                       breaks = bk,
                       gaps_col = cumsum(as.matrix(table(SNF_Cluster))[,1])[-4],
                       gaps_row = 1:2,
                       cellheight = 15,cellwidth = 3,
                       filename = "FPKM_TCGAinferSNF_CDK_hp.pdf",
                       height = 10)
graph2ppt(p,file = "FPKM_TCGAinferSNF_CDK_hp.ppt",height = 10, width = 20)

## stat
CDK.hp = as.data.frame(t(CDK.hp))
CDK.hp$SNF = paste0("SNF",SNF_Cluster )

aov.res = as.data.frame(matrix(nrow = length(CDK), ncol = 2))
rownames(aov.res) = CDK ; colnames(aov.res) = c("p.val","adj.p")

for (i in CDK){
  mat = CDK.hp[,c(i,"SNF")]
  colnames(mat)[1] = "gene"
  shapiro.test(mat[,"gene"])
  bartlett.test(gene~SNF,data=mat)
  res = summary(aov(gene~SNF,data=mat))
  aov.res[i,"p.val"] = res[[1]][["Pr(>F)"]][1]
}
aov.res[,"adj.p"] = p.adjust(aov.res[,"p.val"],method = "fdr")
write.csv(aov.res, file =  "FPKM_TCGAinferSNF_CDK_hp_aov.csv")


#### CNV ------
CDK.hp = cnv_thre_hg19ascat_TCGA[CDK,names(SNF_Cluster)]
for (i in rownames(CDK.hp)){
  tmp = CDK.hp[i,]
  tmp[tmp < 0 ] = 0
  tmp[tmp > 0 ] = 1
  CDK.hp[i,] = tmp
}
annot_col = data.frame(row.names = names(SNF_Cluster), subtype =SNF_Cluster )
tmp = c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A")
names(tmp) = c(paste0("SNF",1:4))
annotation_color = c(list(subtype =tmp))

bk = unique(c(seq(-1,1,length = 50)))
p = pheatmap::pheatmap(CDK.hp[,names(sort(SNF_Cluster))],cluster_rows = F,cluster_cols = F,
                       show_colnames = F,
                       annotation_col =annot_col,
                       border_color = NA,
                       breaks = bk,
                       gaps_col = cumsum(as.matrix(table(SNF_Cluster))[,1])[-4],
                       gaps_row = 1:2,
                       color =  colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))(50),
                       annotation_colors = annotation_color,
                       cellheight = 15,cellwidth = 3,
                       filename = "CNV_SNF_CDK_hp.pdf",
                       height = 10)
graph2ppt(p,file = "CNV_SNF_CDK_hp.ppt",height = 10, width = 20)

## stat
CDK.hp = as.data.frame(t(CDK.hp))
CDK.hp$SNF = paste0("SNF",SNF_Cluster )

fisher.res = as.data.frame(matrix(nrow = length(CDK), ncol = 2))
rownames(fisher.res) = CDK ; colnames(fisher.res) = c("p.val","adj.p")

for (i in CDK){
  mat = CDK.hp[,c(i,"SNF")]
  colnames(mat)[1] = "gene"
  mytab = xtabs( ~ gene+SNF,data=mat)
  res = fisher.test( mytab)
  fisher.res[i,"p.val"] = res[["p.value"]]
}
fisher.res[,"adj.p"] = p.adjust(fisher.res[,"p.val"],method = "fdr")
write.csv(fisher.res, file =  "CNV_TCGAinferSNF_CDK_fisher.csv")



# MDM2 RNA -------
i= "MDM2"
exp.fpkm.TT.log = log2(exp_fpkm_TCGA+1)
gene = as.data.frame(t(exp.fpkm.TT.log[i,]))
colnames(gene ) = "gene"
gene$SNF = SNF_Cluster[colnames(exp.fpkm.TT.log)]
gene$SNF[gene$SNF != "SNF3"] = "Others"
gene$SNF = factor(gene$SNF, levels = c("SNF3","Others"))

shapiro.test(gene$gene)
bartlett.test(gene~SNF,data=gene)

p = ggboxplot(gene,x = "SNF",y= "gene",color = "SNF",
              add = "jitter", 
              palette = c(color[3], "grey") )+
  stat_compare_means(label.y = 3.5) 
ggsave(p, filename = paste0("./SNF3_",i,"_FPKMlog.pdf"),height = 4,width = 4)
graph2ppt(p,file = paste0("./SNF3_",i,"_FPKMlog.ppt"),height = 4,width = 4)

# ATM RNA --------
i= "ATM"
exp.fpkm.TT.log = log2(exp_fpkm_TCGA+1)
gene = as.data.frame(t(exp.fpkm.TT.log[i,]))
colnames(gene ) = "gene"
gene$SNF = SNF_Cluster[colnames(exp.fpkm.TT.log)]
gene$SNF[gene$SNF != "SNF3"] = "Others"
gene$SNF = factor(gene$SNF, levels = c("SNF3","Others"))

shapiro.test(gene$gene)
bartlett.test(gene~SNF,data=gene)

p = ggboxplot(gene,x = "SNF",y= "gene",color = "SNF",
              add = "jitter", 
              palette = c(color[3], "grey") )+
  stat_compare_means(label.y = 3) 
ggsave(p, filename = paste0("./SNF3_",i,"_FPKMlog.pdf"),height = 4,width = 4)
graph2ppt(p,file = paste0("./SNF3_",i,"_FPKMlog.ppt"),height = 4,width = 4)


# MDM2 AMP -----
i = "MDM2"
cnv = as.data.frame(t(cnv_thre_hg19ascat_TCGA[i,names(SNF_Cluster)]))
colnames(cnv) = "gene"
cnv$gene[cnv$gene < 0] = 0
cnv$gene[cnv$gene > 1 ] = 1
cnv$SNF = SNF_Cluster
cnv$SNF[cnv$SNF != "SNF3"] = "Others"
cnv$SNF = factor(cnv$SNF, levels = c("SNF3","Others"))

cluster.merge4plot = as.data.frame(table(cnv$SNF, cnv$gene))
colnames(cluster.merge4plot)[1:2] = c("SNF","gene")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$gene = factor(cluster.merge4plot2$gene,levels = c(0, 1)) 
cnv_color = c( "1" = "#E21F22","0" = "#eaeaea")

tmp = fisher.test(table(cnv$SNF, cnv$gene))
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = gene )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cnv_color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2), 
        axis.ticks = element_line(size = 2),
        axis.text.x=element_text(size = 20),
        axis.text.y=element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename = paste0("./SNF3_",i,"_4plot_percent_barplot.pdf"))
graph2ppt(p, file = paste0("./SNF3_",i,"_4plot_percent_barplot.ppt"))


# ATM Del -----
i = "ATM"
cnv = as.data.frame(t(cnv_thre_hg19ascat_TCGA[i,names(SNF_Cluster)]))
colnames(cnv) = "gene"
cnv$gene[cnv$gene > 0] = 0
cnv$gene[cnv$gene < 0 ] = 1
cnv$SNF = SNF_Cluster
cnv$SNF[cnv$SNF != "SNF3"] = "Others"
cnv$SNF = factor(cnv$SNF, levels = c("SNF3","Others"))

cluster.merge4plot = as.data.frame(table(cnv$SNF, cnv$gene))
colnames(cluster.merge4plot)[1:2] = c("SNF","gene")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$gene = factor(cluster.merge4plot2$gene,levels = c(0, 1)) 
#cnv_color = c("-2" = "#3654A4", "2" = "#E21F22","-1" = "#BBC3D9","1" ="#E4B5B8" ,"0" = "#eaeaea")
cnv_color = c( "1" = "#3654A4","0" = "#eaeaea")

tmp = fisher.test(table(cnv$SNF, cnv$gene))
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = gene )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cnv_color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2), 
        axis.ticks = element_line(size = 2),
        axis.text.x=element_text(size = 20),
        axis.text.y=element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename =  paste0("./SNF3_",i,"_4plot_percent_barplot.pdf"))
graph2ppt(p, file = paste0("./SNF3_",i,"_4plot_percent_barplot.ppt"))



