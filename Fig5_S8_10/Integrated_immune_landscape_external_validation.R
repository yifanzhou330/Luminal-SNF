rm(list = ls()) ; graphics.off()
library(RColorBrewer)
library(pheatmap)
library(export)
library(ggpubr)

load("TCGA_hg38_FPKM_Mut_hg19ascatCNV.Rdata")
SNF.file = "RawNum100_CorThresh0.75_RF_scale_fold0_TCGA"
SNFnew = openxlsx::read.xlsx(paste0(SNF.file,".xlsx"),1, rowNames = T)
SNF_Cluster = SNFnew$inferSNF ; names(SNF_Cluster) = rownames(SNFnew)

exp.fpkm.TT.log = log2(exp_fpkm_TCGA +1)

# PDCD1
PDCD1 = as.data.frame(t(exp.fpkm.TT.log["PDCD1",names(SNF_Cluster)]))
PDCD1$SNF = SNF_Cluster
PDCD1$SNF[PDCD1$SNF != "SNF2"] = "Others"
PDCD1$SNF = factor(PDCD1$SNF,levels = c("SNF2","Others"))

plot1 = ggviolin(data = PDCD1,x = "SNF", y = "PDCD1",fill = "SNF",
         add = "boxplot", add.params = list(fill="white"),
          palette = c(color[2], "grey") ) + stat_compare_means(label.y = 4)
graph2pdf(plot1, file = "PDCD1_SNF2_vs_Others.pdf",height = 4,width = 4)  


## hp
immune  = c("CD8A","GZMA","PRF1","IDO1")

tmp = c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A")
names(tmp) = paste0("SNF",1:4)
annotation_color = c(list(SNF_Cluster =tmp))

annot_col = as.data.frame(SNF_Cluster)

immune.hp = exp_fpkm_TCGA[immune,names(SNF_Cluster)]
bk = unique(c(seq(-1,1,length = 100)))
p = pheatmap::pheatmap(immune.hp[,names(sort(SNF_Cluster))],cluster_rows = F,cluster_cols = F,
                       scale = "row", show_colnames = F,
                       annotation_col =annot_col, annotation_colors = annotation_color,
                       border_color = NA,
                       breaks = bk,
                       gaps_col = cumsum(as.matrix(table(SNF_Cluster))[,1])[-4],
                       gaps_row = 1:3,
                       cellheight = 15,cellwidth = 3,
                       filename = "FPKM_SNF_immune_hp.pdf",
                       height = 10)
graph2ppt(p,file = "FPKM_SNF_immune_hp.ppt",height = 10, width = 25)

