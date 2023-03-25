rm(list = ls()) ; graphics.off()
setwd("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/result/inferSNF_res_EachSNFfixed_V2/TCGA/RawNum100_CorThresh0.75_RF_scale_fold0/RTK_specificGene_SNF1_4")
library(stringi)
library(stringr)
library(ggpubr)
library(export)
load("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/27.Fig重绘/rawdata/TCGA/TCGA_hg38_FPKM_Mut_hg19ascatCNV.Rdata")
SNF.file = "RawNum100_CorThresh0.75_RF_scale_fold0_TCGA"
SNFnew = openxlsx::read.xlsx(paste0("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/rawdata/inferSNF_res_EachSNFfixed_V2/",SNF.file,".xlsx"),1,
                             rowNames = T)
SNF_Cluster = SNFnew$inferSNF ; names(SNF_Cluster) = rownames(SNFnew)


gene = c("EGFR","ALK", "KIT", "PDGFRA","MET")
for (i in gene){
  #  RNA 
  exp.fpkm.TT.log = log2(exp_fpkm_TCGA+1)
  EGFR = as.data.frame(t(exp.fpkm.TT.log[i,]))
  colnames(EGFR) = "gene"
  EGFR$SNF = SNF_Cluster[colnames(exp.fpkm.TT.log)]
  EGFR$SNF[EGFR$SNF != "SNF4"] = "Others"
  EGFR$SNF = factor(EGFR$SNF,levels = c("SNF4","Others"))
  #EGFR$SNF = factor(paste0("SNF",EGFR$SNF),levels = paste0("SNF",c(1:4)))
  
  shapiro.test(EGFR$gene)
  bartlett.test(gene~SNF,data=EGFR)
  
  # p = ggviolin(EGFR,x = "SNF",y= "gene",fill = "SNF",
  #              add = "boxplot", add.params = list(fill="white"),
  #              palette = c(color[4],"grey") )+
  #   #ylim(0,) +
  #   stat_compare_means(label.y = max(EGFR$gene)+1) +
  #   xlab(i) + 
  p = ggboxplot(EGFR,x = "SNF",y= "gene",color = "SNF",
               add = "jitter", 
               palette = c(color[4],"grey") )+
    #ylim(0,) +
    stat_compare_means(label.y = max(EGFR$gene)+1) +
    xlab(i) 
  ggsave(p,filename = paste0("./TCGA_SNF4_Others_",i,"_FPKMlog.pdf"), height = 4, width = 4)
  graph2ppt(p,file = paste0("./TCGA_SNF4_Others_",i,"_FPKMlog.ppt"), height = 4, width = 4)
}

i = "ALK"
EGFR = as.data.frame(t(exp.fpkm.TT.log[i,]))
colnames(EGFR) = "gene"
EGFR$SNF = SNF_Cluster[colnames(exp.fpkm.TT.log)]
EGFR$SNF[EGFR$SNF != "SNF4"] = "Others"
EGFR$SNF = factor(EGFR$SNF,levels = c("SNF4","Others"))
#EGFR$SNF = factor(paste0("SNF",EGFR$SNF),levels = paste0("SNF",c(1:4)))

shapiro.test(EGFR$gene)
bartlett.test(gene~SNF,data=EGFR)

# p = ggviolin(EGFR,x = "SNF",y= "gene",fill = "SNF",
#              add = "boxplot", add.params = list(fill="white"),
#              palette = c(color[4],"grey") )+
#   ylim(0,1) +
#   stat_compare_means(label.y = 1) +
#   xlab(i) + 
p = ggboxplot(EGFR,x = "SNF",y= "gene",color = "SNF",
             add = "jitter",
             palette = c(color[4],"grey") )+ ylim(0,1) +
  stat_compare_means(label.y = 1) +
  xlab(i) 
ggsave(p,filename = paste0("./TCGA_SNF4_Others_",i,"_FPKMlog.pdf"), height = 4, width = 4)
graph2ppt(p,file = paste0("./TCGA_SNF4_Others_",i,"_FPKMlog.ppt"), height = 4, width = 4)
