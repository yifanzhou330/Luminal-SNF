rm(list = ls())  ; graphics.off()
library(xlsx)
library(GSVA)
library(GSEABase)
library(RColorBrewer)
library(data.table)
library(stringi)  ; library(stringr)
library(limma)
library(ggpubr)
library(plyr)
library(export)
library(ComplexHeatmap) ; library(circlize)
test_summary = openxlsx:::read.xlsx("./RawNum100_CorThresh0.75_RF_scale_fold0_CPTAC.xlsx")
rownames(test_summary) = test_summary$X1


clinic = read.csv("./prosp-brca-v5.4-public-sample-annotation.csv",
                  row.names = 1)
clinic = clinic[rownames(test_summary),]
clinic$SNF =  test_summary$inferSNF
clinic$SNF = factor(clinic$SNF, levels = paste0("SNF",1:4))

rna = data.table::fread("./prosp-brca-v5.4-public-rnaseq-fpkm-log2.csv",
                        data.table = F)
rna = rna[!duplicated(rna$V1),]
rownames(rna) = rna[,1]
rna = rna[,-1]
rna = rna[,rownames(test_summary)]

cnv = read.csv("./prosp-brca-v5.4-public-gene-level-cnv-gistic2-all_data_by_genes.gct.csv")
cnv = cnv[!duplicated(cnv$geneSymbol),]
rownames(cnv) = cnv$geneSymbol
cnv = cnv[,-c(1:4)]
cnv = cnv[,rownames(test_summary)]

clinic$Stromal.Score = clinic$xCell.Stromal.Score
clinic$Immune.Score = clinic$xCell.Immune.Score

clinic$CINscore = NA
for (i in rownames(clinic)){
  tmp = cnv[,i]
  tmp = tmp^2
  clinic[i,"CINscore"] = sum(tmp)
}
clinic$CINscore = log2(clinic$CINscore)

gmt = getGmt("./GObp_Hall_Reactome_GOmf_v7_4.gmt") 
tmp = c("HALLMARK_ESTROGEN_RESPONSE_EARLY")#,"HALLMARK_ESTROGEN_RESPONSE_LATE")
gmt_picked = gmt[tmp]
path_sig = gsva(as.matrix(rna),gmt_picked,method = "gsva")#,kcdf = "Gaussian",abs.ranking=F)

#clinic = cbind(clinic,as.data.frame(t(path_sig[,rownames(clinic)])))
clinic = cbind(clinic,ESTROGEN_RESPONSE_EARLY = as.numeric(path_sig[,rownames(clinic)]) )
#########
clinic = clinic[,c("SNF","Stromal.Score", "CINscore","Immune.Score",
                   "ESTROGEN_RESPONSE_EARLY")]
colnames(clinic)[c(2,4)] = c("Stromal.Score", "Immune.Score") 
mat = aggregate(clinic[,-1],list(clinic$SNF),median)
rownames(mat) = mat$Group.1
mat = as.data.frame(t(mat[,-1]))


col_fun = colorRamp2(c(-0.6,-0.5, 0 , 0.5,0.6), rev(brewer.pal(n = 5, name ="RdBu")))
p=ComplexHeatmap::pheatmap(mat,scale = "row",cluster_rows = F,cluster_cols = F,
                           heatmap_legend_param = list(title = "Z-score"),
                           border_color = "white",
                           color = col_fun)
graph2pdf(p,file = "./CPTAC_molecular_feature.pdf",width = 5,height = 7)


######### kruskal.test
kruskal.test.res = as.data.frame(matrix(nrow = nrow(mat),ncol = 2))
rownames(kruskal.test.res) = rownames(mat) ; colnames(kruskal.test.res) = c("p.val",'adj.p')
for ( i in rownames(mat)){
  tmp = data.frame(gene = as.numeric(clinic[,i]),
                   group = clinic$SNF)
  res = kruskal.test(gene~group,data = tmp)
  kruskal.test.res[i,"p.val"] = res[["p.value"]]
  # res = summary(aov(gene~group,data = tmp))
  # aov.res[i,"p.val"] = res[[1]][["Pr(>F)"]][1]
}
kruskal.test.res$adj.p = p.adjust(kruskal.test.res$p.val,method = "fdr")

write.csv(kruskal.test.res,file = "./kruskal.test.res.csv")


