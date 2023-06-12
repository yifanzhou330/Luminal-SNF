rm(list = ls()) ; graphics.off()
setwd("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/result/inferSNF_res_EachSNFfixed_V2/CPTAC_2022/RTK.hp.mean.SNF4phos")
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ComplexHeatmap)
library(export)
library(clusterProfiler)
library(Seurat)
library(tidyverse)
library(paletteer)
library(circlize)
library(RColorBrewer)
########################### PART 1. FUSCC SNF RTK RNA mean hp ---------
load("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/30.NG修改/rawdata/CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
gmt = read.gmt("/Users/ZYF/c5.go.mf.v7.4.symbols.gmt")
RTK = gmt[gmt$term == "GOMF_TRANSMEMBRANE_RECEPTOR_PROTEIN_KINASE_ACTIVITY",]
rownames(RTK) = RTK$gene
kinase = openxlsx::read.xlsx("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/30.NG修改/分型外推/CPTAC_2020/激酶kinase 参考文献：Driver Fusions and Their Implications in the Development and Treatment of Human Cancers.xlsx",1)
colnames(kinase) = kinase[1,] ; kinase = kinase[-1,] 
RTK = intersect(rownames(RTK),kinase$Kinase)

exp.fpkm.TT.log = log2(exp.fpkm.TT +1 )
RTK.hp = exp.fpkm.TT.log[RTK,]
RTK.hp = as.data.frame(t(RTK.hp))
RTK.hp$SNF = paste0("SNF",SNF_Cluster[rownames(RTK.hp)])
RTK.hp.mean = aggregate(RTK.hp[,1:(ncol(RTK.hp)-1)],list(RTK.hp$SNF),mean)

rownames(RTK.hp.mean) = RTK.hp.mean$Group.1 ; RTK.hp.mean = RTK.hp.mean[,-1]
RTK.hp.mean  = as.data.frame(t(RTK.hp.mean))

SNF4 = read.csv("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/30.NG修改/DEG_SNF/DEG.DESeq2.SNF4VsOthers.csv",
                row.names = 1)
SNF4 = na.omit(SNF4)
SNF4 = SNF4[intersect(rownames(RTK.hp.mean),rownames(SNF4)) , ]
SNF4 = SNF4[order(SNF4$log2FoldChange,decreasing = T),]

RTK.hp.mean = RTK.hp.mean[rownames(SNF4),]
ord = rownames(RTK.hp.mean)
RTK.hp.mean.FUSCCRNA.scale = t(scale(t(RTK.hp.mean)))

########################### PART 2. Lumianl neo scRNA ---------
rm(list = ls()[! ls() %in% c("RTK.hp.mean.FUSCCRNA.scale","ord")])
load("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/30.NG修改/rawdata/Luminal_9sample_Cancer鉴定后.Rdata")
raw = raw[ , ! raw@meta.data$major_V3 %in% c("Epithelial","Prolifer","Normal_epithelial")] ; gc()
raw@meta.data$major_V3[raw@meta.data$major_V3 %in% c("CD4T","CD8T","ISG_T","NK")] = "TNK"
# dplot <- DotPlot(object = raw, features = ord, split.by = "SNF_subtype",group.by = "merge_Cell_marjor_V2",
#                  cols = rainbow(100))
dplot <- DotPlot(object = raw, features = ord, group.by = "major_V3",
                 cols = rainbow(100))


ddata <- as.data.frame(dplot[["data"]])
df_spread <- tidyr::spread(ddata[, c(3,4,5)], id, avg.exp.scaled)
rownames(df_spread) <- df_spread[, 1]
df_spread <- df_spread[, -1]
mat <- as.matrix(df_spread)
mat <- na.omit(mat)
mat = mat[,order(colnames(mat))]
scRNA.mean.hp.scale = t(scale(t(mat)))
scRNA.mean.hp.scale = scRNA.mean.hp.scale[ord,]
# phet<- pheatmap(scRNA_mean, fontsize = 6, cellheight = 6, border_color = T,#annotation_row = anno_row,
#                 fontsize_row = 6,cluster_cols = F,cluster_rows = F)
# phet

scRNA.mean.hp.scale.annot = data.frame(row.names = colnames(scRNA.mean.hp.scale),
                                       CellType = str_split_fixed( colnames(scRNA.mean.hp.scale),"_",2)[,1]#,
                                       #SNF = str_split_fixed( colnames(scRNA.mean.hp.scale),"_",2)[,2]
)

########################### PART 3. CPTAC SNF matched RTK phos SNF diff annot ---------
rm(list = ls()[! ls() %in% c("RTK.hp.mean.FUSCCRNA.scale","ord","scRNA.mean.hp.scale","scRNA.mean.hp.scale.annot")])

## input DEphos
DEphos = read.csv("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/result/inferSNF_res_EachSNFfixed_V2/CPTAC_2022/DEPhos/DEPhos.Wilcox.CPTAC.SNF4VsOthers.csv",
                  row.names = 1)
DEphos$index = paste0(DEphos$geneSymbol,"_",DEphos$variableSites)
DEphos = na.omit(DEphos)
## pick best phos site for SNF4 in given gene list 
DEphos_f = DEphos[rownames(DEphos)[DEphos$geneSymbol %in% ord],]
DEphos_f$index2 = paste0(DEphos_f$geneSymbol,DEphos_f$Diff)

tmp <- DEphos_f %>% 
  group_by(geneSymbol) %>%    
  summarise(max = max(Diff,na.rm = T))
tmp$index = paste0(tmp$geneSymbol,tmp$max)

DEphos_f = DEphos_f[match(tmp$index,DEphos_f$index2),]
rownames(DEphos_f) = DEphos_f$geneSymbol

DEphos_f = DEphos_f[match(ord,rownames(DEphos_f)),]
rownames(DEphos_f) = ord 
rownames(DEphos_f)[na.omit(match(DEphos_f$geneSymbol,rownames(DEphos_f)))] = na.omit(DEphos_f[,"index"])
DEphos_f$geneSymbol = ord

DEphos_f = DEphos_f[match(rownames(RTK.hp.mean.FUSCCRNA.scale),DEphos_f$geneSymbol),]
#RTK.hp.mean.FUSCCRNA.scale.annot = DEphos_f[rownames(RTK.hp.mean.FUSCCRNA.scale),"Diff"]
RTK.hp.mean.FUSCCRNA.scale.annot = data.frame(row.names = rownames(DEphos_f),
                                              Diff = DEphos_f[,"Diff"])
########################### PART 4. merge & plot  ---------
rm(list = ls()[!ls() %in% c("RTK.hp.mean.FUSCCRNA.scale","RTK.hp.mean.FUSCCRNA.scale.annot",
                            "scRNA.mean.hp.scale","scRNA.mean.hp.scale.annot")])
color_fun_SNF= c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A")
names(color_fun_SNF) = paste0("SNF",1:4)

color_fun_CellType = as.character(paletteer_d("ggsci::default_nejm")[-4])
names(color_fun_CellType) = c("B", "CAF", "Endothelial", "Myeloid", "NormalEpi", "TNK", "Tumor")

col_fun_FUSCCRNA.scale = colorRamp2(c(-1.2,-0.6, 0 , 0.6,1.2), rev(brewer.pal(n = 5, name ="RdBu")))
col_fun_scRNA.scale = colorRamp2(c(-2,-1, 0 , 1,2), rev(brewer.pal(n = 5, name ="PiYG")))
color_fun_phos_diff = colorRamp2(c(-1.75, -0.875, 0 , 0.875,1.75), rev(brewer.pal(n = 5, name ="PRGn")))

annot = data.frame(row.names = paste0("SNF",1:4),SNF = paste0("SNF",1:4))


p1 = Heatmap(RTK.hp.mean.FUSCCRNA.scale,
             na_col = "#eaeaea",
             #border = T, 
             #show_row_names = F, 
             show_column_names = F,
             show_heatmap_legend = T, 
             row_names_side = "left",
             row_names_gp = gpar(fontsize = 7),
             col = col_fun_FUSCCRNA.scale,
             heatmap_legend_param = list(title = "bulk scaled \n averaged expression"),
             cluster_rows = FALSE, cluster_columns = FALSE,
             top_annotation = HeatmapAnnotation(SNF = annot$SNF,
                                                show_annotation_name = F,
                                                annotation_height = 7,
                                                simple_anno_size = unit(3, "mm"),
                                                #annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 7),
                                                col = list(SNF = color_fun_SNF)) )

{
  CelltypeOrd = c("Tumor","Myeloid","TNK","B","Endothelial","CAF")
  scRNA.mean.hp.scale = scRNA.mean.hp.scale[,CelltypeOrd]
  scRNA.mean.hp.scale.annot = data.frame( row.names = CelltypeOrd,CellType =  CelltypeOrd)
}
p2 = Heatmap(scRNA.mean.hp.scale,
             na_col = "#eaeaea",
             #border = T, 
             show_row_names = F,  show_column_names = F,
             #row_names_side = "left",
             row_names_gp = gpar(fontsize = 7),
             heatmap_legend_param = list(title = "scRNA-seq scaled \n averaged expression"),
             col = col_fun_scRNA.scale,
             cluster_rows = FALSE, cluster_columns = FALSE,
             column_split = factor(scRNA.mean.hp.scale.annot[,"CellType"],levels = c("Tumor","NormalEpi","Myeloid","TNK","B","Endothelial","CAF")),
             column_gap = unit(2, "mm"),
             column_title_gp = gpar(fontsize = 7),
             column_title_rot = 0,
             width = 5,
             top_annotation = HeatmapAnnotation(CellType = scRNA.mean.hp.scale.annot$CellType,
                                                #SNF = scRNA.mean.hp.scale.annot$SNF,
                                                show_annotation_name = F,
                                                #annotation_name_side = "right", 
                                                annotation_name_gp = gpar(fontsize = 7),
                                                simple_anno_size = unit(3, "mm"),
                                                col = list(#SNF = color_fun_SNF,
                                                  CellType = color_fun_CellType)) 
) 
p2
p3 = pheatmap(RTK.hp.mean.FUSCCRNA.scale.annot, 
              na_col = "#eaeaea",
              show_colnames = F,
              border_color = "white", 
              display_numbers = T,
              fontsize = 7,
              cellwidth = 20,
              col = color_fun_phos_diff,
              heatmap_legend_param = list(title = "Diff for SNF4"),
              cluster_rows = FALSE, cluster_cols = FALSE)

p = draw(p1+p2+p3, auto_adjust = F,
         ht_gap = unit(c(5, 2), "mm"))

export::graph2pdf(p,file = "RTK.hp.mean.SNF4phos.pdf",width = 8, height =8)

# ha1 = rowAnnotation(phos_diff = RTK.hp.mean.FUSCCRNA.scale.annot,
#                     annotation_name_gp = gpar(fontsize = 7),
#                     na_col = "#eaeaea",
#                     col = list(phos_diff = color_fun_phos_diff) )
# p = p1+ha1
# draw(p1+ha1, auto_adjust = F)
