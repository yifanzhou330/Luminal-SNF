rm(list = ls())
graphics.off()
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
load("CBCGA_HRposHER2neg354_WES_RNAseq_CNV_Metab_Protein_20211127.Rdata")
exp.fpkm.TT.log = log2(exp.fpkm.TT+1)
#### add infor  -------
MGPS = read.csv("MGPS_gene.csv")[,2]
MGPS = exp.fpkm.TT[intersect(MGPS,rownames(exp.fpkm.TT)),]
MGPS = apply(MGPS,2,mean)

CINscore = as.data.frame(matrix(nrow = length(SNF_Cluster), ncol = 2))
rownames(CINscore) = names(sort(SNF_Cluster))  ;colnames(CINscore) = c("CINscore","SNF_All")
CINscore$SNF_All = paste0("SNF",sort(SNF_Cluster))
for (i in names(sort(SNF_Cluster))){
  tmp = cnv.alldata[,i]
  tmp = tmp^2
  CINscore[i,"CINscore"] = sum(tmp)
}
CINscore$CINscore = log2(CINscore$CINscore)

#### annot bar -------
luminal$PAM50_PC_TCGA = factor(luminal$PAM50_PC_TCGA, levels =c( "LumA" ,"LumB" ,"Her2","Basal" ,"Normal"))
luminal$HRD = as.numeric(luminal$HRD)
luminal$MGPS = MGPS[rownames(luminal)]
luminal$MGPS.scale = scale(luminal$MGPS)
luminal$SNF3orNot = luminal$SNF分型
luminal$SNF3orNot[luminal$SNF3orNot != "SNF3"] = "non_SNF3"
luminal$CINscore = CINscore[rownames(luminal),"CINscore"]

protein_log2_NA = as.data.frame(matrix(NA,nrow = nrow(protein_log2), 
                                       ncol = length(setdiff(rownames(luminal),colnames(protein_log2))) ))
rownames(protein_log2_NA) = rownames(protein_log2)
colnames(protein_log2_NA) = setdiff(rownames(luminal),colnames(protein_log2))


#### mat1 ---------
tmp = c("CCND1","CDK2") #tmp = c("CCND1","CDK1","CDK2")
mat1_RNA = exp.fpkm.TT.log[tmp, names(sort(SNF_Cluster))]
mat1_RNA = t(scale(t(mat1_RNA)))
rownames(mat1_RNA) = paste0(rownames(mat1_RNA),"_RNA")

# mat1_CNV = cnv.thre[tmp, names(sort(SNF_Cluster))]
# mat1_CNV = apply(mat1_CNV, 2, function(x){ifelse(x < 1, 0, 2)} )
mat1_CNV = cnv.alldata[tmp, names(sort(SNF_Cluster))]
mat1_CNV = apply(mat1_CNV, 2, function(x){ifelse(x < log2(4/2), 0, 2)} )
rownames(mat1_CNV) = paste0(rownames(mat1_CNV),"_CNA")

mat1_Pro = protein_log2[intersect(tmp,rownames(protein_log2)), ]
mat1_Pro = t(scale(t(mat1_Pro)))
mat1_Pro = cbind(mat1_Pro,protein_log2_NA[rownames(mat1_Pro),])
mat1_Pro = mat1_Pro[,names(sort(SNF_Cluster))]
rownames(mat1_Pro) = paste0(rownames(mat1_Pro),"_Protein")

mat1 = rbind(mat1_CNV,mat1_RNA,mat1_Pro)
mat1 = mat1[c("CCND1_CNA","CCND1_RNA","CCND1_Protein",
#              "CDK1_CNA","CDK1_RNA","CDK1_Protein",
              "CDK2_CNA","CDK2_RNA","CDK2_Protein"),]

#### mat2 ---------
mat2_RNA = exp.fpkm.TT.log[c("CDK1"), names(sort(SNF_Cluster))]
mat2_RNA = t(scale(t(mat2_RNA)))
rownames(mat2_RNA) = paste0(rownames(mat2_RNA),"_RNA")

mat2_CNV = cnv.alldata[c("CDK1"), names(sort(SNF_Cluster))]
mat2_CNV = t(as.data.frame(apply(mat2_CNV, 2, function(x){ifelse(x < log2(4/2), 0, 2)} )))
rownames(mat2_CNV) = "CDK1_CNA"

mat2_Pro = protein_log2[c("CDK1"), ]
mat2_Pro = t(scale(t(mat2_Pro)))
mat2_Pro = cbind(mat2_Pro,protein_log2_NA[c("CDK1"),])
mat2_Pro = mat2_Pro[,names(sort(SNF_Cluster))]
rownames(mat2_Pro) = paste0(rownames(mat2_Pro),"_Protein")

mat2 = rbind(mat2_RNA,mat2_CNV,mat2_Pro)


#### mat3 ---------
tmp = c("AURKA","MYBL2","TOP2A","ESPL1","DSCC1","GINS4","RAD21","KIF18B","E2F1","E2F2")
mat3_RNA = exp.fpkm.TT.log[tmp, names(sort(SNF_Cluster))]
mat3_RNA = t(scale(t(mat3_RNA)))
rownames(mat3_RNA) = paste0(rownames(mat3_RNA),"_RNA")

mat3_Pro = protein_log2[intersect(rownames(protein_log2),tmp), ]
mat3_Pro = t(scale(t(mat3_Pro)))
mat3_Pro = cbind(mat3_Pro,protein_log2_NA[rownames(mat3_Pro),])
mat3_Pro = mat3_Pro[,names(sort(SNF_Cluster))]
rownames(mat3_Pro) = paste0(rownames(mat3_Pro),"_Protein")

mat3 = rbind(mat3_RNA,mat3_Pro)
mat3 = mat3[sort(rownames(mat3)),]
mat3 = mat3[c(2,1,4,3,5,6,7,9,8,10,11,13,12,15,14) ,]

#### draw  ------
## ord
ord = rownames(luminal[order(luminal$SNF3orNot,luminal$MGPS.scale,decreasing = F),])
luminal = luminal[ord,]

## annot 
color_PAM50 = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
names(color_PAM50) = c("LumA","LumB","Her2","Basal","Normal")
names(color) = paste0("SNF",1:4)
col_fun = colorRamp2(c(-2,-1,0 , 1,2), rev(brewer.pal(n = 5, name ="RdBu")))
#lgd = Legend(col_fun = col_fun,title = "z-score")
# col_fun_ki67 <- colorRamp2(c(3, 30, 98), c('#EBEFF6', '#99AFD2', '#3C74AE'))
col_fun_CIN <- colorRamp2(c(10, 13, 15), c("white", '#DBBEBE', '#B75C5B'))
col_fun_MGPS <-  colorRamp2(c(-1.5,-0.5, 0 , 0.5,1.5), rev(brewer.pal(n = 5, name ="RdBu")))
col_fun_Isig <- colorRamp2(c(1, 1.1 ,1.2), colorRampPalette(c("white", "#08506F"))(20)[c(2,10,20)])
col_fun_SSig <- colorRamp2(c(1, 1.1 ,1.2), colorRampPalette(c("white", "#50A9AC"))(20)[c(2,10,20)])
col_fun_hrd <- colorRamp2(c(0, 16, 75), c('#FDF9E7', '#EFE089', '#D7C801'))

## heatmap
p1 = Heatmap(mat1[,ord],
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             na_col = "#eaeaea",
             border = T,
             show_column_names = FALSE,
             column_split =factor(luminal$SNF3orNot, levels = c("SNF3","non_SNF3")) ,
             column_gap = unit(2, "mm"),
             row_names_side = "left",
             heatmap_legend_param = list(title = "Z-score"),
             top_annotation = HeatmapAnnotation(SNF = luminal$SNF分型,
                                                PAM50 = luminal$PAM50_PC_TCGA,
                                                #ImmuneSignature = luminal$ImmuneSignature,
                                                #StromalSignature = luminal$StromalSignature,
                                                CINscore = luminal$CINscore,
                                                HRDscore = luminal$HRD,
                                                MGPS = luminal$MGPS.scale,
                                                # Age.at.surgery = CBCGAClin$Age,
                                                # Lymph_node_status = CBCGAClin$Lymph_node_status,
                                                # Ki67 = as.numeric(CBCGAClin$Ki67),
                                                #Relapse = luminal$DMFS_status ,
                                                #border = T,
                                                gap = unit(1, "points"),
                                                annotation_name_side = "left",
                                                annotation_name_gp = gpar(fontsize = 9),
                                                col = list(PAM50 = color_PAM50,
                                                           SNF = color,
                                                           #StromalSignature = col_fun_SSig,
                                                           #ImmuneSignature =col_fun_Isig,
                                                           CINscore = col_fun_CIN,
                                                           #Age.at.surgery = col_fun_age,
                                                           #Lymph_node_status = c('0' = 'white', '1' = 'black'),
                                                           #Ki67 = col_fun_ki67,
                                                           HRDscore = col_fun_hrd,
                                                           MGPS = col_fun_MGPS
                                                           #Relapse = c('0' = 'white', '1' = 'black')
                                                ))
)
p2 = Heatmap(mat2[,ord],
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             na_col = "#eaeaea",
             border = T,
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             show_column_names = FALSE,
             column_split = factor(luminal$SNF分型, levels = c('SNF1', 'SNF2', 'SNF3', 'SNF4')),
             column_gap = unit(2, "mm"))

p3 = Heatmap(mat3[,ord],
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             na_col = "#eaeaea",
             border = T,
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             show_column_names = FALSE,
             column_split = factor(luminal$SNF分型, levels = c('SNF1', 'SNF2', 'SNF3', 'SNF4')),
             column_gap = unit(2, "mm"))

p = p1 %v% p2 %v% p3
p

export::graph2pdf(p,file = "Fig4_alldata4_V2.pdf",width = 12,height = 7)
export::graph2ppt(p,file = "Fig4_alldata4_V2.pptx",width = 12,height = 7)




