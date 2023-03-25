rm(list = ls())
graphics.off()
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(RColorBrewer)
library(GSVA) ; library(GSEABase)

load("./CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
exp.fpkm.TT.log = log2(exp.fpkm.TT+1)

luminal$SNF = paste0("SNF",SNF_Cluster[rownames(luminal)])
luminal$SNF2orNot = luminal$SNF
luminal$SNF2orNot[luminal$SNF2orNot != "SNF2" ]  = "nonSNF2"

#### annot bar -------
{
  res_estimate = data.table::fread("./SNF_estimate_estimate_score.gct")
  res_estimate = as.data.frame(t(res_estimate[3:6,-1]))
  colnames(res_estimate) = res_estimate[1,]
  rownames(res_estimate) = res_estimate[,1]
  res_estimate = res_estimate[-1,-1]
  #res_estimate$StromalScore = as.numeric(scale(as.numeric(res_estimate$StromalScore)))
  res_estimate$ImmuneScore = as.numeric(res_estimate$ImmuneScore)
}
luminal$ImmuneSignature = res_estimate[rownames(luminal),"ImmuneScore"]

#
{
  CINscore = as.data.frame(matrix(nrow = length(SNF_Cluster), ncol = 2))
  rownames(CINscore) = names(sort(SNF_Cluster))  ;colnames(CINscore) = c("CINscore","SNF_All")
  CINscore$SNF_All = paste0("SNF",sort(SNF_Cluster))
  for (i in names(sort(SNF_Cluster))){
    tmp = cnv.alldata[,i]
    tmp = tmp^2
    CINscore[i,"CINscore"] = sum(tmp)
  }
  CINscore$CINscore = log2(CINscore$CINscore)
}
luminal$CINscore = CINscore[rownames(luminal),"CINscore"]

#
luminal$PAM50_PC_TCGA = factor(luminal$PAM50_classifier, levels =c( "LumA" ,"LumB" ,"Her2","Basal" ,"Normal"))

#
sTILS = openxlsx::read.xlsx("./CBCGA_TILs_F.xlsx",1)
rownames(sTILS) = sTILS$PatientCode
sTILS = sTILS[rownames(luminal),]
luminal = cbind(luminal,sTILS)

#
Tls.sig = read.csv("./ICB_signature.csv",row.names = 1)
Tls.sig = Tls.sig[rownames(luminal),]
luminal = cbind(luminal,Tls.sig)

#
TCR_BCR = read.csv("./cbcga_927sample_TCRBCR_220719.csv",row.names = 1)
  TCR_BCR = TCR_BCR[str_detect(TCR_BCR$ID , "_T"),] 
  TCR_BCR = TCR_BCR[substring(TCR_BCR$sample_id_detail,1,1) == "B",] #'KTWZ', 'RVAQ' 重复测了RNA
  TCR_BCR = TCR_BCR[!is.na(TCR_BCR$ID),]
  rownames(TCR_BCR) = TCR_BCR$Patient_ID
TCR_BCR = TCR_BCR[rownames(luminal),]
luminal = cbind(luminal,TCR_BCR)

#
protein_log2_NA = as.data.frame(matrix(NA,nrow = nrow(protein_log2), 
                           ncol = length(setdiff(rownames(luminal),colnames(protein_log2))) ))
rownames(protein_log2_NA) = rownames(protein_log2)
colnames(protein_log2_NA) = setdiff(rownames(luminal),colnames(protein_log2))


#### mat1 ---------
tmp = c("ImmuneSignature","TIS.signature","STAT1.signature")
mat1 = t(luminal[,tmp])
mat1 = t(scale(t(mat1)))

#### mat2 ---------
tmp = c("PDCD1","CD274","CTLA4")
mat2_RNA = exp.fpkm.TT.log[tmp, rownames(luminal)]
mat2_RNA = t(scale(t(mat2_RNA)))
rownames(mat2_RNA) = paste0(rownames(mat2_RNA),"_RNA")  

mat2 = mat2_RNA

#### mat3 ---------
tmp = c("TRA_diversity","TRB_diversity")
mat3 =  t(luminal[,tmp])
mat3 = t(scale(t(mat3)))

#### mat4 ---------
tmp = c("CD8A","GZMA","PRF1","IDO1")
mat4_RNA = exp.fpkm.TT.log[tmp,rownames(luminal) ]
mat4_RNA = t(scale(t(mat4_RNA)))
rownames(mat4_RNA) = paste0(rownames(mat4_RNA),"_RNA")

mat4_Pro = protein_log2[intersect(rownames(protein_log2),tmp), ]
mat4_Pro = t(scale(t(mat4_Pro)))
mat4_Pro = cbind(mat4_Pro,protein_log2_NA[rownames(mat4_Pro),])
mat4_Pro = mat4_Pro[,rownames(luminal) ]
rownames(mat4_Pro) = paste0(rownames(mat4_Pro),"_Protein")

mat4_Polar = polar_metabolite_TT_MS2_log2[rownames(polar_metabolite_MS2_mapping)[polar_metabolite_MS2_mapping$metabolite_mapping_name == "L-Kynurenine"],
                                          rownames(luminal)]
mat4_Polar = t(scale(t(mat4_Polar)))
rownames(mat4_Polar) = "L-Kynurenine"
mat4 = rbind(mat4_RNA,mat4_Pro,mat4_Polar)

mat4 = mat4[c("CD8A_RNA","CD8A_Protein",
              "GZMA_RNA","GZMA_Protein",
              "PRF1_RNA","PRF1_Protein",
              "IDO1_RNA","IDO1_Protein","L-Kynurenine"),]

#### mat5 ---------
geneset4imm<-read.csv("./2015-CIBERSORT.csv")
list<- split(as.matrix(geneset4imm)[,1], geneset4imm[,2])

gsva_matrix<- gsva(as.matrix(exp.fpkm.TT), list, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

tmp = pheatmap::pheatmap(gsva_matrix,cluster_rows = T,cluster_cols = F,scale = "row")

mat5 = gsva_matrix[tmp[["tree_row"]][["order"]], rownames(luminal)]
mat5 = t(scale(t(mat5)))

#### draw  ------
## ord
#ord = rownames(luminal[order(luminal$SNF2orNot,luminal$ImmuneSignature,decreasing = T),])
ord = rownames(luminal[order(luminal$SNF,decreasing = F),])
luminal = luminal[ord,]

## annot 
color_PAM50 = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
  names(color_PAM50) = c("LumA","LumB","Her2","Basal","Normal")
names(color) = paste0("SNF",1:4)
col_fun = colorRamp2(c(-2,-1,0 , 1,2), rev(brewer.pal(n = 5, name ="RdBu")))

col_fun_Isig <- colorRamp2(c(-0.5, 0.3 ,2), colorRampPalette(c("white", "#08506F"))(20)[c(2,10,20)])
col_fun_CIN <- colorRamp2(c(10, 13, 15), c("white", '#DBBEBE', '#B75C5B'))
col_fun_sTILS <- colorRamp2(c(0, 0.1, 0.35), c("white", '#FDF2C9', '#F8CA26'))

## heatmap
p1 = Heatmap(mat1[,ord],
        cluster_rows = FALSE,cluster_columns = FALSE,
        col = col_fun,
        na_col = "#eaeaea",
        border = T,
        show_column_names = FALSE,
        column_split =factor(luminal$SNF, levels = paste0("SNF",1:4)) ,
        column_gap = unit(2, "mm"),
        row_names_side = "left",
        heatmap_legend_param = list(title = "Z-score"),
        top_annotation = HeatmapAnnotation(SNF = luminal$SNF,
                                           PAM50 = luminal$PAM50_classifier,
                                           CINscore = luminal$CINscore,
                                           sTils = luminal$sTILs,
                                           na_col = "#eaeaea",
                                           gap = unit(1, "points"),
                                           annotation_name_side = "left",
                                           annotation_name_gp = gpar(fontsize = 9),
                                           col = list(PAM50 = color_PAM50,
                                                      SNF = color,
                                                      CINscore = col_fun_CIN,
                                                      sTils = col_fun_sTILS
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
             column_split =factor(luminal$SNF, levels = paste0("SNF",1:4)) ,
             column_gap = unit(2, "mm"))

p3 = Heatmap(mat3[,ord],
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             border = T,
             na_col = "#eaeaea",
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             show_column_names = FALSE,
             column_split =factor(luminal$SNF, levels = paste0("SNF",1:4)) ,
             column_gap = unit(2, "mm"))

p4 = Heatmap(mat4[,ord],
             na_col = "#eaeaea",
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             border = T,
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             show_column_names = FALSE,
             column_split =factor(luminal$SNF, levels = paste0("SNF",1:4)) ,
             column_gap = unit(2, "mm"))

p5 = Heatmap(mat5[,ord],
             na_col = "#eaeaea",
             cluster_rows = FALSE,cluster_columns = FALSE,
             col = col_fun,
             border = T,
             row_names_side = "left",
             show_heatmap_legend = FALSE,
             show_column_names = FALSE,
             column_split =factor(luminal$SNF, levels = paste0("SNF",1:4)) ,
             column_gap = unit(2, "mm"))


p = p1 %v% p2  %v% p3 %v% p4 %v% p5
p

export::graph2pdf(p,file = "Fig5.pdf",width = 12,height = 8)
export::graph2ppt(p,file = "Fig5.pptx",width = 12,height = 8)


### anova
mat_merge = rbind(as.data.frame(t(luminal[ord,c("CINscore","sTILs")])),mat1[,ord], mat2[,ord], mat3[,ord], mat4[,ord], mat5[,ord])
t.test.res = as.data.frame(matrix(nrow = nrow(mat_merge),ncol = 2))
rownames(t.test.res) = rownames(mat_merge) ; colnames(t.test.res) = c("p.val",'adj.p')
for ( i in rownames(mat_merge)){
  tmp0 = luminal$SNF ; names(tmp0) = row.names(luminal)
  tmp = data.frame(gene = as.numeric(mat_merge[i,rownames(luminal)]),
                   group = luminal$SNF)
  res = summary(aov(gene~group,data = tmp))
  t.test.res[i,"p.val"] = res[[1]][["Pr(>F)"]][1]
}
t.test.res$adj.p = p.adjust(t.test.res$p.val,method = "fdr")

write.csv(t.test.res,file = "aov_res.csv")

