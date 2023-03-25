rm(list = ls()) ; graphics.off()
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
library(pheatmap) ; library(RColorBrewer)
##### fdr val for pheatmap validation---
## fdr_RNA_NMF_351 
tmp1 = read.csv("multiomic_feature.pdf_clinic_aov.csv",
                row.names = 1)
tmp2 = read.csv("multiomic_feature.pdf_mut_chisq.csv",
                row.names = 1)
tmp3 = read.csv("multiomic_feature.pdf_AMP_chisq.csv",
                row.names = 1)
tmp4 = read.csv("multiomic_feature.pdf_path_aov.csv",
                row.names = 1)
tmp5 = read.csv("multiomic_feature.pdf_metab_aov.csv",
                row.names = 1)

fdr_RNA_NMF_351 = -log10(c(tmp1$adj.p, tmp2$adj.p, tmp3$adj.p, tmp4$adj.p,tmp5$adj.p))
names(fdr_RNA_NMF_351) = c(rownames(tmp1) ,rownames(tmp2),rownames(tmp3),rownames(tmp4),rownames(tmp5))
fdr_RNA_NMF_351[fdr_RNA_NMF_351 >= 10] = 10
fdr_RNA_NMF_351[fdr_RNA_NMF_351 <= -log10(0.05) ] = NA

## fdr_CNA_RNA
tmp1 = read.csv("multiomic_feature.pdf_clinic_aov.csv",
                row.names = 1)
tmp2 = read.csv("multiomic_feature.pdf_mut_chisq.csv",
                row.names = 1)
tmp3 = read.csv("multiomic_feature.pdf_AMP_chisq.csv",
                row.names = 1)
tmp4 = read.csv("multiomic_feature.pdf_path_aov.csv",
                row.names = 1)
tmp5 = read.csv("multiomic_feature.pdf_metab_aov.csv",
                row.names = 1)

fdr_RNA_CNA = -log10(c(tmp1$adj.p, tmp2$adj.p, tmp3$adj.p, tmp4$adj.p,tmp5$adj.p))
names(fdr_RNA_CNA) = c(rownames(tmp1) ,rownames(tmp2),rownames(tmp3),rownames(tmp4),rownames(tmp5))
fdr_RNA_CNA[fdr_RNA_CNA >= 10] = 10
fdr_RNA_CNA[fdr_RNA_CNA <= -log10(0.05) ] = NA


## fdr WES_RNA_CNA_Polar
tmp1 = read.csv("multiomic_feature.pdf_clinic_aov.csv",
                row.names = 1)
tmp2 = read.csv("multiomic_feature.pdf_mut_chisq.csv",
                row.names = 1)
tmp3 = read.csv("multiomic_feature.pdf_AMP_chisq.csv",
                row.names = 1)
tmp4 = read.csv("multiomic_feature.pdf_path_aov.csv",
                row.names = 1)
tmp5 = read.csv("multiomic_feature.pdf_metab_aov.csv",
                row.names = 1)

fdr_WES_RNA_CNA_Polar = -log10(c(tmp1$adj.p, tmp2$adj.p, tmp3$adj.p, tmp4$adj.p,tmp5$adj.p))
names(fdr_WES_RNA_CNA_Polar) = c(rownames(tmp1) ,rownames(tmp2),rownames(tmp3),rownames(tmp4),rownames(tmp5))
fdr_WES_RNA_CNA_Polar[fdr_WES_RNA_CNA_Polar >= 10] = 10
fdr_WES_RNA_CNA_Polar[fdr_WES_RNA_CNA_Polar <= -log10(0.05) ] = NA


## fdr_oriSNF
tmp1 = read.csv("Fig2_multiomic_feature.pdf_clinic_aov.csv", row.names = 1)
tmp2 = read.csv("Fig2_multiomic_feature.pdf_mut_chisq.csv", row.names = 1)
tmp3 = read.csv("Fig2_multiomic_feature.pdf_AMP_chisq.csv", row.names = 1)
tmp4 = read.csv("Fig2_multiomic_feature.pdf_path_aov.csv", row.names = 1)
tmp5 = read.csv("Fig2_multiomic_feature.pdf_metab_aov.csv", row.names = 1)


fdr_oriSNF = -log10(c(tmp1$adj.p, tmp2$adj.p, tmp3$adj.p, tmp4$adj.p,tmp5$adj.p))
names(fdr_oriSNF) = c(rownames(tmp1) ,rownames(tmp2),rownames(tmp3),rownames(tmp4),rownames(tmp5))
fdr_oriSNF[fdr_oriSNF >= 10] = 10
fdr_oriSNF[fdr_oriSNF <= -log10(0.05) ] = NA



## fdr Pro_RNA_CNA_Polar
tmp1 = read.csv("multiomic_feature.pdf_clinic_aov.csv",
                row.names = 1)
tmp2 = read.csv("multiomic_feature.pdf_mut_chisq.csv",
                row.names = 1)
tmp3 = read.csv("multiomic_feature.pdf_AMP_chisq.csv",
                row.names = 1)
tmp4 = read.csv("multiomic_feature.pdf_path_aov.csv",
                row.names = 1)
tmp5 = read.csv("multiomic_feature.pdf_metab_aov.csv",
                row.names = 1)

fdr_Pro_RNA_CNA_Polar = -log10(c(tmp1$adj.p, tmp2$adj.p, tmp3$adj.p, tmp4$adj.p,tmp5$adj.p))
names(fdr_Pro_RNA_CNA_Polar) = c(rownames(tmp1) ,rownames(tmp2),rownames(tmp3),rownames(tmp4),rownames(tmp5))
fdr_Pro_RNA_CNA_Polar[fdr_Pro_RNA_CNA_Polar >= 10] = 10
fdr_Pro_RNA_CNA_Polar[fdr_Pro_RNA_CNA_Polar <= -log10(0.05) ] = NA


##
fdr_Pro_RNA_CNA_Polar = fdr_Pro_RNA_CNA_Polar[intersect(names(fdr_Pro_RNA_CNA_Polar),names(fdr_oriSNF))]

fdr_RNA_NMF_351 = as.data.frame(fdr_RNA_NMF_351[names(fdr_Pro_RNA_CNA_Polar)])
fdr_RNA_CNA = as.data.frame(fdr_RNA_CNA[names(fdr_Pro_RNA_CNA_Polar)])
fdr_oriSNF = as.data.frame(fdr_oriSNF[names(fdr_Pro_RNA_CNA_Polar)])
fdr_WES_RNA_CNA_Polar = as.data.frame(fdr_WES_RNA_CNA_Polar[names(fdr_Pro_RNA_CNA_Polar)])
fdr_Pro_RNA_CNA_Polar = as.data.frame(fdr_Pro_RNA_CNA_Polar)

fdr = data.frame(row.names = row.names(fdr_RNA_NMF_351),
                 RNA_NMF_351 = fdr_RNA_NMF_351[,1],
                 RNA_CNA = fdr_RNA_CNA[,1],
                 oriSNF = fdr_oriSNF[,1],
                 WES_RNA_CNA_Polar = fdr_WES_RNA_CNA_Polar[,1],
                 Pro_RNA_CNA_Polar = fdr_Pro_RNA_CNA_Polar[,1])
pheatmap::pheatmap(fdr,cluster_rows = F,cluster_cols = F,display_numbers = T,
                   colorRampPalette(brewer.pal(n = 7, name = "Greys")[1:5])(100), #YlOrRd
                   na_col = "#eaeaea",
                   border  = NA,filename = "./fdr_hp.pdf",height = 8,width = 7.8)



