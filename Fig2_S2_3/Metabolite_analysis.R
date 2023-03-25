#-----------------------------------------------------------------#
# Fig3 
#-----------------------------------------------------------------#
rm(list = ls()) ; graphics.off()
options(stringsAsFactors = F)
set.seed(123)

install.packages("eoffice") 
library(Rtsne)
#library(ggpubr)
library(ggplot2)
library(dplyr)

####Load data
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
load("Fig5_S5.Rdata")

length(unique(metabolic_gene_list$b)) #4134
# write.csv(metabolic_gene_list,"metabolic_gene_list.csv")
mypro <- intersect(unique(metabolic_gene_list$b),rownames(protein_log2)) #2252

#####################################################################################################################
#####################################################################################################################
####------------------Fig 5C-D--------------####
#####################################################################################################################
#####################################################################################################################
## 1. Data filtering and cleaning
#####################################################################################################################

lumimal_SNF <- as.data.frame(cbind(luminal$PatientCode,SNF_Cluster)) 
colnames(lumimal_SNF) <- c("PatientCode","SNF_subtype")

SNF1_ID <- lumimal_SNF[lumimal_SNF$SNF_subtype == "1",1] # 86
SNF2_ID <- lumimal_SNF[lumimal_SNF$SNF_subtype == "2",1] # 89
SNF3_ID <- lumimal_SNF[lumimal_SNF$SNF_subtype == "3",1] # 118
SNF4_ID <- lumimal_SNF[lumimal_SNF$SNF_subtype == "4",1] # 58

SNF1_PRO <- as.data.frame(protein_log2[,intersect(colnames(protein_log2),SNF1_ID)]) # 40
SNF2_PRO <- as.data.frame(protein_log2[,intersect(colnames(protein_log2),SNF2_ID)]) # 46
SNF3_PRO <- as.data.frame(protein_log2[,intersect(colnames(protein_log2),SNF3_ID)]) # 64
SNF4_PRO <- as.data.frame(protein_log2[,intersect(colnames(protein_log2),SNF4_ID)]) # 29
PT_PRO <- protein_PT_log2 # 48

SNF1_PRO <- SNF1_PRO[mypro,] # 40
SNF2_PRO <- SNF2_PRO[mypro,] # 46
SNF3_PRO <- SNF3_PRO[mypro,] # 64
SNF4_PRO <- SNF4_PRO[mypro,] # 29

TT_PRO <- cbind(SNF1_PRO,SNF2_PRO,SNF3_PRO,SNF4_PRO) # 179
PT_PRO <- PT_PRO[mypro,] # 48

#rna
exp.fpkm.TT_meta <- exp.fpkm.TT[rownames(exp.fpkm.TT)%in%KEGG_Meta_signature$Genes,]

SNF1_rna <- exp.fpkm.TT_meta[,intersect(colnames(exp.fpkm.TT_meta),SNF1_ID)] #86
SNF2_rna <- exp.fpkm.TT_meta[,intersect(colnames(exp.fpkm.TT_meta),SNF2_ID)] #89
SNF3_rna <- exp.fpkm.TT_meta[,intersect(colnames(exp.fpkm.TT_meta),SNF3_ID)] #118
SNF4_rna <- exp.fpkm.TT_meta[,intersect(colnames(exp.fpkm.TT_meta),SNF4_ID)] #58

TT_rna <- cbind(SNF1_rna,SNF2_rna,SNF3_rna,SNF4_rna)
PT_rna <- exp.fpkm.PT

#pol
SNF1_pol <- polar_metabolite_TT_MS2_log2[,intersect(colnames(polar_metabolite_TT_MS2_log2),SNF1_ID)] # 86
SNF2_pol <- polar_metabolite_TT_MS2_log2[,intersect(colnames(polar_metabolite_TT_MS2_log2),SNF2_ID)] # 89
SNF3_pol <- polar_metabolite_TT_MS2_log2[,intersect(colnames(polar_metabolite_TT_MS2_log2),SNF3_ID)] # 118
SNF4_pol <- polar_metabolite_TT_MS2_log2[,intersect(colnames(polar_metabolite_TT_MS2_log2),SNF4_ID)] # 58
TT_pol <- cbind(SNF1_pol,SNF2_pol,SNF3_pol,SNF4_pol) # 351
PT_pol <- polar_metabolite_PT_MS2_log2 # 28
write.csv(TT_pol,"TT_pol.CSV")
write.csv(PT_pol,"PT_pol.CSV")


SNF1_lip <- lipid_TT_MS2_log2[,intersect(colnames(lipid_TT_MS2_log2),SNF1_ID)] # 86
SNF2_lip <- lipid_TT_MS2_log2[,intersect(colnames(lipid_TT_MS2_log2),SNF2_ID)] # 89
SNF3_lip <- lipid_TT_MS2_log2[,intersect(colnames(lipid_TT_MS2_log2),SNF3_ID)] # 118
SNF4_lip <- lipid_TT_MS2_log2[,intersect(colnames(lipid_TT_MS2_log2),SNF4_ID)] # 58

TT_lip <- cbind(SNF1_lip,SNF2_lip,SNF3_lip,SNF4_lip) # 351
PT_lip <- lipid_PT_MS2_log2 # 28
# save(TT_PRO,PT_PRO,TT_pol,PT_pol,TT_lip,PT_lip,file = "TT_PT_DATA.Rdata")


#####################################################################################################################
## 2. Metabolic protein and polar metabolite subtype-specific analysis
#####################################################################################################################

##### metabolic protein #####
comparison_matrix <- matrix(ncol=13,nrow=length(mypro))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- mypro

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_PRO[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_PRO[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_PRO[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_PRO[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:4]))-min(as.numeric(comparison_matrix[i,1:4]))
  comparison_matrix[i,"SD"] <- sd(TT_PRO[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(SNF1_PRO[i,]),as.numeric(SNF2_PRO[i,]),as.numeric(SNF3_PRO[i,]),
                         as.numeric(SNF4_PRO[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_PRO[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_PRO[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_PRO[i,]),as.numeric(PT_PRO[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
luminal_PRO_SNF <- comparison_matrix
# write.csv(luminal_PRO_SNF,file="./results/luminal_PRO_SNF.csv")

##### rna #####
comparison_matrix <- matrix(ncol=13,nrow=nrow(TT_rna))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- rownames(TT_rna)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_rna[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_rna[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_rna[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_rna[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:4]))-min(as.numeric(comparison_matrix[i,1:4]))
  comparison_matrix[i,"SD"] <- sd(TT_rna[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(SNF1_rna[i,]),as.numeric(SNF2_rna[i,]),as.numeric(SNF3_rna[i,]),
                         as.numeric(SNF4_rna[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_rna[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_rna[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_rna[i,]),as.numeric(PT_rna[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
luminal_rna_SNF <- comparison_matrix

#filter:T_PT>0
luminal_rna_SNF <- luminal_rna_SNF[luminal_rna_SNF$T_PT>0,]

# write.csv(luminal_rna_SNF,file="./results/luminal_rna_SNF.csv")


##### Polar metabolite #####
comparison_matrix <- matrix(ncol=13,nrow=nrow(polar_metabolite_TT_MS2_log2))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- rownames(polar_metabolite_TT_MS2_log2)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_pol[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_pol[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_pol[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_pol[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:4]))-min(as.numeric(comparison_matrix[i,1:4]))
  comparison_matrix[i,"SD"] <- sd(TT_pol[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(SNF1_pol[i,]),as.numeric(SNF2_pol[i,]),as.numeric(SNF3_pol[i,]),
                         as.numeric(SNF4_pol[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_pol[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_pol[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_pol[i,]),as.numeric(PT_pol[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
luminal_Pol_SNF <- comparison_matrix
write.csv(luminal_Pol_SNF,file="./results/luminal_Pol_SNF.csv")

##### lip #####
comparison_matrix <- matrix(ncol=13,nrow=nrow(lipid_TT_MS2_log2))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- rownames(lipid_TT_MS2_log2)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_lip[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_lip[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_lip[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_lip[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:4]))-min(as.numeric(comparison_matrix[i,1:4]))
  comparison_matrix[i,"SD"] <- sd(TT_lip[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(SNF1_lip[i,]),as.numeric(SNF2_lip[i,]),as.numeric(SNF3_lip[i,]),
                         as.numeric(SNF4_lip[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_lip[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_lip[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_lip[i,]),as.numeric(PT_lip[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
luminal_lip_SNF <- comparison_matrix
# write.csv(luminal_lip_SNF,file="./results/luminal_lip_SNF.csv")

##### Lipid_cat #####
load("CBCGA.Extended_MergedData_V2.5_220722.Rdata")
lipid_cat <- unique(CBCGA_lip_anno$Subclass)

Cus_lipid_cat_TT <- matrix(ncol=ncol(TT_lip),nrow=length(lipid_cat))
rownames(Cus_lipid_cat_TT) <- lipid_cat
colnames(Cus_lipid_cat_TT) <- colnames(TT_lip)

Cus_lipid_cat_PT <- matrix(ncol=ncol(PT_lip),nrow=length(lipid_cat))
rownames(Cus_lipid_cat_PT) <- lipid_cat
colnames(Cus_lipid_cat_PT) <- colnames(PT_lip)

for (i in rownames(Cus_lipid_cat_TT)){
  peak <- rownames(CBCGA_lip_anno)[CBCGA_lip_anno$Subclass==i]
  Cus_lipid_cat_TT[i,] <- c(apply(TT_lip[peak,],2,mean))
}

for (i in rownames(Cus_lipid_cat_PT)){
  peak <- rownames(CBCGA_lip_anno)[CBCGA_lip_anno$Subclass==i]
  Cus_lipid_cat_PT[i,] <- c(apply(PT_lip[peak,],2,mean))
}

SNF1_lip_cat <- Cus_lipid_cat_TT[,intersect(colnames(lipid_TT_MS2_log2),SNF1_ID)] # 86
SNF2_lip_cat <- Cus_lipid_cat_TT[,intersect(colnames(lipid_TT_MS2_log2),SNF2_ID)] # 89
SNF3_lip_cat <- Cus_lipid_cat_TT[,intersect(colnames(lipid_TT_MS2_log2),SNF3_ID)] # 118
SNF4_lip_cat <- Cus_lipid_cat_TT[,intersect(colnames(lipid_TT_MS2_log2),SNF4_ID)] # 58

TT_lip_cat <- cbind(SNF1_lip_cat,SNF2_lip_cat,SNF3_lip_cat,SNF4_lip_cat) # 351
PT_lip_cat <- Cus_lipid_cat_PT

comparison_matrix <- matrix(ncol=13,nrow=nrow(TT_lip_cat))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4","max_min",
                                 "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(comparison_matrix) <- rownames(TT_lip_cat)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"max_min"] <- max(as.numeric(comparison_matrix[i,1:4]))-min(as.numeric(comparison_matrix[i,1:4]))
  comparison_matrix[i,"SD"] <- sd(TT_lip_cat[i,],na.rm = T)
  a <- kruskal.test(list(as.numeric(SNF1_lip_cat[i,]),as.numeric(SNF2_lip_cat[i,]),as.numeric(SNF3_lip_cat[i,]),
                         as.numeric(SNF4_lip_cat[i,])))
  comparison_matrix[i,"P"] <- a$p.value
  comparison_matrix[i,"Tumor"] <- mean(as.numeric(TT_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"PT"] <- mean(as.numeric(PT_lip_cat[i,]),na.rm = T)
  comparison_matrix[i,"T_PT"] <- as.numeric(comparison_matrix[i,"Tumor"])-as.numeric(comparison_matrix[i,"PT"])
  b <- wilcox.test(as.numeric(TT_lip_cat[i,]),as.numeric(PT_lip_cat[i,]))
  comparison_matrix[i,"P_TN"] <- b$p.value
}
comparison_matrix[,"FDR"] <- p.adjust(comparison_matrix[,"P"],method="fdr")
comparison_matrix[,"FDR_TN"] <- p.adjust(comparison_matrix[,"P_TN"],method="fdr")
comparison_matrix <- as.data.frame(comparison_matrix,stringsAsFactors = F)
luminal_lip_cat_SNF <- comparison_matrix
# write.csv(luminal_lip_cat_SNF,file="./results/luminal_lip_cat_SNF.csv")

##### Scale #####
scale_PRO<-t(scale(t(protein_log2),center=T,scale=T))
scale_Pol<-t(scale(t(polar_metabolite_TT_MS2_log2),center=T,scale=T))
scale_lip_cat<-t(scale(t(Cus_lipid_cat_TT),center=T,scale=T))
scale_rna<-t(scale(t(exp.fpkm.TT_meta),center=T,scale=T))
scale_rna <- scale_rna[rownames(scale_rna)%in%rownames(luminal_rna_SNF),]


#PRO
SNF1_matrix <- scale_PRO[,intersect(colnames(scale_PRO),SNF1_ID)]
SNF2_matrix <- scale_PRO[,intersect(colnames(scale_PRO),SNF2_ID)]
SNF3_matrix <- scale_PRO[,intersect(colnames(scale_PRO),SNF3_ID)]
SNF4_matrix <- scale_PRO[,intersect(colnames(scale_PRO),SNF4_ID)]

comparison_matrix <- matrix(ncol=4,nrow=nrow(scale_PRO))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4")
rownames(comparison_matrix) <- rownames(scale_PRO)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_matrix[i,]),na.rm = T)
}
comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
scale_pro_SNF<-comparison_matrix
# write.table(comparison_matrix,file="scale_pro_PAM50.txt",sep="\t")

#rna
SNF1_matrix <- scale_rna[,intersect(colnames(scale_rna),SNF1_ID)]
SNF2_matrix <- scale_rna[,intersect(colnames(scale_rna),SNF2_ID)]
SNF3_matrix <- scale_rna[,intersect(colnames(scale_rna),SNF3_ID)]
SNF4_matrix <- scale_rna[,intersect(colnames(scale_rna),SNF4_ID)]

comparison_matrix <- matrix(ncol=4,nrow=nrow(scale_rna))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4")
rownames(comparison_matrix) <- rownames(scale_rna)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_matrix[i,]),na.rm = T)
}
comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
scale_rna_SNF<-comparison_matrix
# write.table(comparison_matrix,file="scale_rna_PAM50.txt",sep="\t")

#POL
SNF1_matrix <- scale_Pol[,intersect(colnames(scale_Pol),SNF1_ID)]
SNF2_matrix <- scale_Pol[,intersect(colnames(scale_Pol),SNF2_ID)]
SNF3_matrix <- scale_Pol[,intersect(colnames(scale_Pol),SNF3_ID)]
SNF4_matrix <- scale_Pol[,intersect(colnames(scale_Pol),SNF4_ID)]

comparison_matrix <- matrix(ncol=4,nrow=nrow(scale_Pol))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4")
rownames(comparison_matrix) <- rownames(scale_Pol)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_matrix[i,]),na.rm = T)
}
comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
scale_Pol_SNF<-comparison_matrix

#lip
SNF1_matrix <- scale_lip_cat[,intersect(colnames(scale_lip_cat),SNF1_ID)]
SNF2_matrix <- scale_lip_cat[,intersect(colnames(scale_lip_cat),SNF2_ID)]
SNF3_matrix <- scale_lip_cat[,intersect(colnames(scale_lip_cat),SNF3_ID)]
SNF4_matrix <- scale_lip_cat[,intersect(colnames(scale_lip_cat),SNF4_ID)]

comparison_matrix <- matrix(ncol=4,nrow=nrow(scale_lip_cat))
colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4")
rownames(comparison_matrix) <- rownames(scale_lip_cat)

for (i in rownames(comparison_matrix)){
  comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_matrix[i,]),na.rm = T)
  comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_matrix[i,]),na.rm = T)
}
comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
scale_lip_cat_SNF<-comparison_matrix


# node_scale_PAM50<-rbind(scale_pro_PAM50,scale_pol_PAM50)


#####################################################################################################################
## 3. Metabolic protein and metabolite correlation network construction
#####################################################################################################################

# It takes several time to finish part 3/4 and their figures are not drawn by R....... temporarily skip them 

## Metabolic protein correlation network construction
pro_Cor.Res <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
pro_P.val <- matrix(nrow=nrow(TT_PRO),ncol=nrow(TT_PRO))
colnames(pro_Cor.Res) <- rownames(TT_PRO)
rownames(pro_Cor.Res) <- rownames(TT_PRO)
colnames(pro_P.val) <- rownames(TT_PRO)
rownames(pro_P.val) <- rownames(TT_PRO)

TT_PRO <- t(TT_PRO)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_PRO)){
  for (j in 1:ncol(TT_PRO)){
    TEMP_Res <- cor.test(TT_PRO[,i],TT_PRO[,j], method = "spearman")
    pro_Cor.Res[i,j] <- TEMP_Res$estimate
    pro_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_PRO))
}

pro_FDR <- matrix(p.adjust(pro_P.val,method="fdr"),ncol=2252)
colnames(pro_FDR) <- colnames(pro_P.val)
rownames(pro_FDR) <- rownames(pro_P.val)

## rna : Metabolic transcriptomics correlation network construction
TT_rna <- TT_rna[rownames(TT_rna)%in%rownames(luminal_rna_SNF),]
rna_Cor.Res <- matrix(nrow=nrow(TT_rna),ncol=nrow(TT_rna))
rna_P.val <- matrix(nrow=nrow(TT_rna),ncol=nrow(TT_rna))
colnames(rna_Cor.Res) <- rownames(TT_rna)
rownames(rna_Cor.Res) <- rownames(TT_rna)
colnames(rna_P.val) <- rownames(TT_rna)
rownames(rna_P.val) <- rownames(TT_rna)

TT_rna <- t(TT_rna)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_rna)){
  for (j in 1:ncol(TT_rna)){
    TEMP_Res <- cor.test(TT_rna[,i],TT_rna[,j], method = "spearman")
    rna_Cor.Res[i,j] <- TEMP_Res$estimate
    rna_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_rna))
}

rna_FDR <- matrix(p.adjust(rna_P.val,method="fdr"),ncol=836)
colnames(rna_FDR) <- colnames(rna_P.val)
rownames(rna_FDR) <- rownames(rna_P.val)

# save(rna_Cor.Res,rna_P.val,rna_FDR,file="rna_Cor_spearman.Rdata")

##lipid correlation network construction
lip_Cor.Res <- matrix(nrow=nrow(TT_lip_cat),ncol=nrow(TT_lip_cat))
lip_P.val <- matrix(nrow=nrow(TT_lip_cat),ncol=nrow(TT_lip_cat))
colnames(lip_Cor.Res) <- rownames(TT_lip_cat)
rownames(lip_Cor.Res) <- rownames(TT_lip_cat)
colnames(lip_P.val) <- rownames(TT_lip_cat)
rownames(lip_P.val) <- rownames(TT_lip_cat)

TT_lip_cat <- t(TT_lip_cat)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_lip_cat)){
  for (j in 1:ncol(TT_lip_cat)){
    TEMP_Res <- cor.test(TT_lip_cat[,i],TT_lip_cat[,j], method = "spearman")
    lip_Cor.Res[i,j] <- TEMP_Res$estimate
    lip_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_lip_cat))
}

lip_FDR <- matrix(p.adjust(lip_P.val,method="fdr"),ncol=46)
colnames(lip_FDR) <- colnames(lip_P.val)
rownames(lip_FDR) <- rownames(lip_P.val)

# save(lip_Cor.Res,lip_P.val,lip_FDR,file="lip_Cor_spearman.Rdata")


## Polar metabolite correlation network construction
pol_Cor.Res <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
pol_P.val <- matrix(nrow=nrow(TT_pol),ncol=nrow(TT_pol))
colnames(pol_Cor.Res) <- rownames(TT_pol)
rownames(pol_Cor.Res) <- rownames(TT_pol)
colnames(pol_P.val) <- rownames(TT_pol)
rownames(pol_P.val) <- rownames(TT_pol)

TT_pol <- t(TT_pol)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(TT_pol)){
  for (j in 1:ncol(TT_pol)){
    TEMP_Res <- cor.test(TT_pol[,i],TT_pol[,j], method = "spearman")
    pol_Cor.Res[i,j] <- TEMP_Res$estimate
    pol_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(TT_pol))
}

pol_FDR <- matrix(p.adjust(pol_P.val,method="fdr"),ncol=669)
colnames(pol_FDR) <- colnames(pol_P.val)
rownames(pol_FDR) <- rownames(pol_P.val)

# save(pol_Cor.Res,pol_P.val,pol_FDR,file="pol_Cor_spearman.Rdata")


## Metabolic protein and Polar metabolite correlation network construction
Luminal_PRO_POL_ID <- intersect(rownames(TT_PRO),rownames(TT_pol)) # 179

pro_pol_matrix <- cbind(TT_PRO[Luminal_PRO_POL_ID,],TT_pol[Luminal_PRO_POL_ID,]) #179*2921
dim(pro_pol_matrix)

pro_pol_Cor.Res <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
pro_pol_P.val <- matrix(nrow=ncol(pro_pol_matrix),ncol=ncol(pro_pol_matrix))
colnames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
rownames(pro_pol_Cor.Res) <- colnames(pro_pol_matrix)
colnames(pro_pol_P.val) <- colnames(pro_pol_matrix)
rownames(pro_pol_P.val) <- colnames(pro_pol_matrix)

pb <- txtProgressBar(style=3)
for (i in 1:ncol(pro_pol_matrix)){
  for (j in 1:ncol(pro_pol_matrix)){
    TEMP_Res <- cor.test(pro_pol_matrix[,i],pro_pol_matrix[,j], method = "spearman")
    pro_pol_Cor.Res[i,j] <- TEMP_Res$estimate
    pro_pol_P.val[i,j] <- TEMP_Res$p.value
  }
  setTxtProgressBar(pb, i/ncol(pro_pol_matrix))
}

pro_pol_FDR <- matrix(p.adjust(pro_pol_P.val,method="fdr"),ncol=2921)
colnames(pro_pol_FDR) <- colnames(pro_pol_P.val)
rownames(pro_pol_FDR) <- rownames(pro_pol_P.val)

save(pro_pol_Cor.Res,pro_pol_P.val,pro_pol_FDR,file="pro_pol_Cor_spearman.Rdata")
write.csv(pro_pol_Cor.Res,"pro_pol_Cor.Res.csv")
pro_pol_Cor.Res["CAD","M133T340_POS"]
pro_pol_P.val["CAD","M133T340_POS"]
pro_pol_FDR["CAD","M133T340_POS"]

##### [edge] #####
##### Metabolic protein [edge] #####
res.matrix_pro<-data.frame("a","b",0.5,0.05)
colnames(res.matrix_pro)<-c("source","target","weight","FDR")

for (i in 1:2252){
  for (j in (i+1):2252){
    if ((j!=i) & (pro_pol_FDR[i,j]<0.05) & (pro_pol_Cor.Res[i,j] > 0.4)){
      res.matrix_pro<-rbind(res.matrix_pro,c(rownames(pro_pol_Cor.Res)[i],colnames(pro_pol_Cor.Res)[j],
                                             pro_pol_Cor.Res[i,j],pro_pol_FDR[i,j]))
    }
  }
  setTxtProgressBar(pb, i/2252)
}
#19554
res.matrix_pro<-res.matrix_pro[-1,]
res.matrix_pro$type<-rep("undirected",n=nrow(res.matrix_pro))

#Deleting the protein without annotation
anno_manual_pro<-read.csv("anno_manual_pro.csv")
rownames(anno_manual_pro)<-anno_manual_pro$X
anno_manual_pro<-anno_manual_pro[,-1]

res.matrix_pro$anno_source<-anno_manual_pro[res.matrix_pro$source,"recon_kegg"]
res.matrix_pro$anno_target<-anno_manual_pro[res.matrix_pro$target,"recon_kegg"]
res.matrix_pro$anno_source[which(res.matrix_pro$anno_source=="")]<-NA
res.matrix_pro$anno_target[which(res.matrix_pro$anno_target=="")]<-NA
res.matrix_pro_final<-na.omit(res.matrix_pro) #4541
write.csv(res.matrix_pro_final[,c(1,2,3,5)],"res.matrix_pro.csv",row.names=FALSE)

##### rna [edge] #####
res.matrix_rna<-data.frame("a","b",0.5,0.05)
colnames(res.matrix_rna)<-c("source","target","weight","FDR")

load("rna_Cor_spearman.Rdata")
rna_FDR[is.na(rna_FDR)] <- 1
rna_Cor.Res[is.na(rna_Cor.Res)] <- 0
for (i in 1:836){
  for (j in (i+1):836){
    if ((j!=i) & (rna_FDR[i,j]<0.05) & (rna_Cor.Res[i,j] > 0.4)){
      res.matrix_rna<-rbind(res.matrix_rna,c(rownames(rna_Cor.Res)[i],colnames(rna_Cor.Res)[j],
                                             rna_Cor.Res[i,j],rna_FDR[i,j]))
    }
  }
  setTxtProgressBar(pb, i/836)
}
#4645
res.matrix_rna<-res.matrix_rna[-1,]
res.matrix_rna$type<-rep("undirected",n=nrow(res.matrix_rna))
write.csv(res.matrix_rna[,c(1,2,3,5)],"res.matrix_rna.csv",row.names=FALSE)

##### lip [edge] #####
res.matrix_lip<-data.frame("a","b",0.5,0.05)
colnames(res.matrix_lip)<-c("source","target","weight","FDR")

for (i in 1:46){
  for (j in (i+1):46){
    if ((j!=i) & (lip_FDR[i,j]<0.05) & (lip_Cor.Res[i,j] > 0.4)){
      res.matrix_lip<-rbind(res.matrix_lip,c(rownames(lip_Cor.Res)[i],colnames(lip_Cor.Res)[j],
                                             lip_Cor.Res[i,j],lip_FDR[i,j]))
    }
  }
  setTxtProgressBar(pb, i/46)
}
#697
res.matrix_lip<-res.matrix_lip[-1,]
res.matrix_lip$type<-rep("undirected",n=nrow(res.matrix_lip))
write.csv(res.matrix_lip[,c(1,2,3,5)],"res.matrix_lip.csv",row.names=FALSE)

##### Polar Metabolite [edge] #####
res.matrix_pol<-data.frame("a","b",0.5,0.05)
colnames(res.matrix_pol)<-c("source","target","weight","FDR")

for (i in 2253:2921){
  for (j in (i+1):2921){
    if ((j!=i) & (pro_pol_FDR[i,j]<0.05) & (pro_pol_Cor.Res[i,j] > 0.4)){
      res.matrix_pol<-rbind(res.matrix_pol,c(rownames(pro_pol_Cor.Res)[i],colnames(pro_pol_Cor.Res)[j],
                                             pro_pol_Cor.Res[i,j],pro_pol_FDR[i,j]))
    }
  }
  setTxtProgressBar(pb, i/669)
}
#17467
res.matrix_pol<-res.matrix_pol[-1,]
res.matrix_pol$type<-rep("undirected",n=nrow(res.matrix_pol))
write.csv(res.matrix_pol[,c(1,2,3,5)],"res.matrix_pol.csv",row.names=FALSE)


##### [node] #####
# scale_pro_PAM50<-read.table("scale_pro_PAM50.txt",sep="\t")
# scale_pol_PAM50<-read.table("scale_pol_PAM50.txt",sep="\t")

scale_node<-rbind(scale_pro_SNF,scale_Pol_SNF)

scale_node <- scale_lip_cat_SNF
scale_rna_SNF <- na.omit(scale_rna_SNF)
scale_node <- scale_rna_SNF


max(scale_node)
min(scale_node)
for(i in 1:ncol(scale_node)){
  scale_node[,i][which(scale_node[,i] >= 1)]<- 1
  scale_node[,i][which(scale_node[,i] >= 0.8 & scale_node[,i] < 1)]<- 0.8
  scale_node[,i][which(scale_node[,i] >= 0.6 & scale_node[,i] < 0.8)]<- 0.6
  scale_node[,i][which(scale_node[,i] >= 0.4 & scale_node[,i] < 0.6)]<- 0.4
  scale_node[,i][which(scale_node[,i] >= 0.2 & scale_node[,i] < 0.4)]<- 0.2
  scale_node[,i][which(scale_node[,i] >= 0 & scale_node[,i] < 0.2)]<- 0
  scale_node[,i][which(scale_node[,i] >= -0.2 & scale_node[,i] < 0)]<- -0.2
  scale_node[,i][which(scale_node[,i] >= -0.4 & scale_node[,i] < -0.2)]<- -0.4
  scale_node[,i][which(scale_node[,i] >= -0.6 & scale_node[,i] < -0.4)]<- -0.6
  scale_node[,i][which(scale_node[,i] >= -0.8 & scale_node[,i] < -0.6)]<- -0.8
  scale_node[,i][which(scale_node[,i] < -0.8)]<- -1
}


anno_manual<-read.csv("anno_manual_pro.csv")
rownames(anno_manual)<-anno_manual$id
anno_manual<-anno_manual[,-1]
node_SNF<-luminal_rna_SNF

res.matrix_cus<-res.matrix_pol
res.matrix_cus<-res.matrix_pro_final
res.matrix_cus<-res.matrix_lip
res.matrix_cus<-res.matrix_rna

node_cor<-union(unique(res.matrix_cus$source),unique(res.matrix_cus$target))
node_anno<-data.frame(node_cor,"type","Label","pathway_class",
                      "SNF1","SNF2","SNF3","SNF4","max_min",
                      "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
colnames(node_anno)<-c("Id","type","Label","pathway_class",
                       "SNF1","SNF2","SNF3","SNF4","max_min",
                       "SD","P","FDR","Tumor","PT","T_PT","P_TN","FDR_TN")
rownames(node_anno)<-node_anno$Id

for (i in rownames(node_anno)){
  if (i %in% rownames(CBCGA_pol_anno)){
    node_anno[i,"type"]<-"metabolite"
    node_anno[i,"Label"]<-CBCGA_pol_anno[i,"Putative_metabolite_name"]
    node_anno[i,"pathway_class"]<-anno_manual[i,"subclass_final_manual"]
    node_anno[i,"subclass"]<-CBCGA_pol_anno[i,"Metabolite_class"]
  } else {
    node_anno[i,"type"]<-"protein"
    node_anno[i,"Label"]<-node_anno[i,"Id"]
    node_anno[i,"pathway_class"]<-anno_manual_pro[i,"recon_kegg"]}
  node_anno[i,5:8]<-scale_node[i,]
  node_anno[i,9:17]<-node_SNF[i,5:13]
}
node_anno_pol<-node_anno
node_anno_pro<-node_anno
node_anno_lip<-node_anno
node_anno_rna<-node_anno

write.csv(node_anno,"node_anno_pro.csv",row.names = F)
write.csv(node_anno,"node_anno_pol.csv",row.names = F)
write.csv(node_anno,"node_anno_lip.csv",row.names = F)
write.csv(node_anno,"node_anno_rna.csv",row.names = F)

## Refer to Gephi

#################lipboxplot########################
if (F) {
  lipid_main <- unique(CBCGA_lip_anno$Lipid.super.class)
  
  Cus_lipid_main_TT <- matrix(ncol=ncol(TT_lip),nrow=length(lipid_main))
  rownames(Cus_lipid_main_TT) <- lipid_main
  colnames(Cus_lipid_main_TT) <- colnames(TT_lip)
  
  for (i in rownames(Cus_lipid_main_TT)){
    peak <- rownames(CBCGA_lip_anno)[CBCGA_lip_anno$Lipid.super.class==i]
    Cus_lipid_main_TT[i,] <- c(apply(TT_lip[peak,],2,mean))
  }
  
  SNF1_matrix <- Cus_lipid_main_TT[,intersect(colnames(Cus_lipid_main_TT),SNF1_ID)]
  SNF2_matrix <- Cus_lipid_main_TT[,intersect(colnames(Cus_lipid_main_TT),SNF2_ID)]
  SNF3_matrix <- Cus_lipid_main_TT[,intersect(colnames(Cus_lipid_main_TT),SNF3_ID)]
  SNF4_matrix <- Cus_lipid_main_TT[,intersect(colnames(Cus_lipid_main_TT),SNF4_ID)]
  
  comparison_matrix <- matrix(ncol=4,nrow=nrow(Cus_lipid_main_TT))
  colnames(comparison_matrix) <- c("SNF1","SNF2","SNF3","SNF4")
  rownames(comparison_matrix) <- rownames(Cus_lipid_main_TT)
  
  for (i in rownames(comparison_matrix)){
    comparison_matrix[i,"SNF1"] <- mean(as.numeric(SNF1_matrix[i,]),na.rm = T)
    comparison_matrix[i,"SNF2"] <- mean(as.numeric(SNF2_matrix[i,]),na.rm = T)
    comparison_matrix[i,"SNF3"] <- mean(as.numeric(SNF3_matrix[i,]),na.rm = T)
    comparison_matrix[i,"SNF4"] <- mean(as.numeric(SNF4_matrix[i,]),na.rm = T)
  }
  comparison_matrix<-as.data.frame(comparison_matrix,stringsAsFactors = F)
  Cus_lipid_main_TT_SNF<-comparison_matrix
}
library(tidyr)
library(pheatmap)
library(RColorBrewer)

luminal_lip_SNF <- read.csv("./results/luminal_lip_SNF.csv",header = T,row.names = 1)
luminal_lip_SNF$main <- CBCGA_lip_anno$Lipid.super.class
luminal_lip_SNF$SNF1_PT <- luminal_lip_SNF$SNF1-luminal_lip_SNF$PT
luminal_lip_SNF$SNF2_PT <- luminal_lip_SNF$SNF2-luminal_lip_SNF$PT
luminal_lip_SNF$SNF3_PT <- luminal_lip_SNF$SNF3-luminal_lip_SNF$PT
luminal_lip_SNF$SNF4_PT <- luminal_lip_SNF$SNF4-luminal_lip_SNF$PT

plotdata <- luminal_lip_SNF[,c(14:18)]
plotdata$peak <- rownames(plotdata)

plotdata_FA <- subset(plotdata,plotdata$main=="Fatty acyls [FA]")
kruskal.test(list(plotdata_FA$SNF1_PT,
                  plotdata_FA$SNF2_PT,
                  plotdata_FA$SNF3_PT,
                  plotdata_FA$SNF4_PT))#p-value = 0.0004746***
plotdata_GL <- subset(plotdata,plotdata$main=="Glycerolipids [GL]")
kruskal.test(list(plotdata_GL$SNF1_PT,
                  plotdata_GL$SNF2_PT,
                  plotdata_GL$SNF3_PT,
                  plotdata_GL$SNF4_PT))#p-value = 0.008005**
plotdata_GP <- subset(plotdata,plotdata$main=="Glycerophospholipids [GP]")
kruskal.test(list(plotdata_GP$SNF1_PT,
                  plotdata_GP$SNF2_PT,
                  plotdata_GP$SNF3_PT,
                  plotdata_GP$SNF4_PT))#p-value = 9.113e-14***
plotdata_SP <- subset(plotdata,plotdata$main=="Sphingolipids [SP]")
kruskal.test(list(plotdata_SP$SNF1_PT,
                  plotdata_SP$SNF2_PT,
                  plotdata_SP$SNF3_PT,
                  plotdata_SP$SNF4_PT))#0.002252**
plotdata_ST <- subset(plotdata,plotdata$main=="Sterol Lipids [ST]")
kruskal.test(list(plotdata_ST$SNF1_PT,
                  plotdata_ST$SNF2_PT,
                  plotdata_ST$SNF3_PT,
                  plotdata_ST$SNF4_PT))#p-value = 0.761



plotdata <- gather(plotdata,SNF,abundance,-main,-peak)

library(ggplot2)

#please select optimal color and xlim, ylim

ggplot(plotdata,aes(x = main, y = abundance, fill = SNF)) +
  scale_fill_manual(values=c(SNF1_PT="#3D76AE",SNF2_PT="#53AD4A",SNF3_PT="#EDAB3C",SNF4_PT="#C3392A") )+
  geom_boxplot(outlier.shape = NA,linetype="dashed")+
  stat_boxplot(outlier.shape=NA,aes(ymin=..lower..,ymax=..upper..))+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.75)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.75)+ 
  labs(x="lipid abundance",y="Mean |log2FC| between tumor and normal")+
  # ggtitle("Boxplot of hub gene") +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"))
#boxplot of mean log2FC of lipid categories
boxplot(log2FC~catagory,data=Cus_lipid_main_TT_SNF,
        at=c(1:5),xlim=c(0,6),ylim=c(-2,2),
        col= c("#3D76AE" ,"#53AD4A", "#EDAB3C", "#C3392A"), 
        #names=names(table(comparison_TT_PT_matrix_cat$catagory_whole)),
        outline=F,
        main="Fold change"
)
abline(h=0,col = "red", lwd = 2, lty = 2)

#############polar boxplot################

luminal_Pol_SNF <- read.csv("./results/luminal_Pol_SNF.csv",header = T,row.names = 1)
luminal_Pol_SNF$main <- CBCGA_pol_anno$Metabolite_class
luminal_Pol_SNF$SNF1_PT <- luminal_Pol_SNF$SNF1-luminal_Pol_SNF$PT
luminal_Pol_SNF$SNF2_PT <- luminal_Pol_SNF$SNF2-luminal_Pol_SNF$PT
luminal_Pol_SNF$SNF3_PT <- luminal_Pol_SNF$SNF3-luminal_Pol_SNF$PT
luminal_Pol_SNF$SNF4_PT <- luminal_Pol_SNF$SNF4-luminal_Pol_SNF$PT

plotdata <- luminal_Pol_SNF[,c(14:18)]
plotdata$peak <- rownames(plotdata)

plotdata_AA <- subset(plotdata,plotdata$main=="Amino acid")
kruskal.test(list(plotdata_AA$SNF1_PT,
                  plotdata_AA$SNF2_PT,
                  plotdata_AA$SNF3_PT,
                  plotdata_AA$SNF4_PT))#p-value = 4.952e-05***
plotdata_carbon <- subset(plotdata,plotdata$main=="Carbohydrates")
kruskal.test(list(plotdata_carbon$SNF1_PT,
                  plotdata_carbon$SNF2_PT,
                  plotdata_carbon$SNF3_PT,
                  plotdata_carbon$SNF4_PT))#p-value = 0.4084
plotdata_lip <- subset(plotdata,plotdata$main=="Lipid")
kruskal.test(list(plotdata_lip$SNF1_PT,
                  plotdata_lip$SNF2_PT,
                  plotdata_lip$SNF3_PT,
                  plotdata_lip$SNF4_PT))#p-value = 0.02426*
plotdata_nucle <- subset(plotdata,plotdata$main=="Nucleotide")
kruskal.test(list(plotdata_nucle$SNF1_PT,
                  plotdata_nucle$SNF2_PT,
                  plotdata_nucle$SNF3_PT,
                  plotdata_nucle$SNF4_PT))#0.3861
plotdata_Other <- subset(plotdata,plotdata$main=="Other")
kruskal.test(list(plotdata_Other$SNF1_PT,
                  plotdata_Other$SNF2_PT,
                  plotdata_Other$SNF3_PT,
                  plotdata_Other$SNF4_PT))#p-value = 0.6101
plotdata_pep <- subset(plotdata,plotdata$main=="Peptide")
kruskal.test(list(plotdata_pep$SNF1_PT,
                  plotdata_pep$SNF2_PT,
                  plotdata_pep$SNF3_PT,
                  plotdata_pep$SNF4_PT))#p-value = 0.0002545***
plotdata_vita <- subset(plotdata,plotdata$main=="Vitamins and Cofactors")
kruskal.test(list(plotdata_vita$SNF1_PT,
                  plotdata_vita$SNF2_PT,
                  plotdata_vita$SNF3_PT,
                  plotdata_vita$SNF4_PT))#p-value = 0.8432
plotdata_xeno <- subset(plotdata,plotdata$main=="Xenobiotics")
kruskal.test(list(plotdata_xeno$SNF1_PT,
                  plotdata_xeno$SNF2_PT,
                  plotdata_xeno$SNF3_PT,
                  plotdata_xeno$SNF4_PT))#p-value = 0.7095


plotdata <- gather(plotdata,SNF,abundance,-main,-peak)

library(ggplot2)

#please select optimal color and xlim, ylim

ggplot(plotdata,aes(x = main, y = abundance, fill = SNF)) +
  scale_fill_manual(values=c(SNF1_PT="#3D76AE",SNF2_PT="#53AD4A",SNF3_PT="#EDAB3C",SNF4_PT="#C3392A") )+
  geom_boxplot(outlier.shape = NA,linetype="dashed")+
  stat_boxplot(outlier.shape=NA,aes(ymin=..lower..,ymax=..upper..))+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.75)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.75)+ 
  labs(x="lipid abundance",y="lipid abundance")+
  # ggtitle("Boxplot of hub gene") +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_blank())+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid")) 

############## gene heatmap###############
#"node_anno_rna":SNF1~3——SNF>0.6; SNF4——SNF<-0.8
rnalist <- read.csv("figs2_rnalist.csv")

rnaheat <- luminal_rna_SNF[rnalist$genelist,1:4]

library(wesanderson)
n <- rnaheat
ac=data.frame(group=c("SNF1","SNF2","SNF3","SNF4"))
rownames(ac)=colnames(n) 
P <- pheatmap(as.matrix(n),show_colnames =F,
              show_rownames = T,
              scale = "row",
              cluster_rows = F, cluster_cols = F,
              # main = "MGY vs MGN",
              annotation_col=ac,
              annotation_names_col =F,
              annotation_colors = list(group=c(SNF1="#3D76AE",SNF2="#53AD4A",SNF3="#EDAB3C",SNF4="#C3392A")),
              color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
              # legend = T,legend_labels=c("Activity"),
              # gaps_row=c(9,12)                  
) 
P


SNF_gene <- as.character(row.names(exp.fpkm.TT))
KEGG_gene <- as.character(unique(KEGG_Meta_signature$Genes))
othergene <- setdiff(SNF_gene, KEGG_gene)
inter_gene <- intersect(SNF_gene, KEGG_gene)

exp.fpkm.TT_kegg <- exp.fpkm.TT[inter_gene,]
exp.fpkm.TT_normal <- exp.fpkm.PT[inter_gene,]

exp.fpkm.TT_kegg <- log2(cbind(exp.fpkm.TT_kegg,exp.fpkm.TT_normal)+1)


# test <- (exp.fpkm.TT_normal[2,1]-exp.fpkm.TT_normal[2,2])^2


## Calculate the Euclidean distance between Tumor and Normal
RMSD <- rep(0, 3861)
sum <- 0
x <- 0

for (i in 1:351)
{
  for (j in 352:362)
  {
    x <- x + 1
    for (k in 1:nrow(exp.fpkm.TT_kegg))
    {
      sum <- sum + (exp.fpkm.TT_kegg[k,i]-exp.fpkm.TT_kegg[k,j])^2
    }
    RMSD[x] <- sqrt(sum/nrow(exp.fpkm.TT_kegg))
    sum <- 0
  }
}

write.csv(RMSD,"Metabolic_gene_RMSD_Tumor_Normal.csv")

## Calculate the Euclidean distance between Normal and Normal
RMSD <- rep(0, 55)
sum <- 0
x <- 0

for (i in 352:361)
{
  for (j in (i+1):362)
  {
    x <- x + 1
    for (k in 1:nrow(exp.fpkm.TT_kegg))
    {
      sum <- sum + (exp.fpkm.TT_kegg[k,i]-exp.fpkm.TT_kegg[k,j])^2
    }
    RMSD[x] <- sqrt(sum/nrow(exp.fpkm.TT_kegg))
    sum <- 0
  }
}

write.csv(RMSD,"Metabolic_gene_RMSD_Normal_Normal.csv")

## Calculate the Euclidean distance between Tumor and Tumor
RMSD <- rep(0, 61425)
sum <- 0
x <- 0

for (i in 1:350)
{
  for (j in (i+1):351)
  {
    x <- x + 1
    for (k in 1:nrow(exp.fpkm.TT_kegg))
    {
      sum <- sum + (exp.fpkm.TT_kegg[k,i]-exp.fpkm.TT_kegg[k,j])^2
    }
    RMSD[x] <- sqrt(sum/nrow(exp.fpkm.TT_kegg))
    sum <- 0
  }
}

write.csv(RMSD,"Metabolic_gene_RMSD_Tumor_Tumor.csv")


#################plot########################

group <- c("T VS N","T VS T","N VS N")

TVST <- read.csv("Metabolic_gene_RMSD_Tumor_Tumor.csv")
TVST <- as.data.frame(TVST[,2])
TVST$group <- "T VS T"
colnames(TVST) <- c("Euclidean expression distances","group")

TVSN <- read.csv("Metabolic_gene_RMSD_Tumor_Normal.csv")
TVSN <- as.data.frame(TVSN[,2])
TVSN$group <- "T VS N"
colnames(TVSN) <- c("Euclidean expression distances","group")

NVSN <- read.csv("Metabolic_gene_RMSD_Normal_Normal.csv")
NVSN <- as.data.frame(NVSN[,2])
NVSN$group <- "N VS N"
colnames(NVSN) <- c("Euclidean expression distances","group")

plotdata <- rbind(TVST,TVSN,NVSN)
plotdata[,1] <- as.numeric(plotdata[,1])
plotdata[,2] <- factor(plotdata[,2],levels = c("T VS N","T VS T","N VS N"))

########wilcox##########
wilcox.test(TVST$`Euclidean expression distances`,TVSN$`Euclidean expression distances`)
wilcox.test(TVST$`Euclidean expression distances`,NVSN$`Euclidean expression distances`)
wilcox.test(TVSN$`Euclidean expression distances`,NVSN$`Euclidean expression distances`)

mean(TVST$`Euclidean expression distances`)
mean(TVSN$`Euclidean expression distances`)
mean(NVSN$`Euclidean expression distances`)


library(ggplot2)
library(RColorBrewer)
ggplot(plotdata,aes(x = group, y =`Euclidean expression distances`, fill = group)) +
  geom_boxplot(width=0.2,outlier.shape = NA,linetype="dashed")+
  stat_boxplot(outlier.shape=NA,aes(ymin=..lower..,ymax=..upper..))+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.3)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.3)+ 
  scale_fill_manual(values = colorRampPalette(brewer.pal(6, "Spectral"))(3))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.6))+
  labs(x="SNF subtype",y="Distribution distance (r.m.s.d.)")+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid")) 


##############################################################################
{
  ## Calculate the Euclidean distance of all genes
  RMSD <- rep(0, 8680)
  sum <- 0
  x <- 0
  
  for (i in 1:360)
  {
    for (j in 361:448)
    {
      x <- x + 1
      for (k in 1:nrow(FUSCCTNBC_FPKM_log_trans))
      {
        sum <- sum + (FUSCCTNBC_FPKM_log_trans[k,i]-FUSCCTNBC_FPKM_log_trans[k,j])^2
      }
      RMSD[x] <- sqrt(sum/nrow(FUSCCTNBC_FPKM_log_trans))
      sum <- 0
    }
  }
  
  write.csv(RMSD,"All_gene_RMSD_Tumor_Normal.csv")
  
  RMSD <- rep(0, 3828)
  sum <- 0
  x <- 0
  
  for (i in 361:447)
  {
    for (j in (i+1):448)
    {
      x <- x + 1
      for (k in 1:nrow(FUSCCTNBC_FPKM_log_trans))
      {
        sum <- sum + (FUSCCTNBC_FPKM_log_trans[k,i]-FUSCCTNBC_FPKM_log_trans[k,j])^2
      }
      RMSD[x] <- sqrt(sum/nrow(FUSCCTNBC_FPKM_log_trans))
      sum <- 0
    }
  }
  
  write.csv(RMSD,"All_gene_RMSD_Normal_Normal.csv")
  
  RMSD <- rep(0, 64620)
  sum <- 0
  x <- 0
  
  for (i in 1:359)
  {
    for (j in (i+1):360)
    {
      x <- x + 1
      for (k in 1:nrow(FUSCCTNBC_FPKM_log_trans))
      {
        sum <- sum + (FUSCCTNBC_FPKM_log_trans[k,i]-FUSCCTNBC_FPKM_log_trans[k,j])^2
      }
      RMSD[x] <- sqrt(sum/nrow(FUSCCTNBC_FPKM_log_trans))
      sum <- 0
    }
  }
  
  write.csv(RMSD,"All_gene_RMSD_Tumor_Tumor.csv")
}

##############SNF##################
# Calculate Euclidean distance for BCMA expression data
exp.fpkm.TT_kegg

BCMA_RNA_Combat_log2_FPKM_PT <- exp.fpkm.TT_kegg[,352:362] #11
BCMA_RNA_Combat_log2_FPKM_SNF1 <- exp.fpkm.TT_kegg[,intersect(colnames(exp.fpkm.TT_kegg),SNF1_ID)] # 86
BCMA_RNA_Combat_log2_FPKM_SNF2 <- exp.fpkm.TT_kegg[,intersect(colnames(exp.fpkm.TT_kegg),SNF2_ID)] # 89
BCMA_RNA_Combat_log2_FPKM_SNF3 <- exp.fpkm.TT_kegg[,intersect(colnames(exp.fpkm.TT_kegg),SNF3_ID)] #118
BCMA_RNA_Combat_log2_FPKM_SNF4 <- exp.fpkm.TT_kegg[,intersect(colnames(exp.fpkm.TT_kegg),SNF4_ID)] #58

## Calculate the Euclidean distance between Tumor and Normal

comparison_Lum_SNF2E_SNF1_distance <- matrix(ncol=4,nrow=nrow(BCMA_RNA_Combat_log2_FPKM_SNF1))
colnames(comparison_Lum_SNF2E_SNF1_distance) <- c("SNF3","SNF4","SNF2","SNF1")
rownames(comparison_Lum_SNF2E_SNF1_distance) <- rownames(BCMA_RNA_Combat_log2_FPKM_SNF1)

## SNF1 vs PT
RMSD <- rep(0, ncol(BCMA_RNA_Combat_log2_FPKM_SNF1)*11)
sum <- 0
x <- 0

for (i in 1:ncol(BCMA_RNA_Combat_log2_FPKM_SNF1)){
  for (j in 1:11){
    x <- x + 1
    for (k in 1:nrow(comparison_Lum_SNF2E_SNF1_distance)){
      sum <- sum + (BCMA_RNA_Combat_log2_FPKM_SNF1[k,i]-BCMA_RNA_Combat_log2_FPKM_PT[k,j])^2}
    RMSD[x] <- sqrt(sum/nrow(comparison_Lum_SNF2E_SNF1_distance))
    sum <- 0}
  print(i)}
RMSD_SNF1_PT <- c(RMSD)
length(RMSD_SNF1_PT) <- 946


## SNF2 vs PT
RMSD <- rep(0, ncol(BCMA_RNA_Combat_log2_FPKM_SNF2)*11)
sum <- 0
x <- 0

for (i in 1:ncol(BCMA_RNA_Combat_log2_FPKM_SNF2)){
  for (j in 1:11){
    x <- x + 1
    for (k in 1:nrow(comparison_Lum_SNF2E_SNF1_distance)){
      sum <- sum + (BCMA_RNA_Combat_log2_FPKM_SNF2[k,i]-BCMA_RNA_Combat_log2_FPKM_PT[k,j])^2}
    RMSD[x] <- sqrt(sum/nrow(comparison_Lum_SNF2E_SNF1_distance))
    sum <- 0}
  print(i)}
RMSD_SNF2_PT <- c(RMSD)
length(RMSD_SNF2_PT) <- 979


## SNF3 vs PT
RMSD <- rep(0, ncol(BCMA_RNA_Combat_log2_FPKM_SNF3)*11)
sum <- 0
x <- 0

for (i in 1:ncol(BCMA_RNA_Combat_log2_FPKM_SNF3)){
  for (j in 1:11){
    x <- x + 1
    for (k in 1:nrow(comparison_Lum_SNF2E_SNF1_distance)){
      sum <- sum + (BCMA_RNA_Combat_log2_FPKM_SNF3[k,i]-BCMA_RNA_Combat_log2_FPKM_PT[k,j])^2}
    RMSD[x] <- sqrt(sum/nrow(comparison_Lum_SNF2E_SNF1_distance))
    sum <- 0}
  print(i)}
RMSD_SNF3_PT <- c(RMSD)
length(RMSD_SNF3_PT) <- 1298


## SNF4 vs PT
RMSD <- rep(0, ncol(BCMA_RNA_Combat_log2_FPKM_SNF4)*11)
sum <- 0
x <- 0

for (i in 1:ncol(BCMA_RNA_Combat_log2_FPKM_SNF4)){
  for (j in 1:11){
    x <- x + 1
    for (k in 1:nrow(comparison_Lum_SNF2E_SNF1_distance)){
      sum <- sum + (BCMA_RNA_Combat_log2_FPKM_SNF4[k,i]-BCMA_RNA_Combat_log2_FPKM_PT[k,j])^2}
    RMSD[x] <- sqrt(sum/nrow(comparison_Lum_SNF2E_SNF1_distance))
    sum <- 0}
  print(i)}
RMSD_SNF4_PT <- c(RMSD)
length(RMSD_SNF4_PT) <- 638


#### GGPLOT2
SNF3_PT_all <- matrix(ncol=2,nrow=length(RMSD_SNF3_PT))
SNF3_PT_all[,1] <- "SNF3"
SNF3_PT_all[,2] <- RMSD_SNF3_PT

SNF4_PT_all <- matrix(ncol=2,nrow=length(RMSD_SNF4_PT))
SNF4_PT_all[,1] <- "SNF4"
SNF4_PT_all[,2] <- RMSD_SNF4_PT

SNF1_PT_all <- matrix(ncol=2,nrow=length(RMSD_SNF1_PT))
SNF1_PT_all[,1] <- "SNF1"
SNF1_PT_all[,2] <- RMSD_SNF1_PT

SNF2_PT_all <- matrix(ncol=2,nrow=length(RMSD_SNF2_PT))
SNF2_PT_all[,1] <- "SNF2"
SNF2_PT_all[,2] <- RMSD_SNF2_PT

comparison_SNF_BCMA_ggboxplot <- rbind(SNF3_PT_all,SNF4_PT_all,SNF2_PT_all,SNF1_PT_all)
colnames(comparison_SNF_BCMA_ggboxplot) <- c("SNF","RMSD")
comparison_SNF_BCMA_ggboxplot <- as.data.frame(comparison_SNF_BCMA_ggboxplot)
comparison_SNF_BCMA_ggboxplot[,2] <- as.numeric(comparison_SNF_BCMA_ggboxplot[,2])

#
comparison_SNF_BCMA_ggboxplot$SNF <- factor(comparison_SNF_BCMA_ggboxplot$SNF,levels = c("SNF1","SNF2","SNF3","SNF4"))

ggplot(comparison_SNF_BCMA_ggboxplot,
       aes(x=SNF,y=RMSD,fill=SNF))+
  #scale_x_discrete(limits=c("polar","lipid","merge"))+
  scale_fill_manual(values=c(SNF1="#3D76AE",SNF2="#53AD4A",SNF3="#EDAB3C",SNF4="#C3392A") )+
  geom_boxplot(outlier.shape = NA,linetype="dashed")+
  stat_boxplot(outlier.shape=NA,aes(ymin=..lower..,ymax=..upper..))+
  stat_boxplot(geom = "errorbar",aes(ymin=..ymax..),width=0.3)+
  stat_boxplot(geom = "errorbar",aes(ymax=..ymin..),width=0.3)+ 
  scale_y_continuous(expand = c(0,0.1),limits = c(0.3,1.6))+
  labs(x="SNF subtype",y="Distribution distance (r.m.s.d.)")+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background = element_blank())+ 
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"))
#stat_compare_means(hide.ns = TRUE,label = "p.signif",method="wilcox.test")
#geom_signif(comparisons = compaired,step_increase = 0.06,map_signif_level = T,test = wilcox.test) 

mean(as.numeric(SNF1_PT_all[,2]))#0.7590856
mean(as.numeric(SNF2_PT_all[,2]))#0.7675502
mean(as.numeric(SNF3_PT_all[,2]))#0.8377549
mean(as.numeric(SNF4_PT_all[,2]))#0.6759793
#1/2=99%
#1/3=91%
#1/4=89%
#2/3=92%
#2/4=88%
#3/4=81%

wilcox.test(as.numeric(SNF1_PT_all[,2]),as.numeric(SNF2_PT_all[,2]))
#p-value = 0.00412
wilcox.test(as.numeric(SNF1_PT_all[,2]),as.numeric(SNF3_PT_all[,2]))
#p-value < 2.2e-16
wilcox.test(as.numeric(SNF1_PT_all[,2]),as.numeric(SNF4_PT_all[,2]))
#p-value < 2.2e-16
wilcox.test(as.numeric(SNF2_PT_all[,2]),as.numeric(SNF3_PT_all[,2]))
#p-value < 2.2e-16
wilcox.test(as.numeric(SNF2_PT_all[,2]),as.numeric(SNF4_PT_all[,2]))
#p-value < 2.2e-16
wilcox.test(as.numeric(SNF3_PT_all[,2]),as.numeric(SNF4_PT_all[,2]))
#p-value < 2.2e-16

kruskal.test(list(as.numeric(SNF1_PT_all[,2]),
                  as.numeric(SNF2_PT_all[,2]),
                  as.numeric(SNF3_PT_all[,2]),
                  as.numeric(SNF4_PT_all[,2])))
#p-value < 2.2e-16

{
  ggplot(plotdata,aes(x = group, y =plotdata$`Euclidean expression distances`, fill = group)) +
    stat_boxplot(geom = "errorbar", width=0.6)+
    geom_boxplot(width = 0.6, outlier.shape=NA)+
    scale_fill_manual(values = colorRampPalette(brewer.pal(6, "Spectral"))(3))+
    scale_y_continuous(expand = c(0,0),limits = c(0,2))+
    # scale_x_discrete(name = "lipid category") +
    # ggtitle("Boxplot of hub gene") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face =  "bold"),
          text = element_text(size = 12),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 11),
          panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank())
}

mut_met <- read.csv("mut_met.csv")
pol_anno <- CBCGA_pol_anno[,c(1,5,6)]
mut_met <- merge(mut_met,pol_anno,by.x = "metabolite",by.y = "peak",all = F)
write.csv(mut_met,"mut_met_anno.csv")

mut_lip <- read.csv("mut_lip.csv")
pol_anno <- CBCGA_lip_anno[,c(2,4)]
mut_lip <- merge(mut_lip,pol_anno,by.x = "metabolite",by.y = "Subclass",all = F)
write.csv(mut_lip,"mut_lip_anno.csv")

rna_met <- read.csv("rna_met.csv")
pol_anno <- CBCGA_pol_anno[,c(1,5,6)]
rna_met <- merge(rna_met,pol_anno,by.x = "metabolite",by.y = "peak",all = F)
write.csv(rna_met,"rna_met_anno.csv")

cnv_met <- read.csv("cnv_met.csv")
pol_anno <- CBCGA_pol_anno[,c(1,5,6)]
cnv_met <- merge(cnv_met,pol_anno,by.x = "metabolite",by.y = "peak",all = F)
write.csv(cnv_met,"cnv_met_anno.csv")

cyc_met <- read.csv("cyc_met.csv")
pol_anno <- CBCGA_pol_anno[,c(1,5,6)]
cyc_met <- merge(cyc_met,pol_anno,by.x = "metabolite",by.y = "peak",all = F)
write.csv(cyc_met,"cyc_met_anno.csv")

luminal_Pol_SNF$metabolite <- rownames(luminal_Pol_SNF)
luminal_Pol_SNF <- merge(luminal_Pol_SNF,pol_anno,by.x = "metabolite",by.y = "peak",all = F)
write.csv(luminal_Pol_SNF,"luminal_Pol_SNF.csv")

#########################Plot############################################

library(ggplot2)
library(ggpubr)
library(gridExtra)

p1 <- ggplot(Customed_mutation_plus_metabolomics, aes(x=order,y=M613T483_POS, colour = mutation, shape = mutation))+
  geom_point(size = 3) + theme_bw() + ylab("Log2 levels") + scale_color_manual(values = c("black","red")) + scale_shape_manual(values = c(16,17)) +
  theme(panel.grid =element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position="none") 

p2 <- ggplot(Customed_mutation_plus_metabolomics, aes(x = order, y = 0)) + 
  geom_tile(aes(fill = mutation)) + theme_bw() + xlab("Ordered samples") + ylab("mutation mutations") +
  theme(panel.grid =element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position="none") +
  scale_fill_manual(values = c("white","black"))

gp1<- ggplot_gtable(ggplot_build(p1))
gp2<- ggplot_gtable(ggplot_build(p2))
#This identifies the maximum width
maxWidth = unit.pmin(gp1$widths[2:3], gp2$widths[2:3])
#Set each to the maximum width
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
#Put them together
grid.arrange(gp1, gp2)




