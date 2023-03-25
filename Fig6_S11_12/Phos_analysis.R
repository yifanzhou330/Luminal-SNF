rm(list = ls()) ; graphics.off()
setwd("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/result/inferSNF_res_EachSNFfixed_V2/CPTAC_2022/DEPhos")
library(clusterProfiler)
test_summary = openxlsx:::read.xlsx("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/31.NG二修/分型外推/SVM_RF/rawdata/inferSNF_res_EachSNFfixed_V2/RawNum100_CorThresh0.75_RF_scale_fold0_CPTAC.xlsx")
rownames(test_summary) = test_summary$X1
phos = data.table::fread("/Users/ZYF/Desktop/bioinformatics/bioinformatics_resource/外部组学队列/CPTAC/CPTAC 2020/S060_Breast_Cancer_data_freeze_GCTfiles_v5.4-public/prosp-brca-v5.4-public-phosphoproteome-ratio-norm-NArm.gct.csv",
                         data.table = F)
rownames(phos) = phos[,1]
phos_annot = phos[,c(1:22)]
phos = phos[,-c(1:22)]

tmp = intersect(rownames(test_summary),colnames(phos))
phos = phos[,tmp]
test_summary = test_summary[tmp,]

# tmp = apply(phos,1,function(x){ sum(is.na(x))/ncol(phos) < 0.8 })
# phos = phos[tmp,]

########
phos_SNF4 = phos[,rownames(test_summary[test_summary$inferSNF == "SNF4",])]
phos_Others = phos[,rownames(test_summary[test_summary$inferSNF != "SNF4",])]

wil.res = as.data.frame(matrix(nrow = nrow(phos_SNF4), ncol = 5))
rownames(wil.res) = rownames(phos_SNF4)
colnames(wil.res) = c("SNF4","Others","Diff","pvalue","padj")

for ( k in rownames(wil.res)){
  mat = data.frame(gene = c(as.numeric(phos_SNF4[k,]), as.numeric(phos_Others[k,]) ) ,
                   group = c(rep("SNF4", ncol(phos_SNF4)), rep("Others",ncol(phos_Others)))
  )
  mat = mat[!is.na(mat$gene),]
  if(length(unique(mat$group)) != 2){ next }
  
  res = wilcox.test(gene~group, data = mat,exact = FALSE )
  wil.res[k,"pvalue"] = res$p.value
  wil.res[k,"SNF4"] = mean(mat$gene[mat$group == "SNF4"])
  wil.res[k,"Others"] = mean(mat$gene[mat$group == "Others"])
}
wil.res$Diff = wil.res[,"SNF4"] - wil.res[,"Others"]
wil.res$padj = p.adjust(wil.res$pvalue,method = "fdr")

wil.res = cbind(wil.res, phos_annot[rownames(wil.res),])
write.csv(wil.res, file = paste0("./DEPhos.Wilcox.CPTAC.SNF4VsOthers.csv"))

# tmp1 = c("TGFBR3L", "MERTK", "ACVR1C", "CSF1R", "HJV", "EFNA3", "EFNA4", "EFNB3", "EGFR", "EPHA2", "ENG", "EPHA1", "EPHA3", "EPHA4", "EPHA5", "EPHA7", "EPHA8",
#         "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6", "ERBB2", "ERBB3", "ERBB4", "EFEMP1", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "FLT1", "FLT3", "FLT4", "ALK", "SOSTDC1", "AMHR2", "EPHA10", "EPHA6", "IGF1R", "IGF2R", "INSR", "INSRR", "KDR", "KIT", "LTBP1", "LTK", "MET", "MST1R", "MUSK", "NTRK1", "NTRK2", "NTRK3", "ROR2", "DDR2", "CRIM1", "PDGFRA", "PDGFRL", "PDGFRB", "FGFRL1", "AXL", "RET", "ROS1", "BMPR1A", "BMPR1B", "BMPR2", "TEK", "TGFBR1", "TGFBR2", "TGFBR3", "TIE1", "TYRO3", "DDR1", "LTBP4", "NRP2", "NRP1", "ACVR1", "ACVR1B", "ACVR2A", "ACVR2B", "ACVRL1")
# tmp2 = c("MERTK", "CSF1R", "EFNA3", "EFNA4", "EFNB3", "EGFR", "EPHA2", "EPHA1", "EPHA3", "EPHA4", "EPHA5", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "EPHB3", "EPHB4", "EPHB6", "ERBB2", "ERBB3", "ERBB4", "EFEMP1", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "FLT1", "FLT3", "FLT4", "ALK", "EPHA10", "EPHA6", "IGF1R", "IGF2R", "INSR", "INSRR", "KDR", "KIT", "LTK", "MET", "MST1R", "MUSK", "NTRK1", "NTRK2", "NTRK3", "ROR2", "DDR2", "CRIM1", "PDGFRA", "PDGFRL", "PDGFRB", "FGFRL1", "AXL", "RET", "ROS1", "TEK", "TIE1", "TYRO3", "DDR1", "NRP2", "NRP1")
# tmp3 = c("HIPK3", "TNK2", "TESK2", "MERTK", "CAMKK2", "CLK1", "CLK2", "CLK3", "CSF1R", "CSK", "HIPK4", "DYRK1A", "EFNA3", "EFNA4", "EFNB3", "EGFR", "EPHA2", "EPHA1", "EPHA3", "EPHA4", "EPHA5", "EPHA7", "EPHA8", "EPHB1", "EPHB2", "HIPK1", "EPHB3", "EPHB4", "EPHB6", "ERBB2", "ERBB3", "ERBB4", "PTK2B", "EFEMP1", "FER", "FES", "FGFR1", "FGFR3", "FGFR2", "FGFR4", "FGR", "FLT1", "FLT3", "FLT4", "ALK", "FRK", "ABL1", "FYN", "DSTYK", "ABL2", "EPHA10", "EPHA6", "HIPK2", "HCK", "IGF1R", "IGF2R", "INSR", "INSRR", "ITK", "JAK1", "JAK2", "JAK3", "KDR", "KIT", "LCK", "LTK",
#          "LYN", "MATK", "MET", "MST1R", "MUSK", "NEK1", "NTRK1", "NTRK2", "NTRK3", "ROR2", "DDR2", "WEE2", "CRIM1", "PDGFRA", "PDGFRL", "PDGFRB", "FGFRL1", "STYK1", "AXL", "PRKCD", "MAP2K1", "MAP2K2", "MAP2K3", "MAP2K5", "MAP2K6", "MAP2K7", "EIF2AK2", "CLK4", "SCYL1", "PTK2", "PTK6", "RET", "ROS1", "BLK", "MAP2K4", "SLA", "BMX", "SRC", "SRMS", "SYK", "BTK", "TEC", "TEK", "TESK1", "TIE1", "TTK", "TTN", "TXK", "TYK2", "TYRO3", "WEE1", "YES1", "ZAP70", "DDR1", "PEAK1", "DYRK3", "DYRK2", "TTBK1", "STK16", "TNK1", "RIPK2", "DYRK4", "NRP2", "NRP1", "BAZ1B", "PKDCC", "DYRK1B", "AATK", "MELK")
# tmp = union(union(tmp1,tmp2),tmp3)
gmt = read.gmt("/Users/ZYF/c5.go.mf.v7.4.symbols.gmt")
RTK = gmt[gmt$term == "GOMF_TRANSMEMBRANE_RECEPTOR_PROTEIN_KINASE_ACTIVITY",]
rownames(RTK) = RTK$gene
kinase = openxlsx::read.xlsx("/Users/ZYF/Desktop/研究生科研/研一/lumianl平台建立/生信分析/30.NG修改/分型外推/CPTAC_2020/激酶kinase 参考文献：Driver Fusions and Their Implications in the Development and Treatment of Human Cancers.xlsx",1)
colnames(kinase) = kinase[1,] ; kinase = kinase[-1,] 
RTK = intersect(rownames(RTK),kinase$Kinase)

wil.res2 = wil.res[wil.res$GeneSymbol %in% RTK,]
write.csv(wil.res2,file = "DEPhos.RTK.kinase.Wilcox.CPTAC.SNF4VsOthers.csv")

