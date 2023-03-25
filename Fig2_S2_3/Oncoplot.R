rm(list = ls())
graphics.off()
library(maftools)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

load("./CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
luminal$SNF = paste0("SNF",SNF_Cluster[rownames(luminal)])

######## annot prepare -------
luminal$HRD = as.numeric(luminal$HRD)

######## prepre mut plot  --------
luminal_mut = luminal
luminal_mut$Tumor_Sample_Barcode = luminal_mut$PatientCode
luminal_mut = luminal_mut[intersect(rownames(luminal_mut),unique(maf.df$Tumor_Sample_Barcode)),]
SNF_Cluster_mut = SNF_Cluster[rownames(luminal_mut)]

maf.df = maf.df[maf.df$Tumor_Sample_Barcode %in% rownames(luminal_mut),]
SNF_maf = read.maf(maf.df,clinicalData = luminal_mut)

oncoplot(SNF_maf,top = 10,
         clinicalFeatures = "SNF",sortByAnnotation = T,writeMatrix = T,removeNonMutated = F,
         groupAnnotationBySize = F)
mut_mtx0 = read.table("./onco_matrix.txt",sep = "\t", na.strings = c("",0))
tmp = matrix("Unknown",nrow = nrow(mut_mtx0), ncol = nrow(luminal) -ncol(mut_mtx0) )
rownames(tmp) = rownames(mut_mtx0) ; colnames(tmp) = setdiff(rownames(luminal), colnames(mut_mtx0))
mut_mtx0 = cbind(mut_mtx0,tmp)
genelist1 = read.csv("/Census_allMon Jan 11 04_06_52 2021.csv") 
genelist1 = unique(genelist1$Gene.Symbol)
genelist2 = data.table::fread("NCG6_cancergenes.tsv",data.table = F) 
genelist2 = unique(genelist2$symbol)
genelist3 = data.table::fread("./FLAGS supplementary.txt",data.table = F)
genelist3 = genelist3$FLAGS
genelist = setdiff(intersect(genelist1, genelist2),genelist3)
genelist = intersect(rownames(mut_mtx0),genelist)
mut_mtx0 = mut_mtx0[genelist,rownames(luminal)]

mut_mtx =  mutCountMatrix(SNF_maf, includeSyn = F, countOnly = NULL, removeNonMutated = F)
mut_mtx = mut_mtx_stat = apply(mut_mtx, 2, function(x){ifelse(x == 0, "WT", "Mut")} )
tmp = matrix("Unknown",nrow = nrow(mut_mtx), ncol = nrow(luminal) -ncol(mut_mtx) )
rownames(tmp) = rownames(mut_mtx) ; colnames(tmp) = setdiff(rownames(luminal), colnames(mut_mtx))
mut_mtx = cbind(mut_mtx,tmp)
mut_mtx = mut_mtx[ genelist ,rownames(luminal)]
mut_mtx_stat = mut_mtx_stat[genelist,]

######## define order -------
tmp = luminal ; tmp$Tumor_Sample_Barcode = tmp$PatientCode
tmp = tmp[unique(maf.df$Tumor_Sample_Barcode),]
tmp = read.maf(maf.df,clinicalData = tmp)
oncoplot(tmp,genes = rownames(mut_mtx),clinicalFeatures = "SNF",sortByAnnotation = T,writeMatrix = T,removeNonMutated = F,
         groupAnnotationBySize = F)
tmp = read.table("./onco_matrix.txt",sep = "\t")
order = colnames(tmp)
no_mut = luminal[setdiff(rownames(luminal),order),"SNF"]
names(no_mut) = setdiff(rownames(luminal),order)

tmp = luminal[order,"SNF"]
names(tmp) = order

order = c(tmp,no_mut)
order = sort(order)


######## prepare cnv plot  --------
amp_mtx = cnv.alldata[c("CCND1","MDM2","FGFR1"),]
amp_mtx = apply(amp_mtx, 2, function(x){ifelse(x < log2(4/2), "non_Amp", "Amp")} )
amp_mtx = amp_mtx[,names(order)]


######### merge & draw ------
color_PAM50 = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
names(color_PAM50) = c("LumA","LumB","Her2","Basal","Normal")
color_SNF = color
names(color_SNF) = paste0("SNF",1:4)

col_fun_Isig <- colorRamp2(c(-0.5, 0.3 ,2), colorRampPalette(c("white", "#08506F"))(20)[c(2,10,20)])
col_fun_SSig <- colorRamp2(c(-0.1, 0.75 ,1.2), colorRampPalette(c("white", "#50A9AC"))(20)[c(2,10,20)])
col_fun_hrd <- colorRamp2(c(0, 16, 75), c('#FDF9E7', '#EFE089', '#D7C801'))

#### mut ------
col = c("Missense_Mutation" = "#366A9C", "Nonsense_Mutation" = "#BBDE93",
        "In_Frame_Del" = "#EE8632" ,"In_Frame_Ins" = "#D0342B",
        "Multi_Hit" = "black", "Splice_Site" = "#F3C17B", "Frame_Shift_Del" = "#AECDE1",
        "Frame_Shift_Ins" = "#ED9E9B",
        "Nonstop_Mutation" = "#339900",
        "Unknown" = "#eaeaea")
wid = 3
hei = 2

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = "white", col = "white"))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Missense_Mutation"], col = col["Missense_Mutation"]))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Nonsense_Mutation"], col = col["Nonsense_Mutation"]))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["In_Frame_Del"], col = col["In_Frame_Del"]))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["In_Frame_Ins"], col = col["In_Frame_Ins"]))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Multi_Hit"], col = col["Multi_Hit"]))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Splice_Site"], col = col["Splice_Site"]))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Frame_Shift_Del"], col = col["Frame_Shift_Del"]))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Frame_Shift_Ins"], col = col["Frame_Shift_Ins"]))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Nonstop_Mutation"], col = col["Nonstop_Mutation"]))
  },
  Unknown = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col["Unknown"], col = col["Unknown"]))
  }
)

luminal = luminal[names(order),]
p1 <- oncoPrint(mut_mtx0[,names(order)],  
                alter_fun = alter_fun, col = col, show_pct = T, pct_side = 'right',
                column_title = '', column_order = names(order),
                row_names_side = "left", show_column_names = F, remove_empty_columns = F,
                heatmap_legend_param = list(title = "Mut"),border = T,
                row_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 0),
                column_split = factor(luminal[names(order),"SNF"], levels = c('SNF1', 'SNF2', 'SNF3', 'SNF4')),column_gap = unit(2, "mm"),
                top_annotation = HeatmapAnnotation(SNF = luminal$SNF,
                                                   HRD = luminal$HRD,
                                                   annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9),
                                                   col = list( SNF = color_SNF, HRD = col_fun_hrd) ) )

#### cnv ------
amp_mtx2 = amp_mtx
for (i in colnames(amp_mtx2)){
  mat = amp_mtx2[,i]
  mat[mat == "non_Amp"] = ""
  amp_mtx2[,i] = mat
}

col2 = c( "Amp" = "#BC102B")
wid = 3
hei = 2
alter_fun2 <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = "white", col = "white"))
  },
  Amp = function(x, y, w, h) {
    grid.rect(x, y, w-unit(wid, "pt"), h-unit(hei, "pt"),
              gp = gpar(fill = col2["Amp"], col = col2["Amp"]))
  }
  )
p2 = oncoPrint(amp_mtx2[,names(order)],  
               alter_fun = alter_fun2,col = col2,
             show_pct = T, 
             show_column_names = FALSE,
             na_col = "#eaeaea", border = T, 
             row_names_gp = gpar(fontsize = 9),
             heatmap_legend_param = list(title = "CNV"),
             row_names_side = "left",
             pct_side = 'right',
             top_annotation = NULL,
             column_split = factor(luminal[names(order),"SNF"], levels = c('SNF1', 'SNF2', 'SNF3', 'SNF4')),
             column_gap = unit(2, "mm")
)

hp_list = p1 %v% p2
export::graph2pdf(hp_list,file =  "Mut_CNA_plot.pdf",width = 14, height = 7)


######### stat ------------
filename = "Mut_CNA_plot"

## mut
stat_res = as.data.frame(matrix(nrow = nrow(mut_mtx_stat),ncol = 6))
rownames(stat_res) = rownames(mut_mtx_stat)
colnames(stat_res) = c("SNF1","SNF2","SNF3","SNF4","p.val","adj.p")
for( i in rownames(stat_res)){
  mat = data.frame(gene = mut_mtx_stat[i,names(sort(SNF_Cluster_mut))], 
                   group = paste0("SNF",sort(SNF_Cluster_mut)))
  res = chisq.test(table(mat$gene,mat$group))
  tmp = as.matrix.data.frame(table(mat$gene,mat$group))
  
  stat_res[i,"SNF1"] = tmp[1,1]/sum(tmp[,1])
  stat_res[i,"SNF2"] = tmp[1,2]/sum(tmp[,2])
  stat_res[i,"SNF3"] = tmp[1,3]/sum(tmp[,3])
  stat_res[i,"SNF4"] = tmp[1,4]/sum(tmp[,4])
  stat_res[i,"p.val"] = res[["p.value"]]
}
stat_res[,"adj.p"] = p.adjust(stat_res[,"p.val"],method = "fdr")
write.csv(stat_res, file =  paste0(filename,"_mut_chisq.csv"))

## CNA
stat_res = as.data.frame(matrix(nrow = nrow(amp_mtx),ncol = 6))
rownames(stat_res) = rownames(amp_mtx)
colnames(stat_res) = c("SNF1","SNF2","SNF3","SNF4","p.val","adj.p")
for( i in rownames(stat_res)){
  mat = data.frame(gene = amp_mtx[i,names(sort(SNF_Cluster))], 
                   group = paste0("SNF",sort(SNF_Cluster)))
  res = chisq.test(table(mat$gene,mat$group))
  
  tmp = as.matrix.data.frame(table(mat$gene,mat$group))
  stat_res[i,"SNF1"] = tmp[1,1]/sum(tmp[,1])
  stat_res[i,"SNF2"] = tmp[1,2]/sum(tmp[,2])
  stat_res[i,"SNF3"] = tmp[1,3]/sum(tmp[,3])
  stat_res[i,"SNF4"] = tmp[1,4]/sum(tmp[,4])
  stat_res[i,"p.val"] = res[["p.value"]]
}
stat_res[,"adj.p"] = p.adjust(stat_res[,"p.val"],method = "fdr")
write.csv(stat_res, file =  paste0(filename,"_AMP_chisq.csv"))







