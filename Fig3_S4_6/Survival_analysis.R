rm(list = ls()) ; graphics.off()
library(ggplot2)
library(survminer)
library(survival)
library("export")
library(forestplot)
load("CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
luminal$DMFS_status = as.numeric(luminal$DMFS_status)
luminal$DMFS_months = as.numeric(luminal$DMFS_months)
luminal$RFS_status = as.numeric(luminal$RFS_status)
luminal$RFS_months = as.numeric(luminal$RFS_months)
luminal$SNF = paste0("SNF", SNF_Cluster[rownames(luminal)])
luminal$PAM50 = factor(luminal$PAM50_classifier,levels = c("LumA","LumB","Her2","Basal","Normal"))
color_PAM50 = c(Normal = "#CDCFD0" , LumB = "#76CFE6" , LumA = "#1D76BC" , Her2  = "#6E59A6", Basal = "#E11D2E")

############################ km plot -----
for ( i in unique(luminal$PAM50_classifier)){
  mat = luminal[luminal$PAM50_classifier == i,]
  km = survfit(Surv(DMFS_months,DMFS_status) ~ SNF,data = mat)
  p = ggsurvplot(km,pval = T,pval.method = T,palette =  c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A"),
             break.x.by = 12)
  ggsave(p$plot,filename = paste0("surplot_SNF_DMFS_in_",i,".pdf"),height = 6,width = 6)
  graph2ppt(p$plot,file =paste0("surplot_SNF_DMFS_in_",i,".pptx"),height = 6,width = 6)

}

km = survfit(Surv(DMFS_months,DMFS_status) ~ SNF,data = luminal)
p = ggsurvplot(km,pval = T,pval.method = T,palette =  c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A"),
               break.x.by = 12)
ggsave(p$plot,filename = paste0("surplot_SNF_DMFS.pdf"),height = 6,width = 6)
graph2ppt(p$plot,file =paste0("surplot_SNF_DMFS.pptx"),height = 6,width = 6)

km = survfit(Surv(RFS_months,RFS_status) ~ SNF,data = luminal)
p = ggsurvplot(km,pval = T,pval.method = T,palette =  c("#2378B3",  "#1CB038",  "#F8A900", "#D5271A"),
               break.x.by = 12)
ggsave(p$plot,filename = paste0("surplot_SNF_RFS.pdf"),height = 6,width = 6)
graph2ppt(p$plot,file =paste0("surplot_SNF_RFS.pptx"),height = 6,width = 6)



############################ COX ---------
## COX forest --------
colnames(luminal)[colnames(luminal) == "辅助内分泌治疗"  ] ="EndoTherapy" 
luminal = luminal[luminal$EndoTherapy %in% c("Yes","No"),]  

colnames(luminal)[colnames(luminal) == "淋巴结转移数" ] ="LN_count" 
luminal$LN_count = as.numeric(luminal$LN_count)

colnames(luminal)[colnames(luminal) == "辅助化疗" ] ="ChemoTherapy" 
luminal = luminal[luminal$ChemoTherapy %in% c("Yes","No"),]

colnames(luminal)[colnames(luminal) == "肿瘤大小.cm." ] ="TumorSize" 
luminal$TumorSize = as.numeric(luminal$TumorSize)

colnames(luminal)[colnames(luminal) == "组织学分级" ] ="Grade" 
luminal$Grade[luminal$Grade  %in% c("Unknown","Unknown（组织学类型不明）")] = "Unknown"
luminal$Grade[luminal$Grade == "1.5"] = "2"
luminal$Grade[luminal$Grade == "2.5"] = "3"
luminal$Grade = as.numeric(luminal$Grade)

colnames(luminal)[colnames(luminal) == "PAM50_classifier" ] ="PAM50" 
luminal$PAM50[! luminal$PAM50 %in% c("LumA","LumB")] = "nonLum"
#luminal$PAM50 = factor(luminal$PAM50,levels = c("LumA","LumB","Her2","Basal","Normal"))
luminal$PAM50 = factor(luminal$PAM50,levels = c("LumA","LumB","nonLum"))

luminal = luminal[luminal$EndoTherapy == "Yes",]

res.cox = coxph(Surv(DMFS_months,DMFS_status)~ SNF + LN_count + ChemoTherapy+TumorSize+Grade+PAM50 , data = luminal)
res.cox = coxph(Surv(DMFS_months,DMFS_status)~ SNF + LN_count + ChemoTherapy+TumorSize+Grade , data = luminal)
res.cox = summary(res.cox)

tmp1 = res.cox[["coefficients"]]
tmp2 = res.cox[["conf.int"]]
cox.res = as.data.frame(cbind(tmp1[,c(2,5)],tmp2[,3:4]))
colnames(cox.res) = c("HR","p.val","lower.95","upper.95")

write.csv(cox.res,file ='./SNF_TS_LN_CT_COX4plot_pickET.csv' )
######### forest plot
cox.res$p.val = round(cox.res$p.val, 2)
cox.res$factor = rownames(cox.res)
#pdf("SNF_TS_LN_CT_multiCOX_forest_mergeNoLum_pickET.pdf",width = 8,height = 4)
pdf("SNF_TS_LN_CT_multiCOX_forest_pickET.pdf",width = 8,height = 4)
forestplot(labeltext= as.matrix(cox.res[,c(5,2)]), graph.pos=2,
           mean=c(cox.res$HR),
           lower=c(cox.res$lower.95), upper=c(cox.res$upper.95),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"), 
           zero=1, 
           cex=0.5, lineheight = "auto",
           colgap=unit(8,"mm"),
           lwd.ci=2, boxsize=0.3, 
           ci.vertices=TRUE,
           ci.vertices.height = 0.1) 

graph2ppt(file = "SNF_TS_LN_CT_multiCOX_forest_pickET.ppt",width = 8,height = 4)
dev.off()
