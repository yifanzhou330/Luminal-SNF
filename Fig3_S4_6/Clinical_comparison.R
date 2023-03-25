rm(list = ls()) ; graphics.off()
library(ggplot2)
library(export)
library(plyr)

load("/CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20220829.Rdata")
luminal$SNF = paste0("SNF",SNF_Cluster[rownames(luminal)])

# age-------
luminal$Age_2[luminal$Age <= 40] = "Age:<=40"
luminal$Age_2[luminal$Age > 40] = "Age:>40"
cluster.merge4plot = as.data.frame(table(luminal$SNF, luminal$Age_2))
colnames(cluster.merge4plot)[1:2] = c("SNF","Age")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Age = factor(cluster.merge4plot2$Age,levels = c("Age:<=40","Age:>40"))
#color = c("#DBBEBE", "#C88E8D", "#B75C5B", "#A52B29")
color = c("#DBBEBE",  "#A52B29")
names(color) = c("Age:<=40","Age:>40")
set.seed(123)
tmp = fisher.test(table(luminal$SNF, luminal$Age_2),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Age )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename ="./SNF_Age_4plot_percent_barplot_2.pdf")


# Ki_67-----
luminal$Ki_67_revise = luminal$Ki67阳性百分比 
luminal_tmp = luminal[!is.na(luminal$Ki_67_revise),]
luminal_tmp$Ki_67_revise[luminal_tmp$Ki67阳性百分比 < 15] = "ki_67:low"
luminal_tmp$Ki_67_revise[luminal_tmp$Ki67阳性百分比 < 30 & luminal_tmp$Ki67阳性百分比 >= 15 ] = "ki_67:middle"
luminal_tmp$Ki_67_revise[luminal_tmp$Ki67阳性百分比 >= 30] = "ki_67:high"
cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$Ki_67_revise))
colnames(cluster.merge4plot)[1:2] = c("SNF","Ki_67")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Ki_67 = factor(cluster.merge4plot2$Ki_67,levels = c("ki_67:low","ki_67:middle","ki_67:high"))
color = c("#EDD4BC", "#EFA463","#F27900")
names(color) = c("ki_67:low","ki_67:middle","ki_67:high")
set.seed(123)
tmp = fisher.test(table(luminal$SNF, luminal$Ki_67_revise),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Ki_67 )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename ="./SNF_Ki_67_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_Ki_67_4plot_percent_barplot.ppt")

# Menopause_status-----
luminal_tmp = luminal[luminal$绝经状态 %in% c("Yes","No"),]
luminal_tmp$绝经状态[luminal_tmp$绝经状态 %in% c("Yes")] = "TRUE"
luminal_tmp$绝经状态[luminal_tmp$绝经状态 %in% c("No")] = "FALSE"

cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$绝经状态))
colnames(cluster.merge4plot)[1:2] = c("SNF","Menopause_status")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Menopause_status = factor(cluster.merge4plot2$Menopause_status,levels = c("TRUE","FALSE"))
color = c("#7F7E7E", "#C4C4C4")
names(color) = c("TRUE","FALSE")
set.seed(123)
tmp = fisher.test(table(luminal_tmp$SNF, luminal_tmp$绝经状态),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Menopause_status )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename ="./SNF_Menopause_status_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_Menopause_status_4plot_percent_barplot.ppt")

# LVI-----
luminal_tmp = luminal[luminal$脉管癌栓 %in% c("Yes","No"),]
luminal_tmp$脉管癌栓[luminal_tmp$脉管癌栓 %in% "Yes"] = "TRUE"
luminal_tmp$脉管癌栓[luminal_tmp$脉管癌栓 %in% "No"] = "FALSE"
cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$脉管癌栓))
colnames(cluster.merge4plot)[1:2] = c("SNF","LVI")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$LVI = factor(cluster.merge4plot2$LVI,levels = c("TRUE","FALSE"))
color = c("#7F7E7E", "#C4C4C4")
names(color) = c("TRUE","FALSE")
set.seed(123)
tmp = fisher.test(table(luminal_tmp$SNF, luminal_tmp$脉管癌栓),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = LVI )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename ="./SNF_LVI_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_LVI_4plot_percent_barplot.ppt")


# Adjuvant_chemotherapy-----
luminal_tmp = luminal[luminal$辅助化疗 %in% c("Yes","No"),]
luminal_tmp$辅助化疗[luminal_tmp$辅助化疗 %in% "Yes"] = "TRUE"
luminal_tmp$辅助化疗[luminal_tmp$辅助化疗 %in% "No"] = "FALSE"
cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$辅助化疗))
colnames(cluster.merge4plot)[1:2] = c("SNF","Adjuvant_chemotherapy")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Adjuvant_chemotherapy = factor(cluster.merge4plot2$Adjuvant_chemotherapy,levels = c("TRUE","FALSE"))
color = c("#7F7E7E", "#C4C4C4")
names(color) = c("TRUE","FALSE")
set.seed(123)
tmp = fisher.test(table(luminal_tmp$SNF, luminal_tmp$辅助化疗),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Adjuvant_chemotherapy )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75),
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename ="./SNF_Adjuvant_chemotherapy_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_Adjuvant_chemotherapy_4plot_percent_barplot.ppt")


# Radiotherapy-----
luminal_tmp = luminal[luminal$辅助放疗 %in% c("Yes","No"),]
luminal_tmp$辅助放疗[luminal_tmp$辅助放疗 %in% "Yes"] = "TRUE"
luminal_tmp$辅助放疗[luminal_tmp$辅助放疗 %in% "No"] = "FALSE"
cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$辅助放疗))
colnames(cluster.merge4plot)[1:2] = c("SNF","Radiotherapy")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Radiotherapy = factor(cluster.merge4plot2$Radiotherapy,levels = c("TRUE","FALSE"))
color = c("#7F7E7E", "#C4C4C4")
names(color) = c("TRUE","FALSE")
set.seed(123)
tmp = fisher.test(table(luminal_tmp$SNF, luminal_tmp$辅助放疗),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Radiotherapy )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
  ggsave(p,filename ="./SNF_Radiotherapy_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_Radiotherapy_4plot_percent_barplot.ppt")



# Grade -----
luminal_tmp = luminal[luminal$组织学分级 != "Unknown",]
luminal_tmp$Grade = luminal_tmp$组织学分级
luminal_tmp$Grade[luminal_tmp$Grade %in% "2.5"] = 3
luminal_tmp$Grade[luminal_tmp$Grade %in% "1.5"] = 2
luminal_tmp$Grade = as.numeric(luminal_tmp$Grade)
luminal_tmp$Grade[luminal_tmp$Grade > 2] = ">G2"
luminal_tmp$Grade[luminal_tmp$Grade != ">G2" ] = "<=G2"

cluster.merge4plot = as.data.frame(table(luminal_tmp$SNF, luminal_tmp$Grade))
colnames(cluster.merge4plot)[1:2] = c("SNF","Grade")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$Grade = factor(cluster.merge4plot2$Grade,levels = c("<=G2",">G2"))
color = c("#BDD7EE", "#2E75B6")
names(color) = c("<=G2",">G2")
set.seed(123)
tmp = fisher.test(table(luminal_tmp$SNF, luminal_tmp$Grade),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = Grade )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
  ggsave(p,filename ="./SNF_Grade_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_Grade_4plot_percent_barplot.ppt")

  
# pT -----
luminal$pT = luminal$pT.仅计算浸润成分.不计原位癌.
luminal$pT[luminal$pT %in% c("pT1a", "pT1b", "pT1c", "pT1mi")] = "pT1"
cluster.merge4plot = as.data.frame(table(luminal$SNF, luminal$pT))
colnames(cluster.merge4plot)[1:2] = c("SNF","pT")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$pT = factor(cluster.merge4plot2$pT,levels = c("pT1" , "pT2" ,"pT3" ))
color = c("#DBBEBE", "#C27978", "#A52B29")
names(color) = c("pT1" , "pT2" ,"pT3" )
set.seed(123)
tmp = fisher.test(table(luminal$SNF, luminal$pT),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = pT )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
  ggsave(p,filename ="./SNF_pT_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_pT_4plot_percent_barplot.ppt")


# pN -----
cluster.merge4plot = as.data.frame(table(luminal$SNF, luminal$pN))
colnames(cluster.merge4plot)[1:2] = c("SNF","pN")

cluster.merge4plot2 = ddply(cluster.merge4plot,"SNF", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$pN = factor(cluster.merge4plot2$pN,levels = c("pN0" , "pN1", "pN2", "pN3" ))
color = c("#E7DCE6", "#E0C7DC", "#D8B4D0",  "#D1A1C7")
names(color) = c("pN0" , "pN1", "pN2", "pN3" )
set.seed(123)
tmp = fisher.test(table(luminal$SNF, luminal$pN),simulate.p.value=T)
p = ggplot(cluster.merge4plot2, aes(x=SNF, y = percent, fill = pN )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.75), 
        axis.ticks = element_line(size = 0.75),
        axis.text.x=element_text(size = 7),
        axis.text.y=element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],5)) ) 
ggsave(p,filename = "./SNF_pN_4plot_percent_barplot.pdf")
#graph2ppt(p, file = "SNF_pN_4plot_percent_barplot.ppt")













