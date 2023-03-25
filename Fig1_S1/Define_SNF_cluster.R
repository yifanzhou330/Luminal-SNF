rm(list = ls() ) ; graphics.off()
library(paletteer) 
library(export)
library(ggpubr)
load("./CBCGA_HRposHER2neg351_WES_RNAseq_CNV_Metab_Protein_20221010.Rdata")
################################
# data preapre 
################################
# cnv
cnv = cnv.alldata

# exp
exp.fpkm.TT.log = log2(exp.fpkm.TT +1)
exp.fpkm.TT.log.SD2k = exp.fpkm.TT.log[order(apply(exp.fpkm.TT.log,1,sd),decreasing = T)[1:2000],] 

# polar
polar = polar_metabolite_TT_MS2_log2

# data transform
cnv = as.data.frame(t(cnv))
exp.fpkm.TT.log.SD2k = as.data.frame(t(exp.fpkm.TT.log.SD2k))
polar = as.data.frame(t(polar))

all(rownames(cnv) == rownames(exp.fpkm.TT.log.SD2k)) 
all(rownames(cnv) == rownames(polar)) 

rm(list = ls()[! ls() %in% c("luminal","cnv","exp.fpkm.TT.log.SD2k","polar")])
param = "logFPKMsd2k_CNValldata_Polar"

################################
# SNF subtype 
################################
library(SNFtool)
library(tidyr)
library(CancerSubtypes)

# First, set all the parameters:
K = 15;		
alpha = 0.5;  	
T = 20; 	

cnv.m.n = standardNormalization(cnv)
exp.fpkm.m.n = standardNormalization(exp.fpkm.TT.log.SD2k)
polar.m.n = standardNormalization(polar)

# If the data is continuous, we recommend to use the function "dist2" as follows 
dist_cnv = (SNFtool::dist2(cnv.m.n,cnv.m.n))
dist_exp = (SNFtool::dist2(as.matrix(exp.fpkm.m.n),as.matrix(exp.fpkm.m.n)))
dist_polar = (SNFtool::dist2(as.matrix(polar.m.n), as.matrix((polar.m.n))))

# next, construct similarity graphs
W_cnv = affinityMatrix(dist_cnv, K, alpha)
W_exp = affinityMatrix(dist_exp, K, alpha)
W_polar = affinityMatrix(dist_polar, K, alpha)

W = SNF(list(W_cnv, W_exp, W_polar), K, T)
#write.csv(W,file = paste0(param,"_W.csv" ))


################################
# choose the best num of clusters
################################
library(Spectrum)
library(plot3D)
pdf(paste0("./best_num_Spectral_cluster.pdf"))
r <- Spectrum(W,method =1,diffusion = TRUE ,kerneltype = "stsc",NN = 20)
dev.off()
# optimal K: 4

C = 4			# number of clusters
set.seed(1234)
group = spectralClustering(W, C) 	# the final subtypes information

names(group) <- row.names(exp.fpkm.m.n)
SNF_Cluster  <- group
SNF_Cluster[group == 1] = 3
SNF_Cluster[group == 2] = 1
SNF_Cluster[group == 4] = 2
SNF_Cluster[group == 3] = 4


################################
# estimate Silhouette width and visulize
################################
# displayClustersWithHeatmap
set.seed(321)
pdf(paste0(param,"_Re_nonRe_displayClustersWithHeatmap.C",C,".pdf"))
displayClustersWithHeatmap(W = log(W), group = SNF_Cluster, 
                           col = colorRampPalette(c("#4F4F4F","#595959","#7F7F7F","#666666","#999999","#CCCCCC","#D9D9D9","#ECECEC"))(20),
                           labRow = F,labCol = F) 
dev.off()

# Silhouette
palette4silhou = c("#3D76AE","#53AD4A" ,"#EDAB3C" ,"#C3392A",rainbow(n = 5))
pdf(paste0(param,"_Re_nonRe_silhouette_SimilarityMatrix.C",C,".pdf"))
sil <- silhouette_SimilarityMatrix(SNF_Cluster, W)
plot(sil, col = palette4silhou[1:C])
dev.off()








