####演讲的各种图
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)





dir_pic<-'E:/Mirror/lecture/pic/'
dir_data<-'E:/Mirror/lecture/data/'

####3
library(clusterProfiler) # 富集分析R包
library(stringr) # 标签换行
library(AnnotationDbi)
library(org.Hs.eg.db) # 参考物种的基因注释数据库hg19
library(DOSE)
library(ggplot2) # 绘图
library(ggrepel) # 
enrich_data<-read.delim(paste0(dir_data,"hypergeom_re_sig_ues.txt"),stringsAsFactors = F,check.names = F)
GO_BP<-read.delim('E:/Mirror/ST_analysis/data/geneset/pathway_go_bp.txt',stringsAsFactors = F,check.names = F)
GO_BP$length<-lapply(GO_BP$gene_set,function(x){
  aa<-length(unlist(strsplit(x,',')))
}) %>% unlist()
intersect(GO_BP$pathway_name,enrich_data$pathway)
enrich_data$GeneRatio<-GO_BP[match(enrich_data$pathway,GO_BP$pathway_name),'length']
enrich_data$GeneRatio<-enrich_data$inter_num/enrich_data$GeneRatio
enrich_data$logP<-(-log10(enrich_data$p_value))
range(enrich_data$logP)
enrich_data$logP[which(enrich_data$logP>20)]<-5

visual_filter<-read.csv("K:/0_之前的课题/before/Classifier/data/DEG_inter_PML_RARA/data2_M3DEG/DEG_BP/visual_filter.csv")

p <- ggplot(enrich_data,aes(GeneRatio,pathway,colour=logP))+ 
  geom_point(aes(size=inter_num,colour=logP))+
  scale_size_continuous(range=c(4,8))+
  scale_color_gradientn(colours = c(colorRampPalette(c("#FDAE61","#F46D43"))(90),
                                    colorRampPalette(c("#F46D43","#9E0142"))(90)) )+ 
  theme_bw() 
print(p)
pdf(paste0(dir_pic,'GO_BP.pdf'),width = 10,height = 6)
print(p)
dev.off()







