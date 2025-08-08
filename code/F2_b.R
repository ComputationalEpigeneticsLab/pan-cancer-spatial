###所有癌症的各区域干性绘制箱式图
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(magic)



####来源于89_All_slice_stem.R
all_slice_stem<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_stem.txt',
                           stringsAsFactors = F,check.names = F)
all_slice_stem<-all_slice_stem[all_slice_stem$localType%in%c('Boundary','Core','Dispersion'),]
table(all_slice_stem$localType)
all_box<-all_slice_stem
all_box$cancer<-'all'
all_slice_stem<-rbind(all_slice_stem,all_box)
all_slice_stem<-mutate(all_slice_stem,Legend = factor(all_slice_stem$localType, levels = c('Core','Boundary','Dispersion')))


p_box<-ggplot(all_slice_stem, aes(x = cancer, y = CytoTRACE,fill=Legend))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = Legend), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                               "Dispersion"="#EE762D","Imm"="#a9d38a"))+
  scale_color_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                                "Dispersion"="#EE762D","Imm"="#a9d38a"))+
  #ggtitle("Stemness")+
  theme(plot.title = element_text(hjust = 0.4))+
  theme(plot.title = element_text(size = 6))+
  #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
  stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.9)+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(p_box)
pdf('E:/Mirror/ST_analysis/pic/re/3/allCancer_stemness_box2.pdf',width = 14,height = 1.5)
print(p_box)
dev.off()
?geom_boxplot


?stat_compare_means


####每组中core比其他区域的秩和检验

group_cancer<-unique(all_slice_stem$cancer)
all_p_stem<-c()
for(i in 1:length(group_cancer)){#i=1
  cancer_stem<-all_slice_stem[which(all_slice_stem$cancer==group_cancer[i]),]
  p_stem<-wilcox.test(cancer_stem$CytoTRACE[which(cancer_stem$localType=='Core')],
                      cancer_stem$CytoTRACE[which(cancer_stem$localType!='Core')],alternative = 'greater')[["p.value"]]
  all_p_stem<-c(all_p_stem,p_stem)
}

all_p_stem<-data.frame(cancer=group_cancer,p_value=all_p_stem)














