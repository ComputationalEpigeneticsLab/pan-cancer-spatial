library(rlang)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(scales)
library(RColorBrewer)

dir_NMFre<-'/panCancerST/NMF/'
cir_intra_min=23
cir_inter_min=23

MP_list<-read.delim(paste0(dir_NMFre,'final/','intra',cir_intra_min,'_inter',cir_inter_min,'MP_list.txt'),
                    stringsAsFactors = F,check.names = F)
nmf_programs<-read.delim(paste0(dir_NMFre,'final/','intra',cir_intra_min,'_inter',cir_inter_min,'nmf_programm.txt'),
                         stringsAsFactors = F,check.names = F)
nmf_intersect_original<-read.delim(paste0(dir_NMFre,'final/','intra',cir_intra_min,'_inter',cir_inter_min,'nmf_inter.txt'),
                                   stringsAsFactors = F,check.names = F)
inds_new<-read.delim(paste0(dir_NMFre,'final/','intra',cir_intra_min,'_inter',cir_inter_min,'inds_new.txt'),
                     stringsAsFactors = F,check.names = F)
inds_new<-inds_new$x
Cluster_list<-readRDS(paste0(dir_NMFre,'final/','intra',cir_intra_min,'_inter',cir_inter_min,'Cluster_list.rds'))

plot_data<-nmf_intersect_original[inds_new,inds_new]
all_MP_list<-c()
for(j in 1:length(Cluster_list)){
  #i=1
  aaa<-data.frame(dataset_slice=Cluster_list[[j]],MP_list=paste0("MP_",j))
  all_MP_list<-rbind(all_MP_list,aaa)
}

plot_data<-plot_data[all_MP_list$dataset_slice,all_MP_list$dataset_slice]
plot_data = as.matrix(plot_data)
nmf_intersect_meltI_NEW <- reshape2::melt(plot_data[,ncol(plot_data):1])

custom_magma <- c(colorRampPalette(c("white", "white"))(180), rev(magma(180, begin = 0.18)))
p <- ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(0,50), low=custom_magma[1: 150],  mid =custom_magma[151:200], high = custom_magma[201: 360], midpoint = 13.5,oob=squish,  name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(0,50), low=custom_magma[1: 150],  mid =custom_magma[151:200], high = custom_magma[201: 360], midpoint = 13.5, oob=squish,  name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))

pdf('/panCancerST/NMF/final/MP.pdf',width =5,height = 4)
print(p)
dev.off()