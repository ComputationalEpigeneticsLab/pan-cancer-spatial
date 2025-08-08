####TAM 4个亚簇 对切片上的spot的映射
###bdy dis的step1和far_top20的spot映射结果
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)




TAM_subcluster<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/Cytospace2/skcm/TAM_subcluster_infor.txt',
                           stringsAsFactors = F,check.names = F)
TAM_subcluster$seurat_clusters<-paste0('TAM_subCluster',TAM_subcluster$seurat_clusters)
TAM_subcluster$OriginalCID<-rownames(TAM_subcluster)

Cyto_skcm12<-read.csv('E:/Mirror/ST_analysis/data/10X Visium/Cytospace/skcm12_slice1/assigned_locations.csv')
# Cyto_skcm12_scCell<-lapply(Cyto_skcm12$OriginalCID,function(x) unlist(strsplit(x,',')))
# Cyto_skcm12_scCellType<-lapply(Cyto_skcm12$CellType,function(x) unlist(strsplit(x,',')))
# table(lapply(Cyto_skcm12_scCell,length) %>% unlist()==lapply(Cyto_skcm12_scCell,length) %>% unlist())
# Cyto_skcm12<-data.frame(SpotID=rep(Cyto_skcm12$SpotID,times=lapply(Cyto_skcm12_scCell,length) %>% unlist()),
#                         OriginalCID=unlist(Cyto_skcm12_scCell),
#                         CellType=unlist(Cyto_skcm12_scCellType))
#rep(c(1,4,2),times=c(2,5,9))
Cyto_skcm12<-Cyto_skcm12[which(Cyto_skcm12$CellType=='TAM'),c('OriginalCID','SpotID')]
Cyto_skcm12<-merge(Cyto_skcm12,TAM_subcluster[,c('OriginalCID','seurat_clusters')],by='OriginalCID')

Cyto_skcm13<-read.csv('E:/Mirror/ST_analysis/data/10X Visium/Cytospace/skcm13_slice1/assigned_locations.csv')
Cyto_skcm13<-Cyto_skcm13[which(Cyto_skcm13$CellType=='TAM'),c('OriginalCID','SpotID')]
Cyto_skcm13<-merge(Cyto_skcm13,TAM_subcluster[,c('OriginalCID','seurat_clusters')],by='OriginalCID')

Cyto_skcm12$SpotID<-paste0('skcm12_',Cyto_skcm12$SpotID)
Cyto_skcm13$SpotID<-paste0('skcm13_',Cyto_skcm13$SpotID)
Cyto_skcm<-rbind(Cyto_skcm12,Cyto_skcm13)

near_skcm12<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/FarOutNear/skcm12/slice1_BdyDisMinDis.txt',stringsAsFactors = F,check.names = F)
near_skcm13<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/FarOutNear/skcm13/slice1_BdyDisMinDis.txt',stringsAsFactors = F,check.names = F)

skcm12_step1<-lapply(near_skcm12$step_1[which(near_skcm12$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm12[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm12_step3<-lapply(near_skcm12$step_3[which(near_skcm12$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm12[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm12_far<-lapply(c('bdy_min','dis_min'),function(x){#x='bdy_min'
  aa<-near_skcm12[,c('cell_name','FinalLocalType',x)]
  aa<-aa[aa$FinalLocalType%in%c('Normal','Immune'),]
  aa<-aa[order(aa[,3],decreasing = T),]
  aa<-aa$cell_name[1:round(nrow(aa)*0.1)]################前百分之多少
  return(aa)
}) %>% unlist() %>% unique()
skcm12_far<-setdiff(skcm12_far,skcm12_step3)



skcm13_step1<-lapply(near_skcm13$step_1[which(near_skcm13$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm13[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm13_step3<-lapply(near_skcm13$step_3[which(near_skcm13$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm13[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm13_far<-lapply(c('bdy_min','dis_min'),function(x){#x='bdy_min'
  aa<-near_skcm13[,c('cell_name','FinalLocalType',x)]
  aa<-aa[aa$FinalLocalType%in%c('Normal','Immune'),]
  aa<-aa[order(aa[,3],decreasing = T),]
  aa<-aa$cell_name[1:round(nrow(aa)*0.1)]################前百分之多少
  return(aa)
}) %>% unlist() %>% unique()
skcm13_far<-setdiff(skcm13_far,skcm13_step3)

skcm12_step1<-paste0('skcm12_',skcm12_step1)
skcm12_far<-paste0('skcm12_',skcm12_far)

skcm13_step1<-paste0('skcm13_',skcm13_step1)
skcm13_far<-paste0('skcm13_',skcm13_far)


Cyto_skcm$seurat_clusters[Cyto_skcm$SpotID%in%c(skcm12_far,skcm13_far)] %>% table()
Cyto_skcm$seurat_clusters[Cyto_skcm$SpotID%in%c(skcm13_step1,skcm12_step1)] %>% table()


data_1<-data.frame(row.names = c('cluster0','cluster1','cluster2','cluster3'),
                   far_num=c(5,22,6,0),
                   step1_num=c(240,123,19,1))

plot_1<-reshape2::melt(as.matrix(data_1))
plot_1<-mutate(plot_1,Var1 = factor(plot_1$Var1, levels = rev(unique(plot_1$Var1))))

p_compare<-ggplot(plot_1,aes(x=Var2,y=value,fill=Var1)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  scale_fill_manual(values = c('cluster0'='#E48269','cluster1'='#6C7FAA','cluster2'='#94B683','cluster3'='#FDD0B1'))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("")+ylab("proportion")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('')
print(p_compare)
dir_pic11<-'E:/Mirror/ST_analysis/pic/re/4/subTAM_distribution/'
pdf(paste0(dir_pic11,'tam_subcluster_distribution_barplot_dis_new.pdf'),height = 6,width = 4)
print(p_compare)
dev.off()

