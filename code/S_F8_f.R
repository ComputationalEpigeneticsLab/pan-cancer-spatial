#####skcm的TAM和CAF在bdy和dis step_1与far中的点上含量高低
###散点小提琴图表示
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


###RCTD
skcm12_subRCTD<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/RCTD_4/skcm12/slice1_Deconvolution.txt',stringsAsFactors = F,check.names = F)
skcm13_subRCTD<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/RCTD_4/skcm13/slice1_Deconvolution.txt',stringsAsFactors = F,check.names = F)

skcm12_subRCTD<-skcm12_subRCTD[,which(colnames(skcm12_subRCTD)!='unknown')]
skcm13_subRCTD<-skcm13_subRCTD[,which(colnames(skcm13_subRCTD)!='unknown')]

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
skcm12_data<-data.frame(slice='skcm12',
                        distance=c(rep('step1',length(skcm12_step1)),rep('far',length(skcm12_far))),
                        TAM_C0=c(skcm12_subRCTD[skcm12_step1,'TAM_C0'],skcm12_subRCTD[skcm12_far,'TAM_C0']),
                        TAM_C1=c(skcm12_subRCTD[skcm12_step1,'TAM_C1'],skcm12_subRCTD[skcm12_far,'TAM_C1']),
                        TAM_C2=c(skcm12_subRCTD[skcm12_step1,'TAM_C2'],skcm12_subRCTD[skcm12_far,'TAM_C2']),
                        CAF=c(skcm12_subRCTD[skcm12_step1,'CAF'],skcm12_subRCTD[skcm12_far,'CAF']))

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
skcm13_data<-data.frame(slice='skcm13',
                        distance=c(rep('step1',length(skcm13_step1)),rep('far',length(skcm13_far))),
                        TAM_C0=c(skcm13_subRCTD[skcm13_step1,'TAM_C0'],skcm13_subRCTD[skcm13_far,'TAM_C0']),
                        TAM_C1=c(skcm13_subRCTD[skcm13_step1,'TAM_C1'],skcm13_subRCTD[skcm13_far,'TAM_C1']),
                        TAM_C2=c(skcm13_subRCTD[skcm13_step1,'TAM_C2'],skcm13_subRCTD[skcm13_far,'TAM_C2']),
                        CAF=c(skcm13_subRCTD[skcm13_step1,'CAF'],skcm13_subRCTD[skcm13_far,'CAF']))

plot_2<-rbind(skcm12_data,skcm13_data)
plot_2$TAM<-apply(plot_2[,c("TAM_C0","TAM_C1","TAM_C2")],1,sum)
# plot_2$localType<-unlist(lapply(strsplit(plot_2$distance,'_'),function(x)x[1]))
# plot_2$near<-unlist(lapply(strsplit(plot_2$distance,'_'),function(x)x[2]))
# plot_2$slice_Type<-paste0(plot_2$slice,'_',plot_2$localType)

colnames(plot_2)
plot_col<-c("TAM_C0")

#dir_pic22<-'E:/Mirror/ST_analysis/pic/subTAM&CAF_skcm/RCTD_box/'
dir_pic22<-'E:/Mirror/ST_analysis/pic/re/other/RCTD_TAM_CAF_new/'

for(i in 1:length(plot_col)){#i=1
  #plot_pp<-plot_2[,c('slice_Type','near',plot_col[i])]
  plot_pp<-plot_2[,c('slice','distance',plot_col[i])]
  colnames(plot_pp)[3]<-'celltype'
  
  e <- ggplot(plot_pp, aes(x = slice, y = celltype,fill=distance))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    # geom_boxplot(aes(color = near), width = 0.5,position = position_dodge(0.8),alpha=0.5)+
    #geom_sina(size = 1)+
    geom_jitter(aes(color = distance),size = 1,
                position=position_jitterdodge(jitter.width = 0.2, 
                                              jitter.height = 0, 
                                              dodge.width = 0.8)) + # 不重叠的散点图
    stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8),color="black") + # 误差棒，中位数，25%和75%分位数
    stat_summary(aes(fill = distance), fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8),color="black") + # 中位数水平线
    #facet_wrap(~ slice_Type, scales = 'free_y', nrow = 1) +####癌症分组
    theme_classic(base_size = 20)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("step1"="#827CBA","far"="#BBB7D8"))+
    scale_color_manual(values = c("step1"="#827CBA","far"="#BBB7D8"))+
    ggtitle(plot_col[i])+
    theme(plot.title = element_text(hjust = 0.4))+
    theme(plot.title = element_text(size = 12))+
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =max(plot_pp$celltype))+
    theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  pdf(paste0(dir_pic22,plot_col[i],'_bdy_dis.pdf'),width = 6,height = 5)
  print(e)
  dev.off()
}



