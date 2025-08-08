####MP 分5类的各个免疫通路得分绘制
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)

####各种基因集得分 对MP聚类的5个簇可视化热图
##使用恶性spot的均值，再对每个MP类的所有切片算均值
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gghalves)
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


###转移相关基因集
geneset1<-read.delim('E:/Mirror/ST_analysis/data/geneset/genest_mp.txt',stringsAsFactors = F,check.names = F)
all_geneset<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(all_geneset)<-geneset1$pathway
####可以用AddModuleScore计算的
all_geneset_use<-all_geneset

dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dir_EMT<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'

for(i in 1:length(file_bdy)){
  #i=2
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  
  st_rds<-st_rds[,rownames(st_bdy)]
  
  all_geneset_use<-all_geneset
  inter_gene<-lapply(all_geneset,function(x) length(intersect(x,rownames(st_rds))))%>%unlist()
  if(length(which(inter_gene<=0))>0){
    for(j in 1:length(which(inter_gene<=0))){
      all_geneset_use[[which(inter_gene<=0)[j]]]<-rownames(st_rds)[1:3]
    }
  }
  
  st_rds<-AddModuleScore(st_rds,all_geneset_use,name = names(all_geneset_use))
  score_re<-st_rds@meta.data
  score_re<-score_re[,paste0(names(all_geneset_use),1:length(names(all_geneset_use)))]
  colnames(score_re)<-names(all_geneset_use)
  score_re[,which(inter_gene<=0)]<-NA######没有交集的基因集得分为NA
  
  score_re$LocalType<-st_bdy$FinalLocalType
  score_re$imagerow<-st_bdy$imagerow
  score_re$imagecol<-st_bdy$imagecol
  score_re$cell_name<-rownames(score_re)
  
  write.table(score_re,paste0(dir_EMT,dataset_slice[i],'_MP5_score.txt'),quote = F,sep = '\t')
  print(dataset_slice[i])
}



dir_pic<-'E:/Mirror/ST_analysis/pic/NMF_module/Imm_path/'
dir_EMT<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'
file_metas<-list.files(pattern = '_MP5_score.txt',path = dir_EMT,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_metas,'_'),function(x) x[1]))

pathway<-colnames(st_MPscore)[1:18]

for(j in 1:length(pathway)){
  #j=1
  pdf(paste0(dir_pic,pathway[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(file_MPscore)){
    #i=1
    st_MPscore<-read.delim(paste0(dir_EMT,file_metas[i]),stringsAsFactors = F,check.names = F)
    
    plot_data<-data.frame(imagerow=st_MPscore$imagerow,imagecol=st_MPscore$imagecol,
                          MPscore=st_MPscore[,pathway[j]])
    plot_data$MPscore<-plot_data$MPscore-min(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    plot_data$MPscore<-plot_data$MPscore/max(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    range(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    plot_data$MPscore[which(st_MPscore$LocalType!='Core')]<--0.001
    p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
      geom_point(aes(colour=MPscore),size=.6) +
      # scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
      #                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
      scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(1),
                                        colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                        colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                        colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
      )+
      theme_classic()+
      labs(title = paste0(dataset_slice[i],'_of_',pathway[j]),x = "",y = "",col='score')
    print(p)
  }
  dev.off()
}
















