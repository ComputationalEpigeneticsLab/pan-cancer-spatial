###MP得分
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)

####文章提供的函数
source('E:/Mirror/ST_analysis/program/moduleScoreToCell.R')


dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))
dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'

MP_list<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/MP_list_intersect_initia25_intersect_cluster24.txt',
                    stringsAsFactors = F,check.names = F)
MP_list<-MP_list[,-1]
modules<-lapply(1:ncol(MP_list),function(x) MP_list[,x])
names(modules)<-colnames(MP_list)


for(i in 1:2){
  #i=1
  #if(!dir.exists(paste0(dir_moduleScore,dataset[i]))) dir.create(paste0(dir_moduleScore,dataset[i]))
  data_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  #data = GetData(data_rds, slot = 'data')
  
  srt = GeneToEnrichment(srt=data_rds, db = modules)
  #range(enrichment.profile[7,])
  MP_score<-srt@meta.data
  MP_score<-MP_score[,names(modules)]
  
  write.table(MP_score,paste0(dir_bdy,dataset_slice[i],'_MP_score.txt'),quote = F,sep = '\t')
  print(dataset_slice[i])
}


dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_MP<-list.files(pattern = '_MP_score.txt',recursive = T,path = dir_bdy)
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
patient<-unlist(lapply(strsplit(file_MP,'/'),function(x)x[1]))



pdf('E:/Mirror/ST_analysis/pic/ESCC/MPscore/assign_MP.pdf',width = 7,height = 7)
for(i in 1:2){#i=1
  MP_score<-read.delim(paste0(dir_bdy,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  spot_mal<-c('Core','Boundary','Dispersion')
  table(MP_score$MP_top1)
  # MP_score$MP95<-apply(MP_score,1,function(x){#x=MP_score[1,]
  #   names(x)<-colnames(MP_score)
  #   x<-x[order(x,decreasing = T)]
  #   y<-'unresolved'
  #   if(x[1]*0.95>=x[2]) y<-names(x)[1]
  #   return(y)
  # })
  # MP_score$MP_top1<-apply(MP_score[,1:13],1,function(x){
  #   x<-x[order(x,decreasing = T)]
  #   return(names(x)[1])
  # })
  # write.table(MP_score,paste0(dir_bdy,file_MP[i]),quote = F,sep = '\t')
  plot_data<-as.data.frame.array(table(MP_score[rownames(st_bdy)[st_bdy$FinalLocalType%in%spot_mal],'MP_top1']))
  plot_data$assign_MP95<-rownames(plot_data)
  colnames(plot_data)[1]<-'num'
  plot_data$value<-plot_data$num/sum(plot_data$num)
  #plot_data<-plot_data[order(plot_data$value,decreasing = T),]
  plot_data$assign_MP95<-factor(plot_data$assign_MP95, levels = c(paste0("MP_",2:14)))
  plot_data<-plot_data[order(plot_data$assign_MP95,decreasing = T),]
  plot_data$ymax<-cumsum(plot_data$value)
  plot_data$ymin<-c(0,head(plot_data$ymax,n=-1))
  labelPosition<-(plot_data$ymax + plot_data$ymin)/2
  
  p5 <- ggplot(plot_data,aes(x = 1, y = value, fill = assign_MP95)) +
    geom_col(colour = "white")+ 
    ggtitle(patient[i])+
    coord_polar(theta = "y", start = 1.65) +
    geom_text(aes(label = paste0(round(value * 100, 2), "%"),x=1.6,y=labelPosition),
              #position = position_fill(vjust = 0.5),
              #hjust = 3,vjust = 0.5,
              size=3) +
    scale_fill_manual(values=c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
                               "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
                               "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
                               'MP_14'='#FF9694','unresolved'='#7F7F7F'))+
    xlim(c(-0.2, 2)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  
  print(p5)
  
}
dev.off()




