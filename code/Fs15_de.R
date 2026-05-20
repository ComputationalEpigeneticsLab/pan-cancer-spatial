###MP得分
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
source('/moduleScoreToCell.R')
dir_rds<-'/new_st/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('.rds','',file_rds)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))
dir_mp<-'/MPscore/MPscore/'

# MP_list<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/MP_list_intersect_initia25_intersect_cluster24.txt',
#                     stringsAsFactors = F,check.names = F)
# MP_list<-MP_list[,-1]
# modules<-lapply(1:ncol(MP_list),function(x) MP_list[,x])
# names(modules)<-colnames(MP_list)


dir_mp <- "\\MPscore\\all_MPscore\\"
file_MP<-list.files(pattern = '_MP_score.txt',recursive = T,path = dir_mp)
pdf('\\MPscore\\all_MPscore\\assign_MP_public.pdf',width = 7,height = 7)
for(i in c(1,3,4)){
  #i=2
  MP_score<-read.delim(paste0(dir_mp,file_MP[i]),stringsAsFactors = F,check.names = F)

  plot_data<-as.data.frame.array(table(MP_score[MP_score$st_bdy.FinalLocalType%in%spot_mal,'MP_top1']))
  plot_data$assign_MP<-rownames(plot_data)
  colnames(plot_data)[1]<-'num'
  plot_data$value<-plot_data$num/sum(plot_data$num)
  #plot_data<-plot_data[order(plot_data$value,decreasing = T),]
  plot_data$assign_MP<-factor(plot_data$assign_MP, levels = c(paste0("MP_",1:17)))
  plot_data<-plot_data[order(plot_data$assign_MP,decreasing = T),]
  plot_data$ymax<-cumsum(plot_data$value)
  plot_data$ymin<-c(0,head(plot_data$ymax,n=-1))
  labelPosition<-(plot_data$ymax + plot_data$ymin)/2
  
  p5 <- ggplot(plot_data,aes(x = 1, y = value, fill = assign_MP)) +
    geom_col(colour = "white")+ 
    ggtitle(file_MP[i])+
    coord_polar(theta = "y", start = 1.65) +
    geom_text(aes(label = paste0(round(value * 100, 2), "%"),x=1.6,y=labelPosition),
              #position = position_fill(vjust = 0.5),
              #hjust = 3,vjust = 0.5,
              size=3) +
    scale_fill_manual(values=c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1","MP_5"="#6A8EC9",
                               "MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861","MP_9"="#652884","MP_10"="#444577",
                               "MP_11"="#8A7355","MP_12"="#B3BB61","MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
                               "MP_16"="#399335",'MP_17'='#96DD88'))+
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

