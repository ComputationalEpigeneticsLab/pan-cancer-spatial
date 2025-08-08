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


dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_MP<-list.files(pattern = '_MP_score.txt',recursive = T,path = dir_bdy)
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
patient<-unlist(lapply(strsplit(file_MP,'/'),function(x)x[1]))


for(i in 1:2){#i=1
  MP_score<-read.delim(paste0(dir_bdy,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  spot_mal<-c('Core','Boundary','Budding')
  st_bdy<-st_bdy[rownames(MP_score),]
  
  pdf(paste0('E:/Mirror/ST_analysis/pic/ESCC/MPscore/',patient[i],'_spot_MPscore2.pdf'),width = 4.8, height = 4)
  for(j in 1:13){#j=3
    p_score<-data.frame(score=MP_score[,j],
                        localType=st_bdy$FinalLocalType,
                        MP_top1=MP_score$MP_top1,
                        imagerow=st_bdy$imagerow,
                        imagecol=st_bdy$imagecol)
    p_score$score[!p_score$localType%in%spot_mal]<-(-0.001)
    p_score$score[!p_score$MP_top1%in%colnames(MP_score)[j]]<-(-0.001)
    p_score$MP_type<-'other'
    p_score$MP_type[which(p_score$score>=0)]<-colnames(MP_score)[j]
    
    p<-ggplot(p_score, aes(x=imagerow, y=imagecol)) +
      geom_point(aes(colour=score),size=1) +
      scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(1),
                                        colorRampPalette(c("#5b2177","#b93c6d"))(5),
                                        colorRampPalette(c("#b93c6d","#eb7d60"))(40),
                                        colorRampPalette(c("#eb7d60","#f6f0b7"))(40))
      )+
      theme_classic()+
      labs(title = paste0(patient[i],'_',colnames(MP_score)[j]),x = "",y = "",col='score')
    print(p)
    
    # p_spot2<-ggplot(data = p_score,aes(x = imagerow, y = imagecol,color=MP_type)) + 
    #   geom_point(size=1)+
    #   scale_color_manual(values =  c('#d62d28','grey80'))+
    #   theme_classic()+
    #   labs(title = paste0(patient[i],'_',colnames(MP_score)[j]),x = "",y = "",col='score')+
    #   guides(colour = guide_legend(override.aes = list(size=2)))
    # print(p_spot2)
    
  }
  dev.off()
}
