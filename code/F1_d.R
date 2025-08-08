####spot标签改名字
library(tidyverse)
library(ggplot2)
library(ggpubr)

dir_bdy<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"
file_bdy<-list.files(pattern = "_BdyTumorCore.txt",path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_bdy,"_"),function(x) x[1]))

pdf("E:/Mirror/ST_analysis/pic/copykat/TumorCoreBdy_4.pdf",width = 5, height = 4)
for(i in 1:length(file_bdy)){
  #i=168
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  #st_bdy<-st_nearSpot
  # table(st_bdy$LocationType)
  # st_bdy$FinalLocalType<-st_bdy$LocationType
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="Tumor")]<-"Dispersion"
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="diploid")]<-"Normal"
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
  table(st_bdy$FinalLocalType)
  #write.table(st_bdy,paste0(dir_bdy,file_bdy[i]),quote = F,sep = '\t')
  
  p_spot2<-ggplot(data = st_bdy,aes(x = imagerow, y = imagecol,color=FinalLocalType)) + 
    geom_point(size=1)+
    # scale_color_manual(values =  c("Dispersion"="#CD5A5A","Boundary"="#FF9900","Normal"="#669933","Core"="#990033",
    #                                "Immune"="#5477AF","not.defined"="#CCCCCC"))+
    scale_color_manual(values =  c('Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d',
                                   'Immune'='#a9d38a','Normal'='#2375ae'))+
    theme_classic()+
    ggtitle(dataset_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=2)))
  print(p_spot2)
  # pdf('E:/Mirror/ST_analysis/pic/re/4/bdy_core_dis.pdf',width = 4.6, height = 4)
  # print(p_spot2)
  # dev.off()
  
}
dev.off()

write.table(dataset_slice,'E:/Mirror/ST_analysis/data/dataset_slice_order.txt',quote = F,sep = '\t',row.names = F)


pdf("E:/Mirror/ST_analysis/pic/re/1/TumorCoreBdy_3.pdf",width = 5, height = 4)
print(p_spot2)
dev.off()


dir_bdy<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"
select_slice<-c('brca01/slice3','luad03/slice1',
                'brca08/slice4','brca12/slice1',
                'gbm04/slice7','lihc02/slice1',
                'lihc02/slice12','lihc02/slice14',
                'lihc03/slice6','luad01/slice2',
                'ovca03/slice1','pdac03/slice1')
file_bdy<-paste0(select_slice,'_BdyTumorCore.txt')
dataset_slice<-unlist(lapply(strsplit(file_bdy,"_"),function(x) x[1]))

pdf("E:/Mirror/ST_analysis/pic/copykat/select_CAFTAM.pdf",width = 5, height = 4)
for(i in 1:length(file_bdy)){
  #i=168
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  #st_bdy<-st_nearSpot
  # table(st_bdy$LocationType)
  # st_bdy$FinalLocalType<-st_bdy$LocationType
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="Tumor")]<-"Dispersion"
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="diploid")]<-"Normal"
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
  table(st_bdy$FinalLocalType)
  #write.table(st_bdy,paste0(dir_bdy,file_bdy[i]),quote = F,sep = '\t')
  
  p_spot2<-ggplot(data = st_bdy,aes(x = imagerow, y = imagecol,color=FinalLocalType)) + 
    geom_point(size=0.6)+
    # scale_color_manual(values =  c("Dispersion"="#CD5A5A","Boundary"="#FF9900","Normal"="#669933","Core"="#990033",
    #                                "Immune"="#5477AF","not.defined"="#CCCCCC"))+
    scale_color_manual(values =  c('Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d',
                                   'Immune'='#a9d38a','Normal'='#2375ae'))+
    theme_classic()+
    ggtitle(dataset_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=2)))
  print(p_spot2)
  # pdf('E:/Mirror/ST_analysis/pic/re/4/bdy_core_dis.pdf',width = 4.6, height = 4)
  # print(p_spot2)
  # dev.off()
  
}
dev.off()

