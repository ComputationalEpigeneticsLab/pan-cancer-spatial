###
library(dplyr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(monocle)
library(igraph)
library(RobustRankAggreg)
library(jsonlite)
library(reshape2)


dir_MPdist<-'/10X Visium/new_st_MPdist/'
dir_MPdist<-'/10X_Visium/new_st_MPdist/'
file_dist<-list.files(pattern = '_MPdistance.txt',path = dir_MPdist,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_dist,'_MPdist'),function(x)x[1]))
cancer<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[1]))
dataSet<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[2]))
slices<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[3]))
grep('OVCA/ovca01/slice1',dataSlice)

all_MP_compare_dist<-c()
all_eco_compare_dist<-c()

for(i in 1:length(dataSlice)){
  # i=1
  st_dist<-read.delim(paste0(dir_MPdist,file_dist[i]),stringsAsFactors = F,check.names = F)
  st_dist$scale_dist<-st_dist$distance-min(st_dist$distance)
  st_dist$scale_dist<-st_dist$scale_dist/max(st_dist$scale_dist)
  st_dist<-st_dist[st_dist$spot1_localType%in%c('Core','Budding','Boundary'),]
  st_dist<-st_dist[st_dist$spot2_localType%in%c('Core','Budding','Boundary'),]
  
  table(st_dist$spot1_localType)
  table(st_dist$spot2_localType)
  
  region_Type_all<-intersect(st_dist$spot1_localType,st_dist$spot2_localType)
  
  st_dist$spot1_eco<-'ecosystem1'
  st_dist$spot1_eco[st_dist$spot1_MP%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
  st_dist$spot1_eco[st_dist$spot1_MP%in%c('MP_7','MP_17')]<-'ecosystem3'
  st_dist$spot1_eco[st_dist$spot1_MP%in%c('MP_1','MP_4')]<-'ecosystem4'
  st_dist$spot1_eco[st_dist$spot1_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
  st_dist$spot1_eco[st_dist$spot1_MP%in%c('MP_15','MP_11')]<-'ecosystem6'
  
  st_dist$spot2_eco<-'ecosystem1'
  st_dist$spot2_eco[st_dist$spot2_MP%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
  st_dist$spot2_eco[st_dist$spot2_MP%in%c('MP_7','MP_17')]<-'ecosystem3'
  st_dist$spot2_eco[st_dist$spot2_MP%in%c('MP_1','MP_4')]<-'ecosystem4'
  st_dist$spot2_eco[st_dist$spot2_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
  st_dist$spot2_eco[st_dist$spot2_MP%in%c('MP_15','MP_11')]<-'ecosystem6'
  
  
  aa<-unique(st_dist$spot1[which(st_dist$spot1_localType=='Core')])
  bb<-unique(st_dist$spot2[which(st_dist$spot2_localType=='Core')])
  cc<-union(aa,bb)
  
  
  
  ###每个MP内部的距离和外部的距离
  MP_compare_dist<-c()
  for(jj in region_Type_all){
    region_Type<-jj
    st_dist_region<-st_dist[st_dist$spot1_localType%in%region_Type,]
    st_dist_region<-st_dist_region[st_dist_region$spot2_localType%in%region_Type,]
    
    MP_c_dist<-lapply(unique(c(st_dist_region$spot1_MP,st_dist_region$spot2_MP)),function(x){#x='MP_1'
      if(length(which(st_dist_region$spot1_MP==x&st_dist_region$spot2_MP==x))>3){
        MP_intra<-mean(st_dist_region$scale_dist[which(st_dist_region$spot1_MP==x&st_dist_region$spot2_MP==x)])
      }else{
        MP_intra<-NA
      }
      
      if(length(which(st_dist_region$spot1_MP==x&st_dist_region$spot2_MP!=x))>3){
        MP_inter<-mean(st_dist_region$scale_dist[which(st_dist_region$spot1_MP==x&st_dist_region$spot2_MP!=x)])
      }else{
        MP_inter<-NA
      }
      
      c_dist<-c(MP_intra,MP_inter)
      names(c_dist)<-c('MP_intra','MP_inter')
      return(c_dist)
    })
    MP_c_dist<-do.call(rbind,MP_c_dist)
    rownames(MP_c_dist)<-unique(c(st_dist_region$spot1_MP,st_dist_region$spot2_MP))
    MP_c_dist<-reshape2::melt(as.matrix(MP_c_dist))
    colnames(MP_c_dist)[1:3]<-c('MP','group','scale_dist')
    MP_c_dist$dataSlice<-dataSlice[i]
    MP_c_dist$region<-region_Type
    
    MP_compare_dist<-rbind(MP_compare_dist,MP_c_dist)
  }
  
  
  
  ###每个ecosystem内部的距离和外部的距离
  eco_compare_dist<-c()
  for(jj in region_Type_all){
    region_Type<-jj
    st_dist_region<-st_dist[st_dist$spot1_localType%in%region_Type,]
    st_dist_region<-st_dist_region[st_dist_region$spot2_localType%in%region_Type,]
    
    eco_c_dist<-lapply(unique(c(st_dist_region$spot1_eco,st_dist_region$spot2_eco)),function(x){#x='ecosystem1'
      if(length(which(st_dist_region$spot1_eco==x&st_dist_region$spot2_eco==x))>3){
        eco_intra<-mean(st_dist_region$scale_dist[which(st_dist_region$spot1_eco==x&st_dist_region$spot2_eco==x)])
      }else{
        eco_intra<-NA
      }
      
      if(length(which(st_dist_region$spot1_eco==x&st_dist_region$spot2_eco!=x))>3){
        eco_inter<-mean(st_dist_region$scale_dist[which(st_dist_region$spot1_eco==x&st_dist_region$spot2_eco!=x)])
      }else{
        eco_inter<-NA
      }
      
      c_dist<-c(eco_intra,eco_inter)
      names(c_dist)<-c('eco_intra','eco_inter')
      return(c_dist)
    })
    eco_c_dist<-do.call(rbind,eco_c_dist)
    rownames(eco_c_dist)<-unique(c(st_dist_region$spot1_eco,st_dist_region$spot2_eco))
    eco_c_dist<-reshape2::melt(as.matrix(eco_c_dist))
    colnames(eco_c_dist)[1:3]<-c('eco','group','scale_dist')
    eco_c_dist$dataSlice<-dataSlice[i]
    eco_c_dist$region<-region_Type
    
    eco_compare_dist<-rbind(eco_compare_dist,eco_c_dist)
  }
  
  all_MP_compare_dist<-rbind(all_MP_compare_dist,MP_compare_dist)
  all_eco_compare_dist<-rbind(all_eco_compare_dist,eco_compare_dist)
  print(dataSlice[i])
}
write.table(all_MP_compare_dist,paste0(dir_MPdist,'00data/all_MP_compare_dist.txt'),
            quote = F,sep = '\t',row.names = F)
write.table(all_eco_compare_dist,paste0(dir_MPdist,'00data/all_eco_compare_dist.txt'),
            quote = F,sep = '\t',row.names = F)




###可视化
all_MP_compare_dist<-read.delim('/10X Visium/new_st_MPdist/00data/all_MP_compare_dist.txt',
                                stringsAsFactors = F,check.names = F)
all_eco_compare_dist<-read.delim('/10X Visium/new_st_MPdist/00data/all_eco_compare_dist.txt',
                                 stringsAsFactors = F,check.names = F)

dir_pic<-'/10X Visium/new_st_MPdist/00pic/'

all_MP_compare_dist<-all_MP_compare_dist[!is.na(all_MP_compare_dist$scale_dist),]

pdf(paste0(dir_pic,'all_region_MP.pdf'),width = 8,height = 3)
pp_data<-all_MP_compare_dist
colnames(pp_data)
table(pp_data$group)
pp_data$group<-factor(pp_data$group,levels = c("MP_intra","MP_inter"))
pp_data$MP<-factor(pp_data$MP,levels = paste0('MP_',1:17))
range(pp_data$scale_dist)

p_box<-ggplot(pp_data, aes(x = MP, y = scale_dist,fill=group))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = group), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("MP_intra"="#99CCCC","MP_inter"="#336699"))+
  scale_color_manual(values = c("MP_intra"="#99CCCC","MP_inter"="#336699"))+
  ggtitle(paste0('MP_',"all_region"))+
  # theme(plot.title = element_text(hjust = 0.4))+
  # theme(plot.title = element_text(size = 6))+
  #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
  stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
  #theme_bw()+
  ylim(c(0,0.7))+
  theme(axis.title.x = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(p_box)

for(i in c('Core','Boundary','Budding')){
  pp_data<-all_MP_compare_dist[all_MP_compare_dist$region%in%i,]
  colnames(pp_data)
  table(pp_data$group)
  pp_data$group<-factor(pp_data$group,levels = c("MP_intra","MP_inter"))
  pp_data$MP<-factor(pp_data$MP,levels = paste0('MP_',1:17))
  range(pp_data$scale_dist)
  
  p_box<-ggplot(pp_data, aes(x = MP, y = scale_dist,fill=group))+ 
    stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = group), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("MP_intra"="#99CCCC","MP_inter"="#336699"))+
    scale_color_manual(values = c("MP_intra"="#99CCCC","MP_inter"="#336699"))+
    ggtitle(paste0('MP_',i))+
    stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
    ylim(c(0,0.7))+
    theme(axis.title.x = element_text(size=12),
          axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  print(p_box)
}
dev.off()



all_eco_compare_dist<-all_eco_compare_dist[!is.na(all_eco_compare_dist$scale_dist),]

pdf(paste0(dir_pic,'all_region_eco.pdf'),width = 4,height = 3)
pp_data<-all_eco_compare_dist
colnames(pp_data)
table(pp_data$group)
pp_data$group<-factor(pp_data$group,levels = c("eco_intra","eco_inter"))
pp_data$eco<-factor(pp_data$eco,levels = paste0('ecosystem',1:6))
range(pp_data$scale_dist)

p_box<-ggplot(pp_data, aes(x = eco, y = scale_dist,fill=group))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = group), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("eco_intra"="#99CCCC","eco_inter"="#336699"))+
  scale_color_manual(values = c("eco_intra"="#99CCCC","eco_inter"="#336699"))+
  ggtitle(paste0('eco_',"all_region"))+
  # theme(plot.title = element_text(hjust = 0.4))+
  # theme(plot.title = element_text(size = 6))+
  #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
  stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
  #theme_bw()+
  ylim(c(0,0.7))+
  theme(axis.title.x = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(p_box)

for(i in c('Core','Boundary','Budding')){
  pp_data<-all_eco_compare_dist[all_eco_compare_dist$region%in%i,]
  colnames(pp_data)
  table(pp_data$group)
  pp_data$group<-factor(pp_data$group,levels = c("eco_intra","eco_inter"))
  pp_data$eco<-factor(pp_data$eco,levels = paste0('ecosystem',1:6))
  range(pp_data$scale_dist)
  
  p_box<-ggplot(pp_data, aes(x = eco, y = scale_dist,fill=group))+ 
    stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = group), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("eco_intra"="#99CCCC","eco_inter"="#336699"))+
    scale_color_manual(values = c("eco_intra"="#99CCCC","eco_inter"="#336699"))+
    ggtitle(paste0('eco_',i))+
    stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
    ylim(c(0,0.7))+
    theme(axis.title.x = element_text(size=12),
          axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  print(p_box)
}
dev.off()






