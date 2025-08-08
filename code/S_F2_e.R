#####B细胞亲和力 和 T细胞耗竭得分与core占比的关系
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(psych)


B_avidity<-c('IGHG1','IGHG2','IGHG3','IGHG4','IGHM','cc11')
T_exhausted<-read.delim('E:/Mirror/ST_analysis/data/geneset/other.txt',stringsAsFactors = F)
T_exhausted<-T_exhausted[6:7,]
name_exh<-T_exhausted$geneset_name
T_exhausted<-lapply(T_exhausted$symbol,function(x)unlist(strsplit(x,',')))
names(T_exhausted)<-name_exh

# load('E:/Mirror/ST_analysis/program/TCellSI-main/data/markers.rda')
# TcellSI_exh<-markers[8:9]
# T_exhausted<-c(T_exhausted,TcellSI_exh)


dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
table(unlist(lapply(strsplit(file_bdy,'_'),function(x) x[1]))==dataset_slice)

dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD3/'
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
table(unlist(lapply(strsplit(file_RCTD,'_'),function(x) x[1]))==dataset_slice)

all_score_re<-c()
for(i in 1:229){
  #i=2
  # st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[rownames(st_RCTD),]
  # st_RCTD<-st_RCTD[colnames(st_rds),]
  st_RCTD$cellType<-apply(st_RCTD,1,function(x){#x=st_RCTD[1,]
    names(x)<-colnames(st_RCTD)
    id<-which(x>0.5)
    aa<-'no'
    if(length(id)>=1){
      aa<-names(x)[which.max(x)]
    } 
    return(aa)
  })
  table(st_RCTD$cellType)
  st_RCTD$localType<-st_bdy$FinalLocalType
  #st_RCTD<-st_RCTD[which(st_RCTD$localType!='not.defined'),]
  
  
  # all_geneset_use<-T_exhausted
  # inter_gene<-lapply(T_exhausted,function(x) length(intersect(x,rownames(st_rds))))%>%unlist()
  # if(length(which(inter_gene<=0))>0){
  #   for(j in 1:length(which(inter_gene<=0))){
  #     all_geneset_use[[which(inter_gene<=0)[j]]]<-rownames(st_rds)[1:3]
  #   }
  # }
  # st_rds<-AddModuleScore(st_rds,all_geneset_use,name = names(all_geneset_use))
  # score_re<-st_rds@meta.data
  # score_re<-score_re[,paste0(names(all_geneset_use),1:length(names(all_geneset_use)))]
  # colnames(score_re)<-names(all_geneset_use)
  # 
  # 
  # B_score<-as.matrix(st_rds@assays[["Spatial"]]@counts) %>% as.data.frame()
  # B_score<-B_score[B_avidity,] %>% as.matrix()
  # rownames(B_score)<-B_avidity
  # B_score[is.na(B_score)]<-0
  # B_score<-apply(B_score,2,function(x){
  #   aa<-sum(x[1:4])-x[5]
  #   return(aa)
  # })
  # 
  # score_re$B_avidity<-B_score
  # score_re<-cbind(score_re,st_RCTD[,c('cellType','localType')])
  # score_re$cell_name<-rownames(score_re)
  # score_re$slice<-dataset_slice[i]
  st_RCTD$cell_name<-rownames(st_RCTD)
  st_RCTD$slice<-dataset_slice[i]
  
  all_score_re<-rbind(all_score_re,st_RCTD[,c('cellType','localType','cell_name','slice')])
  
  print(dataset_slice[i])
}
write.table(all_score_re,'E:/Mirror/ST_analysis/data/pic_data/T_exhausted_B_avidity.txt',quote = F,sep = '\t',row.names = F)

T_exhausted_B_avidity<-read.delim('E:/Mirror/ST_analysis/data/pic_data/T_exhausted_B_avidity.txt',
                                  stringsAsFactors = F,check.names = F)
table(paste0(all_score_re$slice,all_score_re$cell_name)==paste0(T_exhausted_B_avidity$slice,T_exhausted_B_avidity$cell_name))
T_exhausted_B_avidity$cellType2<-all_score_re$cellType
write.table(T_exhausted_B_avidity,'E:/Mirror/ST_analysis/data/pic_data/T_exhausted_B_avidity.txt',quote = F,sep = '\t',row.names = F)

T_exhausted_B_avidity$N_M<-T_exhausted_B_avidity$localType
T_exhausted_B_avidity$N_M[T_exhausted_B_avidity$localType%in%c('Normal','Immune')]<-'Nor_Imm'

dataSlice<-unique(T_exhausted_B_avidity$slice)

core_ratio<-c()
T_score1<-c()
T_score2<-c()
B_score<-c()
for(i in 1:length(dataSlice)){
  #i=1
  slice_score<-T_exhausted_B_avidity[T_exhausted_B_avidity$slice%in%dataSlice[i],]
  slice_score<-slice_score[which(slice_score$localType!='not.defined'),]
  table(slice_score$cellType2)
  
  core_ratio<-c(core_ratio,length(which(slice_score$localType=='Core'))/nrow(slice_score))
  
  T_site<-which(slice_score$cellType2=='T lymphocytes'&slice_score$N_M=='Nor_Imm')
  if(length(T_site)>2){
    T_score1<-c(T_score1,mean(slice_score$T_exhausted_Zhangcl[T_site]))
    T_score2<-c(T_score2,mean(slice_score$T_exhausted_Xug[T_site]))
  }else{
    T_score1<-c(T_score1,0)
    T_score2<-c(T_score2,0)
  }
  
  B_site<-which(slice_score$cellType2=='B lymphocytes'&slice_score$N_M=='Nor_Imm')
  if(length(B_site)>2){
    B_score<-c(B_score,mean(slice_score$B_avidity[B_site]))
  }else{
    B_score<-c(B_score,0)
  }
}

slice_re<-data.frame(slice=dataSlice,
                     T_exhausted_Zhangcl=T_score1,
                     T_exhausted_Xug=T_score2,
                     B_avidity=B_score,
                     core_ratio=core_ratio)

slice_re<-slice_re[which(slice_re$core_ratio!=0),]

dir_picc<-'E:/Mirror/ST_analysis/pic/T_exhaust_B_avidity/'

plot_data<-data.frame(core_ratio=rep(slice_re$core_ratio,2),
                      T_exhausted=c(slice_re$T_exhausted_Zhangcl,slice_re$T_exhausted_Xug),
                      exhausted_type=rep(c('T_exhausted_Zhangcl','T_exhausted_Xug'),each=nrow(slice_re)))
plot_data<-plot_data[which(plot_data$T_exhausted!=0),]
plot_data<-plot_data[which(plot_data$exhausted_type=='T_exhausted_Xug'),]

pdf(paste0(dir_picc,'T_exhausted_core2.pdf'),width = 5.5,height = 3.5)
b <- ggplot(plot_data, aes(x = core_ratio, y = T_exhausted))
pp<-b + geom_point(aes(color = exhausted_type),size=0.01)+
  geom_smooth(aes(color = exhausted_type, fill = exhausted_type), method = "lm",se=T) +
  #geom_rug(aes(color =exhausted_type)) +
  scale_color_manual(values = c("#53b654", "#16afe2"))+
  scale_fill_manual(values = c("#53b654", "#16afe2"))+
  ggpubr::stat_cor(aes(color = exhausted_type), label.x = 0.1)+
  theme_classic()+
  ggtitle('T_exhausted')
print(pp)
dev.off()
?stat_cor


plot_data<-data.frame(core_ratio=rep(slice_re$core_ratio,1),
                      avidity=c(slice_re$B_avidity),
                      B_avidity=rep(c('B_avidity'),each=nrow(slice_re)))
plot_data<-plot_data[which(plot_data$avidity!=0),]

pdf(paste0(dir_picc,'B_avidity_core.pdf'),width = 5,height = 3.5)
b <- ggplot(plot_data, aes(x = core_ratio, y = avidity))
pp<-b + geom_point(aes(color = B_avidity),size=0.01)+
  geom_smooth(aes(color = B_avidity, fill = B_avidity), method = "lm",se=T) +
  #geom_rug(aes(color =B_avidity)) +
  scale_color_manual(values = c("#f08634"))+
  scale_fill_manual(values = c("#f08634"))+
  ggpubr::stat_cor(aes(color = B_avidity), label.x = 0.1)+
  theme_classic()+
  ggtitle('B_avidity')
print(pp)
dev.off()



plot_data<-data.frame(avidity=rep(slice_re$B_avidity,2),
                      T_exhausted=c(slice_re$T_exhausted_Zhangcl,slice_re$T_exhausted_Xug),
                      exhausted_type=rep(c('T_exhausted_Zhangcl','T_exhausted_Xug'),each=nrow(slice_re)))
plot_data<-plot_data[which(plot_data$T_exhausted!=0&plot_data$avidity!=0),]
plot_data<-plot_data[which(plot_data$exhausted_type=='T_exhausted_Xug'),]

pdf(paste0(dir_picc,'T_exhausted_B_avidity2.pdf'),width = 5.5,height = 3.5)
b <- ggplot(plot_data, aes(x = avidity, y = T_exhausted))
pp<-b + geom_point(aes(color = exhausted_type),size=0.01)+
  geom_smooth(aes(color = exhausted_type, fill = exhausted_type), method = "lm",se=T) +
  #geom_rug(aes(color =exhausted_type)) +
  scale_color_manual(values = c("#53b654", "#16afe2"))+
  scale_fill_manual(values = c("#53b654", "#16afe2"))+
  ggpubr::stat_cor(aes(color = exhausted_type), label.x = 0.1)+
  theme_classic()+
  ggtitle('T_exhausted')
print(pp)
dev.off()





