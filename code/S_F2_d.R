####每个切片的恶性spot的TLS得分与干性得分算相关
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



dir_TLS<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'
file_TLS<-list.files(pattern = '_TLS_score.txt',path = dir_TLS,recursive = T)
dataSlice<-gsub('_TLS_score.txt','',file_TLS)

all_stem<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_stem.txt',stringsAsFactors = F,check.names = F)


all_slice_cor<-c()

for(i in 1:length(dataSlice)){
  #i=2
  stem_data<-all_stem[all_stem$dataset_slice%in%dataSlice[i],]
  rownames(stem_data)<-stem_data$cell_name
  TLS_data<-read.delim(paste0(dir_TLS,file_TLS[i]),stringsAsFactors = F,check.names = F)
  TLS_data<-TLS_data[intersect(rownames(TLS_data),rownames(stem_data)),]
  stem_data<-stem_data[intersect(rownames(TLS_data),rownames(stem_data)),]
  table(TLS_data$FinalLocalType)
  Mal_site<-rownames(TLS_data)[TLS_data$FinalLocalType%in%c('Core','Boundary','Dispersion')]
  Core_site<-rownames(TLS_data)[TLS_data$FinalLocalType%in%c('Core')]
  Bdy_site<-rownames(TLS_data)[TLS_data$FinalLocalType%in%c('Boundary')]
  Dis_site<-rownames(TLS_data)[TLS_data$FinalLocalType%in%c('Dispersion')]
  
  
  if(length(Mal_site)>2){
    Mal_cor<-apply(TLS_data[Mal_site,1:3],2,function(x){#x=TLS_data[Mal_site,2]
      cor.test(stem_data[Mal_site,'CytoTRACE'],x,method = 'spearman')[["estimate"]][["rho"]]
    })
    Mal_corP<-apply(TLS_data[Mal_site,1:3],2,function(x){
      cor.test(stem_data[Mal_site,'CytoTRACE'],x,method = 'spearman')[["p.value"]]
    })
  }else{
    Mal_cor<-rep(1000,3)
    Mal_corP<-rep(1000,3)
  }
  
  if(length(Core_site)>2){
    Core_cor<-apply(TLS_data[Core_site,1:3],2,function(x){
      cor.test(stem_data[Core_site,'CytoTRACE'],x,method = 'spearman')[["estimate"]][["rho"]]
    })
    Core_corP<-apply(TLS_data[Core_site,1:3],2,function(x){
      cor.test(stem_data[Core_site,'CytoTRACE'],x,method = 'spearman')[["p.value"]]
    })
  }else{
    Core_cor<-rep(1000,3)
    Core_corP<-rep(1000,3)
  }
  
  if(length(Bdy_site)>2){
    Bdy_cor<-apply(TLS_data[Bdy_site,1:3],2,function(x){
      cor.test(stem_data[Bdy_site,'CytoTRACE'],x,method = 'spearman')[["estimate"]][["rho"]]
    })
    Bdy_corP<-apply(TLS_data[Bdy_site,1:3],2,function(x){
      cor.test(stem_data[Bdy_site,'CytoTRACE'],x,method = 'spearman')[["p.value"]]
    })
  }else{
    Bdy_cor<-rep(1000,3)
    Bdy_corP<-rep(1000,3)
  }
  
  if(length(Dis_site)>2){
    Dis_cor<-apply(TLS_data[Dis_site,1:3],2,function(x){
      cor.test(stem_data[Dis_site,'CytoTRACE'],x,method = 'spearman')[["estimate"]][["rho"]]
    })
    Dis_corP<-apply(TLS_data[Dis_site,1:3],2,function(x){
      cor.test(stem_data[Dis_site,'CytoTRACE'],x,method = 'spearman')[["p.value"]]
    })
  }else{
    Dis_cor<-rep(1000,3)
    Dis_corP<-rep(1000,3)
  }
  
  slice_cor<-data.frame(dataSlice=dataSlice[i],
                        TLS=rep(colnames(TLS_data)[1:3],4),
                        localType=rep(c('Mal','Core','Boundary','Dispersion'),each=3),
                        R=c(Mal_cor,Core_cor,Bdy_cor,Dis_cor),
                        pValue=c(Mal_corP,Core_corP,Bdy_corP,Dis_corP))
  all_slice_cor<-rbind(all_slice_cor,slice_cor)
  print(dataSlice[i])
}
all_slice_cor$R[is.na(all_slice_cor$R)]<-1000
all_slice_cor$pValue[is.na(all_slice_cor$pValue)]<-1000
write.table(all_slice_cor,'E:/Mirror/ST_analysis/data/pic_data/TLS_stem_corr.txt',quote = F,sep = '\t',row.names = F)

all_slice_cor<-read.delim('E:/Mirror/ST_analysis/data/pic_data/TLS_stem_corr.txt',stringsAsFactors = F,check.names = F)

all_slice_cor$cancer<-unlist(lapply(strsplit(all_slice_cor$dataSlice,'/'),function(x)x[1]))
all_slice_cor$cancer<-substr(all_slice_cor$cancer,1,nchar(all_slice_cor$cancer)-2) %>% toupper()




####################################################################################
###所有切片放一起分三个得分画

plot_stem<-all_stem[all_stem$localType%in%'Core',]
plot_TLS<-all_TLS_data[all_TLS_data$FinalLocalType%in%'Core',]
rownames(plot_stem)<-paste0(plot_stem$dataset_slice,'_',plot_stem$cell_name)
rownames(plot_TLS)<-paste0(plot_TLS$dataset_slice,'_',plot_TLS$cell_name)
length(intersect(rownames(plot_stem),rownames(plot_TLS)))
plot_TLS<-plot_TLS[rownames(plot_stem),]
table(paste0(plot_TLS$dataset_slice,'_',plot_TLS$cell_name)==paste0(plot_stem$dataset_slice,'_',plot_stem$cell_name))
plot_TLS$CytoTRACE<-plot_stem$CytoTRACE
# rep(1:3,3)
# rep(1:3,each=3)

TLS_type<-c('TLS_2021','TLS_31942071','TLS_CancerSRT')

pdf(paste0(dir_picc,'panCancer_TLS_cor2.pdf'),width = 4,height = 3.5)
for(i in TLS_type){#i='TLS_CancerSRT'
  p_data<-plot_TLS[,c(i,'CytoTRACE')]
  colnames(p_data)<-c('TLS','stem')
  
  b <- ggplot(p_data, aes(x = stem, y = TLS))
  pp<-b + geom_point(color = "#f6c046",size=0.01,alpha=0.5)+
    geom_smooth(method = "lm", color = "black", fill = "gray60",size=0.7,se=T)+
    #geom_rug(aes(color =TLS_type)) +
    # scale_color_manual(values = c("#f6c046"))+
    # scale_fill_manual(values = c("#f6c046"))+
    stat_cor(method = "pearson", label.x = 0, label.y = 0.5)+
    theme_classic()+
    ggtitle(paste0('pan_cancer_',i))
  print(pp)
}
dev.off()







