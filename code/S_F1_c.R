# devtools::install_github("data2intelligence/SpaCET")
# devtools::install_github("JEFworks/MUDAN")
# remotes::install_local("E:/Mirror/ST_analysis/program/SpaCET/MUDAN-master.zip", upgrade=F, dependencies=T)
# remotes::install_local("E:/Mirror/ST_analysis/program/SpaCET/SpaCET-main.zip", upgrade=F, dependencies=T)
#remotes::install_local("/data/zhouweiwei/1-project/spatial_analysis/code/SpaCET-main.zip", upgrade=F, dependencies=T)
library(SpaCET)
library(Seurat)
library(tidyverse)

# install.packages("/data/zhouweiwei/spatial/code/SpaCET/sva_3.52.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/spatial/code/SpaCET/genefilter_1.86.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/spatial/code/SpaCET/BiRewire_3.36.0.tar.gz", repos = NULL, type = "source")
# remotes::install_local("/data/zhouweiwei/spatial/code/SpaCET/MUDAN-master.zip", upgrade=F, dependencies=T)
# remotes::install_local("/data/zhouweiwei/spatial/code/SpaCET/SpaCET-main.zip", upgrade=F, dependencies=T)

# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/matrixStats_1.4.1.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/MatrixGenerics_1.18.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/statmod_1.5.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/limma_3.62.1.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/genefilter_1.86.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/edgeR_4.4.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/sva_3.52.0.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/slam_0.1-55.tar.gz", repos = NULL, type = "source")
# install.packages("/data/zhouweiwei/1-project/spatial_analysis/code/BiRewire_3.36.0.tar.gz", repos = NULL, type = "source")
# remotes::install_local("/data/zhouweiwei/1-project/spatial_analysis/code/MUDAN-master.zip", upgrade=F, dependencies=T)
# remotes::install_local("/data/zhouweiwei/1-project/spatial_analysis/code/SpaCET-main.zip", upgrade=F, dependencies=T)


library(SpaCET)
library(Seurat)
library(tidyverse)

dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
data_slice<-gsub('ST/expression_position/','',file_rds)
data_slice<-gsub('.rds','',data_slice)
dir_out<-'/data/zhouweiwei/ST_analysis/data/10X Visium/SpaCET/'
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x)x[1]))
cancer<-substr(dataSet,1,nchar(dataSet)-2) %>% toupper()
unique(cancer)
cancer[grep('CRC',cancer)]<-'COAD'
cancer[grep('GIST',cancer)]<-'STAD'
cancer[grep('HGSC',cancer)]<-'OVCA'
cancer[grep('IPMN',cancer)]<-'PAAD'
cancer[grep('MIBC',cancer)]<-'BLCA'
cancer[grep('OSCC',cancer)]<-'HNSC'
cancer[grep('PDAC',cancer)]<-'PAAD'
cancer[grep('RCC',cancer)]<-'KIRC'
cancer[grep('CSCC',cancer)]<-'HNSC'

load( system.file("extdata", 'cancerDictionary.rda', package = 'SpaCET') )
cancerTypes <- unique(c(names(cancerDictionary$CNA),names(cancerDictionary$expr)))
cancerTypes <- sapply(strsplit(cancerTypes,"_",fixed=T),function(x) return(x[2]))

file_already<-list.files(pattern = 'SpaCET.txt',path = dir_out,recursive = T)
dataSlice_already<-unlist(lapply(strsplit(file_already,'_'),function(x)x[1]))

for(i in 1:229){
  #i=4
  if(!(data_slice[i]%in%dataSlice_already)){
    if(!dir.exists(paste0(dir_out,dataSet[i]))) dir.create(paste0(dir_out,dataSet[i]))
    
    if(cancer[i]%in%cancerTypes){
      cancer_lab<-cancer[i]
    }else{
      cancer_lab<-'PANCAN'
    }
    
    Seurat_obj<-readRDS(paste0(dir_rds,file_rds[i]))
    
    if(!'row'%in%colnames(Seurat_obj@images[["slice1"]]@coordinates)){
      Seurat_obj@images[["slice1"]]@coordinates$row<-Seurat_obj@images[["slice1"]]@coordinates$imagerow
      Seurat_obj@images[["slice1"]]@coordinates$col<-Seurat_obj@images[["slice1"]]@coordinates$imagecol
    }
    
    SpaCET_obj<-convert.Seurat(Seurat_obj)
    
    SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType=cancer_lab, coreNo=1)
    aa<-as.data.frame(t(SpaCET_obj@results$deconvolution$propMat))
    rownames(aa)<-SpaCET_obj@input[["spotCoordinates"]][["barcode"]]
    write.table(aa,paste0(dir_out,data_slice[i],'_SpaCET.txt'),quote = F,sep = '\t')
    print(data_slice[i])
  }
  
}



############在85上跑
library(SpaCET)
library(Seurat)
library(tidyverse)

dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
data_slice<-gsub('ST/expression_position/','',file_rds)
data_slice<-gsub('.rds','',data_slice)
dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/SpaCET/'
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x)x[1]))
cancer<-substr(dataSet,1,nchar(dataSet)-2) %>% toupper()
unique(cancer)
cancer[grep('CRC',cancer)]<-'COAD'
cancer[grep('GIST',cancer)]<-'STAD'
cancer[grep('HGSC',cancer)]<-'OVCA'
cancer[grep('IPMN',cancer)]<-'PAAD'
cancer[grep('MIBC',cancer)]<-'BLCA'
cancer[grep('OSCC',cancer)]<-'HNSC'
cancer[grep('PDAC',cancer)]<-'PAAD'
cancer[grep('RCC',cancer)]<-'KIRC'
cancer[grep('CSCC',cancer)]<-'HNSC'

load( system.file("extdata", 'cancerDictionary.rda', package = 'SpaCET') )
cancerTypes <- unique(c(names(cancerDictionary$CNA),names(cancerDictionary$expr)))
cancerTypes <- sapply(strsplit(cancerTypes,"_",fixed=T),function(x) return(x[2]))

file_already<-list.files(pattern = 'SpaCET.txt',path = dir_out,recursive = T)
dataSlice_already<-unlist(lapply(strsplit(file_already,'_'),function(x)x[1]))

for(i in 1:100){
  #i=4
  if(!(data_slice[i]%in%dataSlice_already)){
    if(!dir.exists(paste0(dir_out,dataSet[i]))) dir.create(paste0(dir_out,dataSet[i]))
    
    if(cancer[i]%in%cancerTypes){
      cancer_lab<-cancer[i]
    }else{
      cancer_lab<-'PANCAN'
    }
    
    Seurat_obj<-readRDS(paste0(dir_rds,file_rds[i]))
    
    if(!'row'%in%colnames(Seurat_obj@images[["slice1"]]@coordinates)){
      Seurat_obj@images[["slice1"]]@coordinates$row<-Seurat_obj@images[["slice1"]]@coordinates$imagerow
      Seurat_obj@images[["slice1"]]@coordinates$col<-Seurat_obj@images[["slice1"]]@coordinates$imagecol
    }
    
    SpaCET_obj<-convert.Seurat(Seurat_obj)
    
    SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType=cancer_lab, coreNo=1)
    aa<-as.data.frame(t(SpaCET_obj@results$deconvolution$propMat))
    rownames(aa)<-SpaCET_obj@input[["spotCoordinates"]][["barcode"]]
    write.table(aa,paste0(dir_out,data_slice[i],'_SpaCET.txt'),quote = F,sep = '\t')
    print(data_slice[i])
  }
  
}


####基于SpaCET的结果评判copykat注释的准确性
library(SpaCET)
library(Seurat)
library(pROC)
library(ggplot2)
library(ggpubr)
library(glmnet)
library(foreign)
library(caret)
library(rlang)
library(xgboost)
library(Matrix)


dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dir_SpaCET<-'E:/Mirror/ST_analysis/data/10X Visium/SpaCET/'
file_CET<-list.files(pattern = 'SpaCET.txt',path = dir_SpaCET,recursive = T)

dataSlice<-unlist(lapply(strsplit(file_CET,'_'),function(x)x[1]))
table(dataSlice==unlist(lapply(strsplit(file_bdy,'_'),function(x)x[1])))



pdf('E:/Mirror/ST_analysis/pic/AUC/SpaCET_AUC_pin2.pdf',width = 4,height = 4)

i=1
st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
st_CET<-read.delim(paste0(dir_SpaCET,file_CET[i]),stringsAsFactors = F,check.names = F)
st_bdy<-st_bdy[rownames(st_CET),]
st_bdy<-st_bdy[which(st_bdy$copykat.pred!='not.defined'),]
st_CET<-st_CET[rownames(st_bdy),]

AUC_data<-data.frame(Mal_score=st_CET$Malignant,
                     copykat_lab=st_bdy$copykat.pred)
AUC_data$copykat_lab<-factor(AUC_data$copykat_lab,levels = c('diploid','aneuploid'))

plot(roc(AUC_data$copykat_lab,AUC_data$Mal_score),col='#B2533E')

for(i in 2:length(file_CET)){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_CET<-read.delim(paste0(dir_SpaCET,file_CET[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[rownames(st_CET),]
  st_bdy<-st_bdy[which(st_bdy$copykat.pred!='not.defined'),]
  st_CET<-st_CET[rownames(st_bdy),]
  
  AUC_data<-data.frame(Mal_score=st_CET$Malignant,
                       copykat_lab=st_bdy$copykat.pred)
  AUC_data$copykat_lab<-factor(AUC_data$copykat_lab,levels = c('diploid','aneuploid'))
  AUC_score<-as.numeric(roc(AUC_data$copykat_lab,AUC_data$Mal_score)$auc)
  AUC_score<-c(0.7,0.8)
  
  plot(roc(AUC_data$copykat_lab,AUC_data$Mal_score),add=T,col='#B2533E')
  
  # legend_lab<-paste0(dataSlice[1:2],' (',round(AUC_score,3),')')
  # legend("bottomright",
  #        legend=legend_lab,
  #        col=c('blue4','dodgerblue'),
  #        lwd=2,bty='n',inset=c(0.28,0))
}
dev.off()



####AUC大于0.75的算个比例
plot_AUC<-all_AUC[which(all_AUC$AUC>0.7),]

####环形条状图
plot_data<-data.frame(AUC=c('>0.7','<=0.7'),
                      num=c(211,18))###c(205,24)
plot_data$value<-plot_data$num/sum(plot_data$num)
plot_data$AUC<-factor(plot_data$AUC,levels = c('>0.7','<=0.7'))

p5 <- ggplot(plot_data,aes(x = 1, y = value, fill = AUC)) +
  geom_col(colour = "white")+ 
  coord_polar(theta = "y", start = 1.65) +
  geom_text(aes(label = paste0(round(value * 100, 2), "%")),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c(">0.7"="#f4b184","<=0.7"="#8fabdd"))+
  xlim(c(-0.2, 2)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
print(p5)
pdf(paste0("E:/Mirror/ST_analysis/pic/AUC/",'AUC_pie_0.7.pdf'),width = 7,height = 7)
print(p5)
dev.off()






