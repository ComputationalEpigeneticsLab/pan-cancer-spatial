###计算SVG marker
library(ggplot2)
library(ggpubr)
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
#install.packages('Rfast2')
library(Rfast2)



dir_marker<-'E:/Mirror/ST_analysis/data/10X Visium/marker_SVG/'

dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
aa<-gsub('ST/expression_position/','',file_rds)
aa<-gsub('.rds','',aa)

dir_bdy<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"
file_bdy<-list.files(pattern = "_BdyTumorCore.txt",path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_bdy,"_"),function(x) x[1]))
table(aa==dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))


for(i in 1:229){
  #i=1
  if(!dir.exists(paste0(dir_marker,dataSet[i]))) dir.create(paste0(dir_marker,dataSet[i]))
  
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy$FinalLocalType[c(which(st_bdy$FinalLocalType=='Core'),
                          which(st_bdy$FinalLocalType=='Boundary'),
                          which(st_bdy$FinalLocalType=='Dispersion'))]<-'Mal'
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-st_bdy[colnames(st_rds),]
  table(st_bdy$FinalLocalType)
  
  ###marker计算
  Idents_new=st_bdy$FinalLocalType
  names(Idents_new)=colnames(st_rds)
  Seurat::Idents(object = st_rds)=Idents_new
  st_rds$seurat_clusters=Idents_new
  ST_marker<-FindAllMarkers(st_rds,only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25)
  length(unique(ST_marker$gene))
  table(ST_marker$cluster)
  write.table(ST_marker,paste0(dir_marker,dataset_slice[i],'_marker_Mal_Imm.txt'),quote = F,sep = '\t',row.names = F)
  
  ###SVG计算
  # st_rds<-FindSpatiallyVariableFeatures(st_rds, assay = "SCT", features=unique(ST_marker$gene), selection.method = "moransi")
  # aaa<-st_rds@assays[["SCT"]]@meta.features
  # aaa<-aaa[-which(aaa$MoransI_observed%in%NA),]
  # aaa$gene<-rownames(aaa)
  # aaa$MoransI_FDR<-p.adjust(aaa$MoransI_p.value,method = "BH")
  # aaa<-aaa[which(aaa$MoransI_observed>0&aaa$MoransI_FDR<0.05),]
  # aaa<-aaa[intersect(aaa$gene,unique(ST_marker$gene)),]
  # write.table(aaa,paste0(dir_marker,dataset_slice[i],'_SVG.txt'),quote = F,sep = '\t',row.names = F)
  print(dataset_slice[i])
}








