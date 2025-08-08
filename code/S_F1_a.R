####癌症类型的切片数目和spot数码的环形堆砌图
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(circlize)


dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_bdy,'_Bdy'),function(x) x[1]))
cancer<-unlist(lapply(strsplit(dataSlice,'/'),function(x) x[1]))
cancer<-substr(cancer,1,(nchar(cancer)-2))%>%unique()

dir_TCGA<-'E:/Mirror/ST_analysis/data/TCGA/Core_MP/'
file_TCGA<-list.files(pattern = 'MP_1_',path = dir_TCGA)

dir_SC<-'E:/Mirror/ST_analysis/data/SC_data/'
file_SC<-list.files(pattern = 'SC_cellType.txt',path = dir_SC,recursive = T)

slice_num<-c()
spot_num<-c()
sc_cell_num<-c()
bulk_sample_num<-c()
for(i in 1:length(cancer)){#i=1
  SC_num<-read.delim(paste0(dir_SC,file_SC[i]),stringsAsFactors = F,check.names = F)
  sc_cell_num<-c(sc_cell_num,nrow(SC_num))
  
  bulk_num<-read.delim(paste0(dir_TCGA,file_TCGA[i]),stringsAsFactors = F,check.names = F)
  bulk_sample_num<-c(bulk_sample_num,nrow(bulk_num))
  
  cancer_file<-file_bdy[grep(cancer[i],file_bdy)]
  slice_num<-c(slice_num,length(cancer_file))
  
  cancer_spot_num<-c()
  for(j in 1:length(cancer_file)){#j=1
    st_bdy<-read.delim(paste0(dir_bdy,cancer_file[j]),stringsAsFactors = F,check.names = F)
    cancer_spot_num<-c(cancer_spot_num,nrow(st_bdy))
  }
  spot_num<-c(spot_num,sum(cancer_spot_num))
}

plot_bar<-data.frame(bulk_sample_num=bulk_sample_num,
                     sc_cell_num=sc_cell_num,
                     slice_num=slice_num,
                     spot_num=spot_num,
                     aaa=10,
                     row.names = toupper(cancer))
sum(plot_bar$spot_num)
sum(plot_bar$sc_cell_num[setdiff(1:19,c(4,7,9))])
sum(plot_bar$bulk_sample_num)
sum(plot_bar$slice_num)
write.table(plot_bar,'E:/Mirror/ST_analysis/data/pic_data/data_statistics.txt',quote = F,sep = '\t')

