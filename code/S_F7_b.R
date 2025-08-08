####映射8中CAF亚型到空转切片上
library(tidyverse)
library(ggplot2)
library(ggpubr)


dir_Cyto<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_Cyto<-list.files(pattern = 'new_assigned_locations.csv',path = dir_Cyto,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Cyto,'/'),function(x)x[1]))
###"brca09/slice1"  "gbm02/slice1"   "ipmn01/slice12"
lapply(c("brca09_slice1","gbm02_slice1","ipmn01_slice12"),function(x) grep(x,dataset_slice)) %>% unlist()
file_Cyto<-file_Cyto[-c(10,61,128)]
dataset_slice<-dataset_slice[-c(10,61,128)]

cancer<-unlist(lapply(strsplit(dataset_slice,'_'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()%>% toupper()


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)

for(i in 5:length(cancer)){#i=1
  cancer_file<-file_Cyto[grep(cancer[i],file_Cyto,ignore.case = T)]
  cancer_slice<-unlist(lapply(strsplit(cancer_file,'/'),function(x)x[1]))
  
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  meta_data<-sc_rds@meta.data
  rm(sc_rds)
  
  for(j in 1:length(cancer_file)){#j=1
    st_cyto<-read.csv(paste0(dir_Cyto,cancer_file[j]),stringsAsFactors = F,check.names = F)
    
    st_cyto$subCAF_celltype<-lapply(st_cyto$OriginalCID,function(x){#x=st_cyto$OriginalCID[1]
      cell_name<-unlist(strsplit(x,','))
      cell_celltype<-paste0(meta_data[cell_name,'subCAF'],collapse = ',')
    }) %>% unlist()
    
    write.table(st_cyto,paste0(dir_Cyto,cancer_slice[j],'/new_assigned_subCAF.txt'),quote = F,sep = '\t')
    print(cancer_slice[j])
  }
  Sys.sleep(5)
  
}



######绘制比例条形图
dir_subCAF<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_sub<-list.files(pattern = 'new_assigned_subCAF.txt',path=dir_subCAF,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_sub,'/'),function(x)x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-paste0(gsub('_','/',dataset_slice),'_BdyTumorCore.txt')

all_sliceCAF<-c()
for(i in 1:length(file_sub)){
  #i=1
  st_subCAF<-read.delim(paste0(dir_subCAF,file_sub[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[st_subCAF$SpotID,]
  st_bdy$subCAF_celltype<-st_subCAF$subCAF_celltype
  st_bdy<-st_bdy[st_bdy$FinalLocalType%in%c('Immune','Normal'),]
  
  cellType<-lapply(st_bdy$subCAF_celltype,function(x){
    strsplit(x,',') %>% unlist()
  }) %>% unlist()
  cellType<-cellType[grep('CAF',cellType)]
  table(cellType)
  
  cellType<-as.data.frame(table(cellType))
  cellType$ratio<-cellType$Freq/sum(cellType$Freq)
  cellType$slice<-dataset_slice[i]
  
  all_sliceCAF<-rbind(all_sliceCAF,cellType)
  print(dataset_slice[i])
}

all_sliceCAF$cancer<-unlist(lapply(strsplit(all_sliceCAF$slice,'_'),function(x)x[1]))
all_sliceCAF$cancer<-substr(all_sliceCAF$cancer,1,nchar(all_sliceCAF$cancer)-2) %>% toupper()
unique(all_sliceCAF$cellType)
all_sliceCAF$cellType<-factor(all_sliceCAF$cellType,levels = c('apCAF','iCAF',"ifnCAF","rCAF","mCAF","vCAF","dCAF","tCAF"))

p_compare<-ggplot(all_sliceCAF,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("apCAF"="#c26701",'iCAF'='#e2a503',"ifnCAF"="#ffc000","rCAF"="#bbc100",
                               "mCAF"="#007182","vCAF"="#00827c","dCAF"="#4a806c","tCAF"="#737f6b"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subCAF')
print(p_compare)

pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subCAF_bar.pdf',height = 4,width = 6.5)
print(p_compare)
dev.off()