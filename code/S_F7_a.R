####映射5中TAM亚型到空转切片上
library(tidyverse)
library(ggplot2)
library(ggpubr)


dir_Cyto<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_Cyto<-list.files(pattern = 'new_assigned_subCAF.txt',path = dir_Cyto,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Cyto,'/'),function(x)x[1]))
###"brca09/slice1"  "gbm02/slice1"   "ipmn01/slice12"
lapply(c("brca09_slice1","gbm02_slice1","ipmn01_slice12"),function(x) grep(x,dataset_slice)) %>% unlist()
# file_Cyto<-file_Cyto[-c(10,61,128)]
# dataset_slice<-dataset_slice[-c(10,61,128)]

cancer<-unlist(lapply(strsplit(dataset_slice,'_'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()%>% toupper()


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)
meta_file<-gsub('.rds','_meta.txt',file_sc)

for(i in 3:length(cancer)){#i=2
  cancer_file<-file_Cyto[grep(cancer[i],file_Cyto,ignore.case = T)]
  cancer_slice<-unlist(lapply(strsplit(cancer_file,'/'),function(x)x[1]))
  
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  meta_data<-sc_rds@meta.data
  rm(sc_rds)
  write.table(meta_data,paste0(dir_sc,meta_file[i]),quote = F,sep = '\t')
  
  for(j in 1:length(cancer_file)){#j=1
    st_cyto<-read.delim(paste0(dir_Cyto,cancer_file[j]),stringsAsFactors = F,check.names = F)
    st_cyto<-st_cyto[,-4]
    
    st_cyto$subCAFTAM<-lapply(st_cyto$OriginalCID,function(x){#x=st_cyto$OriginalCID[1]
      cell_name<-unlist(strsplit(x,','))
      cell_celltype<-paste0(meta_data[cell_name,'subCAFTAM'],collapse = ',')
    }) %>% unlist()
    
    write.table(st_cyto,paste0(dir_Cyto,cancer_slice[j],'/new_assigned_subCAFTAM.txt'),quote = F,sep = '\t')
    print(cancer_slice[j])
  }
  Sys.sleep(5)
  
}



######绘制比例条形图
dir_subTAM<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_sub<-list.files(pattern = 'new_assigned_subCAFTAM.txt',path=dir_subTAM,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_sub,'/'),function(x)x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-paste0(gsub('_','/',dataset_slice),'_BdyTumorCore.txt')

all_sliceTAM<-c()
for(i in 1:length(file_sub)){
  #i=1
  st_subTAM<-read.delim(paste0(dir_subTAM,file_sub[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[st_subTAM$SpotID,]
  st_bdy$subCAFTAM<-st_subTAM$subCAFTAM
  st_bdy<-st_bdy[st_bdy$FinalLocalType%in%c('Immune','Normal'),]
  
  cellType<-lapply(st_bdy$subCAFTAM,function(x){
    strsplit(x,',') %>% unlist()
  }) %>% unlist()
  cellType<-cellType[grep('TAM',cellType)]
  cellType[grep('FCN1_TAM',cellType)]<-'SPP1_TAM'
  table(cellType)
  
  cellType<-as.data.frame(table(cellType))
  cellType$ratio<-cellType$Freq/sum(cellType$Freq)
  cellType$slice<-dataset_slice[i]
  
  all_sliceTAM<-rbind(all_sliceTAM,cellType)
  print(dataset_slice[i])
}

all_sliceTAM$cancer<-unlist(lapply(strsplit(all_sliceTAM$slice,'_'),function(x)x[1]))
all_sliceTAM$cancer<-substr(all_sliceTAM$cancer,1,nchar(all_sliceTAM$cancer)-2) %>% toupper()
unique(all_sliceTAM$cellType)
write.table(all_sliceTAM,'E:/Mirror/ST_analysis/data/pic_data/subTAM_all.txt',quote = F,sep = '\t',row.names = F)
all_sliceTAM<-read.delim('E:/Mirror/ST_analysis/data/pic_data/subTAM_all.txt',stringsAsFactors = F,check.names = F)

all_sliceTAM$cellType<-factor(all_sliceTAM$cellType,levels = c("SPP1_TAM",'C1Q_TAM',"Microglial_TAM","blood_TAM"))
ppp_data<-aggregate(all_sliceTAM$Freq,by=list(all_sliceTAM$cellType,all_sliceTAM$cancer),sum)
site_cancer<-lapply(unique(ppp_data$Group.2), function(x){
  aaa<-ppp_data[ppp_data$Group.2%in%x,]
  return(aaa[which(aaa[,1]=='SPP1_TAM'),3]-aaa[which(aaa[,1]=='C1Q_TAM'),3])
}) %>% unlist()
names(site_cancer)<-unique(ppp_data$Group.2)
site_cancer<-site_cancer[order(abs(site_cancer),decreasing = F)]

all_sliceTAM$cancer<-factor(all_sliceTAM$cancer,levels = names(site_cancer))

p_compare<-ggplot(all_sliceTAM,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("C1Q_TAM"="#71a3c5",'FCN1_TAM'='#cc889f',"SPP1_TAM"="#b83231",
                               "Microglial_TAM"="#80c598","blood_TAM"="#1d804e"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subTAM')
print(p_compare)

pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subTAM_bar_SPP1&FN1.pdf',height = 4,width = 6.5)
print(p_compare)
dev.off()

