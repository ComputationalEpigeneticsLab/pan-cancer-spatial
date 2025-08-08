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

######bdy邻居的TAM变化
dir_subTAM<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_sub<-list.files(pattern = 'new_assigned_subCAFTAM.txt',path=dir_subTAM,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_sub,'/'),function(x)x[1]))

dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(gsub('_','/',dataset_slice),'_BdyNearCore.txt')

all_step_TAM<-c()
for(i in 1:length(file_sub)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_subTAM<-read.delim(paste0(dir_subTAM,file_sub[i]),stringsAsFactors = F,check.names = F)
  rownames(st_subTAM)<-st_subTAM$SpotID
  
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  spot_step1<-intersect(spot_step1,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
  #spot_step2<-union(spot_step1,spot_step2)
  spot_step2<-intersect(spot_step2,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
  #spot_step3<-union(spot_step2,spot_step3)
  spot_step3<-intersect(spot_step3,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
  #spot_step4<-union(spot_step3,spot_step4)
  spot_step4<-intersect(spot_step4,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
  #spot_step5<-union(spot_step4,spot_step5)
  spot_step5<-intersect(spot_step5,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_list<-list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5)
  names(spot_list)<-c('step1','step2','step3','step4','step5')
  step_TAM<-lapply(1:5,function(x){#x=1
    spot_TAM<-st_subTAM[intersect(spot_list[[x]],rownames(st_subTAM)),]
    cellType<-lapply(spot_TAM$subCAFTAM,function(x){
      strsplit(x,',') %>% unlist()
    }) %>% unlist()
    cellType<-cellType[grep('TAM',cellType)]
    cellType[grep('FCN1_TAM',cellType)]<-'SPP1_TAM'
    if(length(cellType)>5){
      table(cellType)
      cellType<-as.data.frame(table(cellType))
      cellType$ratio<-cellType$Freq/sum(cellType$Freq)
      cellType$step<-names(spot_list)[x]
    }else{
      cellType<-data.frame(cellType='',Freq=0,ratio=0,step='no_step')
    }
    return(cellType)
  })
  step_TAM<-do.call(rbind,step_TAM)
  step_TAM$slice<-dataset_slice[i]
  
  all_step_TAM<-rbind(all_step_TAM,step_TAM)
  print(dataset_slice[i])
}

all_step_TAM$cancer<-unlist(lapply(strsplit(all_step_TAM$slice,'_'),function(x)x[1]))
all_step_TAM$cancer<-substr(all_step_TAM$cancer,1,nchar(all_step_TAM$cancer)-2) %>% toupper()
all_step_TAM<-all_step_TAM[which(all_step_TAM$Freq!=0),]
cancer<-unique(all_step_TAM$cancer)
write.table(all_step_TAM,'E:/Mirror/ST_analysis/data/pic_data/all_step_subTAM.txt',quote = F,sep = '\t',row.names = F)

####绘制成热图 每个step画一个图
###一行一个切片 一列是一类TAM
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

all_step_TAM$cellType<-as.vector(all_step_TAM$cellType)
unique(all_step_TAM$cellType)

col_data<-c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
            'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
            'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
            'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
range(all_step_TAM$ratio)


###################################################
####按每个step画
pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subTAM_step_SPP1&FN1.pdf',width = 4,height = 8)
for(i in 1:length(unique(all_step_TAM$step))){
  #i=1
  stat_data<-all_step_TAM[all_step_TAM$step%in%unique(all_step_TAM$step)[i],]
  p_data<-as.data.frame(reshape2::acast(stat_data[,c('cellType','slice','ratio')],slice~cellType))
  colnames(p_data)
  p_data<-p_data[,c('SPP1_TAM','C1Q_TAM','Microglial_TAM','blood_TAM')]
  p_data$cancer<-unlist(lapply(strsplit(rownames(p_data),'_'),function(x)x[1]))
  p_data$cancer<-substr(p_data$cancer,1,nchar(p_data$cancer)-2) %>% toupper()
  unique(p_data$cancer)
  
  
  la <- rowAnnotation(df = data.frame(cancer=p_data$cancer),
                      col = list(cancer=col_data[unique(p_data$cancer)])
  )
  heat_data<-as.matrix(p_data[,1:4])
  # range(heat_data,na.rm = T)
  # seq(1,5,length.out=5)
  # ?seq()
  #up_thre<-max(range(heat_data,na.rm = T))
  col_fun <- colorRamp2(
    seq(0, 0.4,length.out=10), 
    rev(c('#ab2c74','#B84687','#C5619B','#D98AB8','#e7a5cc','#ECCBD6','#EEDADA','#f0eadf','#bbdf8e','#82AF5D')) 
    #c("#16160e","#707020",'#b0b018','#f2e500')
  )
  p1<-Heatmap(heat_data,
              col = col_fun,
              #top_annotation =ha,#####顶部注释
              left_annotation=la,#####左侧注释
              show_row_names = F,#####不显示行名
              show_column_names =T,#####不显示列名
              cluster_rows =F,
              cluster_columns = F,
              name="ratio",
              column_title=unique(all_step_TAM$step)[i],
              row_names_gp = gpar(fontsize = 10))
  print(p1)
}
dev.off()
