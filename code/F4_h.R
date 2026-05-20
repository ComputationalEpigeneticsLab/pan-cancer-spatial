###对每个切片的T细胞主导的spot计算8种T细胞状态
library(TCellSI)
library(dplyr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(monocle)
library(igraph)
library(RobustRankAggreg)
library(jsonlite)

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_TcellSI/'

dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_copykat/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/RCTD3/'

file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = '_Deconvolution.txt',path = dir_RCTD,recursive = T)
file_bdy<-file_bdy[-grep('slice',file_bdy)]
file_rds<-file_rds[-grep('slice',file_rds)]
file_RCTD<-file_RCTD[-grep('slice',file_RCTD)]
length(file_bdy)
file_bdy[1:10]
length(file_rds)
file_rds[1:10]
length(file_RCTD)
file_RCTD[1:10]

dataslice<-unlist(lapply(strsplit(file_bdy,'_BdyCoreBud'),function(x)x[1]))
table(dataslice==gsub('processed_|.rds','',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_Deconvolution'),function(x)x[1])))
cancer<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))
slice<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[3]))
cancer[1:10]
dataset[1:10]
slice[1:10]


for(i in 1:length(file_rds)){
  #i=1
  if(!dir.exists(paste0(dir_out,cancer[i]))) dir.create(paste0(dir_out,cancer[i]))
  if(!dir.exists(paste0(dir_out,cancer[i],'/',dataset[i]))) dir.create(paste0(dir_out,cancer[i],'/',dataset[i]))
  
  
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
  st_bdy<-st_bdy[rownames(st_RCTD),]
  st_RCTD$celltype<-apply(st_RCTD,1,function(x){
    names(x)<-colnames(st_RCTD)
    x<-x[order(x,decreasing = T)]
    return(names(x)[1])
  })
  st_RCTD$FinalLocalType<-st_bdy$FinalLocalType
  st_RCTD<-st_RCTD[which(st_RCTD$celltype=='T lymphocytes'),]
  T_spot<-rownames(st_RCTD)[which(st_RCTD$FinalLocalType=='Immune'|st_RCTD$FinalLocalType=='Normal')]
  
  tryCatch({
    if(length(T_spot)>=5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_count<-st_rds@assays[["Spatial"]]@counts[,T_spot] %>% as.matrix()
      st_scores <- TCSS_Calculate(st_count) #%>% t()
      st_scores<-t(st_scores)
      write.table(st_scores,paste0(dir_out,cancer[i],'/',dataset[i],'/',slice[i],'_TCellSI.txt'),
                  quote = F,sep = '\t')
    }
    print(dataslice[i])
  },error = function(e){
    stopMessage_aaa<-"Unable to convert data"
    print(paste0(dataslice[i],'_stop'))
  })
  
}
# paste0(dir_out,cancer[i],'/',dataset[i],'/',slice[i],'_startRes_CAF.txt')


dir_near<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_hop10/'
file_near<-list.files(pattern = '_nearSpotStep1to10.txt',path = dir_near,recursive = T)
dir_Tcell<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_TcellSI/'
file_Tcell<-list.files(pattern = 'TCellSI.txt',path = dir_Tcell,recursive = T)
length(file_Tcell)
file_Tcell[1:10]

dataSlice<-unlist(lapply(strsplit(file_Tcell,'_T'),function(x)x[1]))
file_near<-file_near[match(paste0(dataSlice,'_nearSpotStep1to10.txt'),file_near)]
file_near[1:10]
dataSlice[1:10]

all_step_Tcell<-c()
for(i in 1:length(file_Tcell)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  st_Tcell[1:3,]
  colnames(st_near)
  
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  length(intersect(spot_step1,rownames(st_Tcell)))
  spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
  spot_step2<-union(spot_step1,spot_step2)
  spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
  spot_step3<-union(spot_step2,spot_step3)
  spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
  spot_step4<-union(spot_step3,spot_step4)
  spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
  spot_step5<-union(spot_step4,spot_step5)
  
  
  step_Tcell<-lapply(list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5),function(x){
    spot_Tcell<-intersect(x,rownames(st_Tcell))
    near_Tcell<-rep(0,8)
    names(near_Tcell)<-colnames(st_Tcell)
    if(length(spot_Tcell)>1){
      near_Tcell<-apply(st_Tcell[spot_Tcell,],2,mean)
    }
    return(near_Tcell)
  })
  step_Tcell<-do.call(rbind,step_Tcell)
  rownames(step_Tcell)<-paste0('step_',1:5)
  step_Tcell<-reshape2::melt(as.matrix(step_Tcell))
  step_Tcell$slice<-dataSlice[i]
  
  all_step_Tcell<-rbind(all_step_Tcell,step_Tcell)
}

all_step_Tcell$cancer<-unlist(lapply(strsplit(all_step_Tcell$slice,'/'),function(x)x[1]))
# all_step_Tcell$cancer<-substr(all_step_Tcell$cancer,1,nchar(all_step_Tcell$cancer)-2) %>% toupper()
all_step_Tcell<-all_step_Tcell[which(all_step_Tcell$value!=0),]



write.table(all_step_Tcell,paste0(dir_Tcell,'all_step_TcellSI.txt'),quote = F,sep = '\t',row.names = F)
all_step_Tcell<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/new_st_TcellSI/all_step_TcellSI.txt',stringsAsFactors = F,check.names = F)
cancer<-unique(all_step_Tcell$cancer)


cancer_col<-read.delim('E:/Mirror/ST_vascul/35cancer_color.txt',stringsAsFactors = F,check.names = F)
cancer_col2<-cancer_col$color
names(cancer_col2)<-cancer_col$cancer
####按每个step画
#pdf('E:/Mirror/ST_analysis/pic/lihc_CAF_isolation/TcellSI/step_state_cluster.pdf',width = 4,height = 8)
pdf('E:/Mirror/ST_analysis/data/10X Visium/new_st_TcellSI/pic/step_state.pdf',width = 8,height = 3)
for(i in 1:length(unique(all_step_Tcell$Var1))){
  #i=1
  stat_data<-all_step_Tcell[all_step_Tcell$Var1%in%unique(all_step_Tcell$Var1)[i],]
  p_data<-as.data.frame(reshape2::acast(stat_data[,c('Var2','slice','value')],slice~Var2))
  p_data$cancer<-unlist(lapply(strsplit(rownames(p_data),'/'),function(x)x[1]))
  # p_data$cancer<-substr(p_data$cancer,1,nchar(p_data$cancer)-2) %>% toupper()
  unique(p_data$cancer)
  
  
  # la <- rowAnnotation(df = data.frame(cancer=p_data$cancer),
  #                     col = list(cancer=col_data[unique(p_data$cancer)])
  # )
  
  
  heat_data<-t(as.matrix(p_data[,c(1,2,4,6,5,7,3,8)]))
  
  h_order<-apply(heat_data,2,function(x) sum(x,na.rm = T))
  h_order<-h_order[order(h_order,decreasing = T)]
  heat_data<-heat_data[,names(h_order)]
  
  p_data<-p_data[names(h_order),]
  ha <- HeatmapAnnotation(df = data.frame(cancer=p_data$cancer),
                          col = list(cancer=cancer_col2[unique(p_data$cancer)])
  )
  # range(heat_data,na.rm = T)
  # seq(1,5,length.out=5)
  # ?seq()
  up_thre<-max(range(heat_data,na.rm = T))
  col_fun <- colorRamp2(
    seq(0, ifelse(up_thre>0.2,0.2,up_thre),length.out=10), 
    rev(c('#ab2c74','#B84687','#C5619B','#D98AB8','#e7a5cc','#ECCBD6','#EEDADA','#f0eadf','#bbdf8e','#82AF5D')) 
    #c("#16160e","#707020",'#b0b018','#f2e500')
  )
  p1<-Heatmap(heat_data,
              col = col_fun,
              top_annotation =ha,#####顶部注释
              #left_annotation=la,#####左侧注释
              show_row_names = T,#####不显示行名
              show_column_names =F,#####不显示列名
              cluster_rows =F,
              show_row_dend = F,
              show_column_dend = F,
              cluster_columns = F,
              name="score",
              column_title=unique(all_step_Tcell$Var1)[i],
              row_names_gp = gpar(fontsize = 10))
  print(p1)
}
dev.off()








