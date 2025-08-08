#####克隆型
####所有恶性区域的克隆型个数
####克隆型小于等于3和大于3的切片 比较其干性均值高低
####分core bdy dis区域
library(coda)
library(Seurat)
library(infercnv)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(magic)
library(ggpubr)


###参考基因组 
geneFile <- read.table("E:/Mirror/ST_analysis/data/gencode.v43.annotation.txt",header = T,sep = "\t",stringsAsFactors = F)
geneFile<-geneFile[!duplicated(geneFile$gene_name),]
rownames(geneFile)=geneFile$gene_name

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
#file_bdy<-file_bdy[-grep('ovca02',file_bdy)]  ####还在跑
#file_bdy<-file_bdy[-grep('luad01/slice2',file_bdy)] ####inferCNV未知原因报错
data_slice<-unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x) x[1]))
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x) x[1]))

dir_CNV<-'E:/Mirror/ST_analysis/data/10X Visium/inferCNV_R2/'
file_clone<-list.files(pattern = 'hmm_mode-subclusters.cell_groupings',path = dir_CNV,recursive = T)
table(unlist(lapply(strsplit(file_clone,'/'),function(x) x[1]))==gsub('/','_',data_slice))
###统计克隆型需要用到cell_groupings文件
###绘制染色体克隆扩增热图需要用run.final.infercnv_obj文件

###统计克隆型
core_clone_num<-c()
core_num<-c()
bdy_clone_num<-c()
bdy_num<-c()
dis_clone_num<-c()
dis_num<-c()
for(i in 1:length(file_bdy)){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  table(st_bdy$FinalLocalType)
  
  clone_type<-read.delim(paste0(dir_CNV,file_clone[i]),stringsAsFactors = F,check.names = F)
  rownames(clone_type)<-clone_type$cell
  st_bdy<-st_bdy[rownames(clone_type),]
  clone_type$FinalLocalType<-st_bdy$FinalLocalType
  clone_type<-clone_type[st_bdy$cell_name[st_bdy$FinalLocalType%in%c('Core','Boundary','Dispersion')],]
  clone_type$cell_group_name<-unlist(lapply(strsplit(clone_type$cell_group_name,'.all_observations.'),function(x) x[2]))
  table(clone_type$cell_group_name)
  
  core_clone_num<-c(core_clone_num,unique(clone_type$cell_group_name[which(clone_type$FinalLocalType=='Core')]) %>% length())
  core_num<-c(core_num,length(which(clone_type$FinalLocalType=="Core")))
  bdy_clone_num<-c(bdy_clone_num,unique(clone_type$cell_group_name[which(clone_type$FinalLocalType=='Boundary')]) %>% length())
  bdy_num<-c(bdy_num,length(which(clone_type$FinalLocalType=="Boundary")))
  dis_clone_num<-c(dis_clone_num,unique(clone_type$cell_group_name[which(clone_type$FinalLocalType=='Dispersion')]) %>% length())
  dis_num<-c(dis_num,length(which(clone_type$FinalLocalType=="Dispersion")))
  
}

clone_data<-data.frame(data_slice=data_slice,
                       core_num=core_num,
                       core_clone_num=core_clone_num,
                       bdy_num=bdy_num,
                       bdy_clone_num=bdy_clone_num,
                       dis_num=dis_num,
                       dis_clone_num=dis_clone_num)
clone_data<-clone_data[which(clone_data$core_num!=0&clone_data$bdy_num!=0&clone_data$dis_num!=0),]

###每个切片的干性
all_slice_stem<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_stem.txt',
                           stringsAsFactors = F,check.names = F)
all_slice_stem_core<-all_slice_stem[which(all_slice_stem$localType=='Core'),]
stem_mean<-aggregate(all_slice_stem_core$CytoTRACE,by=list(all_slice_stem_core$dataset_slice),mean)
clone_data$core_stem<-stem_mean$x[match(clone_data$data_slice,stem_mean$Group.1)]
all_slice_stem_bdy<-all_slice_stem[which(all_slice_stem$localType=='Boundary'),]
stem_mean<-aggregate(all_slice_stem_bdy$CytoTRACE,by=list(all_slice_stem_bdy$dataset_slice),mean)
clone_data$bdy_stem<-stem_mean$x[match(clone_data$data_slice,stem_mean$Group.1)]
all_slice_stem_dis<-all_slice_stem[which(all_slice_stem$localType=='Dispersion'),]
stem_mean<-aggregate(all_slice_stem_dis$CytoTRACE,by=list(all_slice_stem_dis$dataset_slice),mean)
clone_data$dis_stem<-stem_mean$x[match(clone_data$data_slice,stem_mean$Group.1)]


clone_data$core_group<-'high'
clone_data$core_group[which(clone_data$core_clone_num<=3)]<-'low'
clone_data$bdy_group<-'high'
clone_data$bdy_group[which(clone_data$bdy_clone_num<=3)]<-'low'
clone_data$dis_group<-'high'
clone_data$dis_group[which(clone_data$dis_clone_num<=3)]<-'low'
clone_data$cancer<-unlist(lapply(strsplit(clone_data$data_slice,'/'),function(x)x[1]))
clone_data$cancer<-substr(clone_data$cancer,1,nchar(clone_data$cancer)-2) %>% toupper()

pdf('E:/Mirror/ST_analysis/pic/clone_num/cloneHiLoStem_CoreBdyDis2.pdf',width = 4,height = 4)
for(i in c('core','bdy','dis')){#i='core'
  plot_data<-clone_data[,grep(i,colnames(clone_data))]
  colnames(plot_data)[3:4]<-c('stemness','clone_group')
  
  p_box<-ggplot(plot_data, aes(x = clone_group, y = stemness,fill=clone_group))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = clone_group), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("high"="#B79D71","low"="#2D59A6"))+
    scale_color_manual(values = c("high"="#B79D71","low"="#2D59A6"))+
    ggtitle(i)+
    theme(plot.title = element_text(hjust = 0.4))+
    theme(plot.title = element_text(size = 15))+
    ylab('stemness')+
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 1.5,label.y =1)+
    #stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
    #theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  print(p_box)
  
}
dev.off()





