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



select_gene<-'SPP1'###TAGLN

for(i in 1:229){
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  plot_marker1<-data.frame(imagerow=st_rds@images[["slice1"]]@coordinates$imagerow,
                           imagecol=st_rds@images[["slice1"]]@coordinates$imagecol,
                           gene=st_rds@assays[["Spatial"]]@counts[select_gene,])
  
  
  p<-ggplot(plot_marker1, aes(x=imagerow, y=imagecol)) +
    geom_point(aes(colour=gene),size=.6) +
    scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(20),
                                      colorRampPalette(c("#F39C67","#B20A1C"))(90)) )+ #设置填充颜色
    labs(title = select_gene,
         x = "",
         y = "") +
    theme_bw()
  print(p)
  
}



dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD3/'
file_RCTD<-list.files(pattern = 'Deconvolution.txt',path = dir_RCTD,recursive = T)
dir_bdy<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"
file_bdy<-list.files(pattern = "_BdyTumorCore.txt",path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_bdy,"_"),function(x) x[1]))
table(dataset_slice==unlist(lapply(strsplit(file_RCTD,"_"),function(x) x[1])))

select_slice<-c('brca01/slice3','brca08/slice4','brca10/slice1','brca12/slice1','brca14/slice2','brca15/slice2',
                'gbm03/slice12','gbm04/slice7','hn-as02/slice4','lihc02/slice1','lihc02/slice12','lihc02/slice14',
                'lihc02/slice8','lihc02/slice9','lihc03/slice6','luad01/slice2','luad03/slice1','ovca01/slice1',
                'ovca03/slice1','ovca07/slice5','ovca07/slice8','pdac03/slice1','prad02/slice1','rcc01/slice4')
select_sit<-match(paste0(select_slice,'_Deconvolution.txt'),file_RCTD)
file_RCTD[select_sit]

for(i in select_sit){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[rownames(st_RCTD),]
  
  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c('Endothelial','Epithelial'))]
  
  # st_RCTD$Immune<-apply(st_RCTD[,setdiff(colnames(st_RCTD),c('CAF','TAM'))],1,sum)
  # st_RCTD<-st_RCTD[,c('CAF','TAM','Immune')]
  ####Imm合不合一起
  
  plot_data<-st_RCTD
  plot_data[st_bdy$FinalLocalType%in%c('Core',"Boundary",'Dispersion'),]<-0
  plot_data$Core<-0
  plot_data$Core[st_bdy$FinalLocalType%in%c('Core')]<-1
  plot_data$Boundary<-0
  plot_data$Boundary[st_bdy$FinalLocalType%in%c('Boundary')]<-1
  plot_data$Dispersion<-0
  plot_data$Dispersion[st_bdy$FinalLocalType%in%c('Dispersion')]<-1
  #plot_data$LocalType<-st_bdy$FinalLocalType
  plot_data$imagerow<-st_bdy$imagerow
  plot_data$imagecol<-st_bdy$imagecol
  plot_data<-plot_data[which(st_bdy$FinalLocalType!='not.defined'),]
  # table(apply(plot_data[1000:2000,1:12],1,sum)==1)
  # sum(plot_data[1,1:12])==1
  cellname <- colnames(plot_data)[1:(ncol(plot_data)-2)]
  
  p_pie<-ggplot() + scatterpie::geom_scatterpie(data = plot_data, 
                                                aes(x = imagerow, y = imagecol), col = cellname, color = NA,
                                                pie_scale = 0.35) + 
    coord_fixed(ratio = 1) + 
    scale_fill_manual(values = c('Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d','Immune'='#a9d38a',
                                 "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                                 "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                                 'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                                 "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',
                                 
                                 "GC B cells in the DZ"='#CC9933',
                                 "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',
                                 "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
                                 "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                                 "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066' )) + 
    theme_classic()+
    labs(x = "imagerow", y = "imagecol")+
    ggtitle(dataset_slice[i])
  pdf(paste0("E:/Mirror/ST_analysis/pic/re/Deconvolution/CAFTAM_decon/",gsub('/','_',dataset_slice[i]),"_scatterpie.pdf"),width = 5, height = 4)
  print(p_pie)
  dev.off()
  print(dataset_slice[i])
}




