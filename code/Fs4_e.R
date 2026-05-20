####TF RBP与干性的关系
###将干性相关性高的排上面
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)


all_RBP<-read.delim('/pic_data/TF_RBP/all_RBP_quantile_mean.txt',stringsAsFactors = F,check.names = F)
all_TF<-read.delim('/pic_data/TF_RBP/all_TF_quantile_mean.txt',stringsAsFactors = F,check.names = F)

all_RBP_TF<-rbind(all_TF,all_RBP)
all_RBP_TF<-all_RBP_TF[,c(grep('core',colnames(all_RBP_TF)))]

###每个切片的干性
all_slice_stem<-read.delim('/pic_data/all_slice_stem.txt',
                           stringsAsFactors = F,check.names = F)
all_slice_stem_core<-all_slice_stem[which(all_slice_stem$localType=='Core'),]
stem_mean<-aggregate(all_slice_stem_core$CytoTRACE,by=list(all_slice_stem_core$dataset_slice),mean)
rownames(stem_mean)<-paste0(stem_mean$Group.1,'_core')
all_RBP_TF_cor<-all_RBP_TF[,paste0(stem_mean$Group.1,'_core')]
cor_order<-apply(all_RBP_TF_cor,1,function(x){#x=as.vector(as.matrix(all_RBP_TF_cor[1,1:227]))
  cor.test(x,stem_mean$x)[["estimate"]][["cor"]]
})
cor_order<-cor_order[order(cor_order,decreasing = T)]

all_RBP_TF<-all_RBP_TF[c(names(cor_order)[names(cor_order)%in%TF],names(cor_order)[names(cor_order)%in%RBP]),]
all_RBP_TF<-all_RBP_TF[c(c('MYC','SOX9'),setdiff(rownames(all_RBP_TF),c('MYC','SOX9'))),]



meta_data<-data.frame(dataSlice=colnames(all_RBP_TF),
                      order=apply(all_RBP_TF,2,function(x){sum(x,na.rm = T)}))
meta_data$region<-unlist(lapply(strsplit(meta_data$dataSlice,'_'),function(x)x[2]))
meta_data$cancer<-unlist(lapply(strsplit(meta_data$dataSlice,'/'),function(x)x[1]))
meta_data$cancer<-substr(meta_data$cancer,1,nchar(meta_data$cancer)-2) %>% toupper()
meta_data<-meta_data[order(meta_data$order,decreasing = T),]
meta_data<-meta_data[order(meta_data$region),]
#meta_data<-rbind(meta_data[1:229,],meta_data[match(paste0(dataset_slice,'_other'),meta_data$dataSlice),])
all_RBP_TF<-all_RBP_TF[,meta_data$dataSlice]

# stem_mean<-stem_mean[colnames(all_RBP_TF),]
# stem_mean$x[is.na(stem_mean$x)]<-0

ha <- HeatmapAnnotation(df = data.frame(cancer = meta_data$cancer),
                        col = cancer_data)

la <- rowAnnotation(df = data.frame(geneType=c(rep('TF',length(TF)),
                                               rep('RBP',length(RBP)))),
                    col = list(geneType=c("TF"="#99CC99","RBP"="#F39C67")
                    ))

heat_data<-t(scale(t(all_RBP_TF)))
heat_data<-all_RBP_TF
range(heat_data,na.rm = T)
col_fun <- circlize::colorRamp2(
  seq(from=0,to=2,length.out=3), 
  c("white","#FCAD87","#9F001C")
#?seq

p_Heatmap<-Heatmap(heat_data,
                   col = col_fun,
                   top_annotation =ha,#####顶部注释
                   left_annotation=la,#####左侧注释
                   show_row_names = T,#####不显示行名
                   show_column_names =F,#####不显示列名
                   cluster_rows =F,
                   cluster_columns = F,
                   #column_split = factor(meta_data$region, unique(meta_data$region)),
                   #column_gap = unit(1, "mm"),
                   row_split = c(rep('1TF',length(TF)),
                                 rep('RBP',length(RBP))),
                   row_gap = unit(1, "mm"),
                   name="scale of exp",
                   column_title='TF&RBP',
                   use_raster=F,
                   row_names_gp = gpar(fontsize = 7))
print(p_Heatmap)
#?Heatmap

pdf('/pic/RBP_TF/all_slice_heatmap.pdf',width = 10,height = 7.1)
print(p_Heatmap)
dev.off()











