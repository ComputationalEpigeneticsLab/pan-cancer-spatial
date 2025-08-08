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
###数据来源148_allSlice_TF&RBP.R

RBP<-c('ELAVL1', 'HNRNPF', 'HNRNPH1', 'HNRNPH3', 'MBNL1', 'PUM1', 'TIAL1')
TF<-c('FLI1','FOXA1','MXI1', 'RARA','FOXA2','RUNX2','E2F3',
      'HDAC2','FOS','CTCF','JUN','MYC','E2F4','REST','HES1','MAFB','JUND','ELK1','ETV4',
      'AHR','IRF2','MAFF','CEBPD','ETV7','TBP','TCF21','TRIM28','GABPA','MLX','MAX','SREBF2',
      'CEBPA','GATA2','FOXO3','ETV1','SOX9','BCLAF1','TCF12', 'TEAD4','ZNF274',
      'YY1','ETS2','PHF8','FOXO4','KLF5','FOXO1','RXRB','CREB1','SP4','GATA6','POLR3A')
all_RBP<-read.delim('E:/Mirror/ST_analysis/data/pic_data/TF_RBP/all_RBP_quantile_mean.txt',stringsAsFactors = F,check.names = F)
all_TF<-read.delim('E:/Mirror/ST_analysis/data/pic_data/TF_RBP/all_TF_quantile_mean.txt',stringsAsFactors = F,check.names = F)

all_RBP_TF<-rbind(all_TF,all_RBP)
all_RBP_TF<-all_RBP_TF[,c(grep('core',colnames(all_RBP_TF)))]

###每个切片的干性
all_slice_stem<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_stem.txt',
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
                        col = list(cancer=c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                                            'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                                            'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                                            'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
                                   
                        ))
# 'Boundary'='#FF9900','Dispersion'='#CD5A5A',
# 'Immune'='#5477AF','Normal'='#669933'

la <- rowAnnotation(df = data.frame(geneType=c(rep('TF',length(TF)),
                                               rep('RBP',length(RBP)))),
                    col = list(geneType=c("TF"="#99CC99","RBP"="#F39C67")
                    ))

heat_data<-t(scale(t(all_RBP_TF)))
heat_data<-all_RBP_TF
range(heat_data,na.rm = T)
col_fun <- circlize::colorRamp2(
  seq(from=0,to=2,length.out=3), 
  c("white","#FCAD87","#9F001C")##"#0053A6","#7CB9DB", 
)
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

pdf('E:/Mirror/ST_analysis/pic/RBP_TF/all_slice_heatmap.pdf',width = 10,height = 7.1)
print(p_Heatmap)
dev.off()





####单独绘制LUAD的
LUAD_meta<-meta_data[grep('luad',meta_data$dataSlice),]
LUAD_data<-all_RBP_TF[,LUAD_meta$dataSlice]
colnames(LUAD_data)<-unlist(lapply(strsplit(colnames(LUAD_data),'_co'),function(x)x[1]))
NA_row<-apply(LUAD_data,1,function(x){
  length(which(is.na(x)==T))
})
LUAD_data<-LUAD_data[-which(NA_row>0),]
low_row<-match(c('ZNF274','MXI1','PHF8','TBP','POLR3A','FOXA1','REST','IRF2','FOXO1','GABPA','TCF12','E2F3','FLI1'),
               rownames(LUAD_data))
LUAD_data<-LUAD_data[-low_row,]

gene_order<-apply(LUAD_data,1,function(x){sum(x,na.rm = T)})
gene_order<-data.frame(type=c(rep('TF',length(intersect(TF,rownames(LUAD_data)))),
                              rep('RBP',length(intersect(RBP,rownames(LUAD_data))))),
                       order=gene_order,
                       gene=names(gene_order))
gene_order<-gene_order[order(gene_order$order,decreasing = T),]
gene_order<-gene_order[order(gene_order$type,decreasing = T),]
LUAD_data<-LUAD_data[gene_order$gene,]

la <- rowAnnotation(df = data.frame(geneType=c(rep('TF',length(intersect(TF,rownames(LUAD_data)))),
                                               rep('RBP',length(intersect(RBP,rownames(LUAD_data)))))),
                    col = list(geneType=c("TF"="#99CC99","RBP"="#F39C67")
                    ))

heat_data2<-t(scale(t(LUAD_data)))
heat_data2<-LUAD_data
range(heat_data2,na.rm = T)
col_fun <- circlize::colorRamp2(
  c(0,0.5,2), 
  c("white","#FCAD87","#9F001C")##"#0053A6","#7CB9DB", 
)
#?seq

p_Heatmap<-Heatmap(heat_data2,
                   col = col_fun,
                   #top_annotation =ha,#####顶部注释
                   left_annotation=la,#####左侧注释
                   show_row_names = T,#####不显示行名
                   show_column_names =T,#####不显示列名
                   cluster_rows =F,
                   cluster_columns = F,
                   #column_split = factor(meta_data$region, unique(meta_data$region)),
                   #column_gap = unit(1, "mm"),
                   row_split = c(rep('1TF',length(intersect(TF,rownames(LUAD_data)))),
                                 rep('RBP',length(intersect(RBP,rownames(LUAD_data))))),
                   row_gap = unit(1, "mm"),
                   name="exp",
                   column_title='TF&RBP',
                   use_raster=F,
                   row_names_gp = gpar(fontsize = 9))
print(p_Heatmap)

pdf('E:/Mirror/ST_analysis/pic/RBP_TF/all_slice_heatmap_LUAD.pdf',width = 4,height = 8)
print(p_Heatmap)
dev.off()







