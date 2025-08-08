####core bdy dis区的干性得分热图
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)



all_slice_stem<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_stem.txt',
                           stringsAsFactors = F,check.names = F)

slice_stem<-aggregate(all_slice_stem$CytoTRACE,by=list(all_slice_stem$dataset_slice,all_slice_stem$localType),mean)
slice_stem<-slice_stem[which(slice_stem$Group.2=='Core'|slice_stem$Group.2=='Boundary'|slice_stem$Group.2=='Dispersion'),]
slice_stem<-reshape2::dcast(slice_stem,fun.aggregate = sum,Group.1~Group.2)
rownames(slice_stem)<-slice_stem$Group.1
slice_stem<-slice_stem[,-1]
slice_stem<-slice_stem[,c('Core','Boundary','Dispersion')]
slice_stem$cancer<-unlist(lapply(strsplit(rownames(slice_stem),'/'),function(x) x[1]))
slice_stem$cancer<-substr(slice_stem$cancer,1,(nchar(slice_stem$cancer)-2))%>%toupper()
unique(slice_stem$cancer)
slice_stem<-slice_stem[order(slice_stem$Core,decreasing = T),]
slice_stem<-slice_stem[order(slice_stem$cancer),]

la <- rowAnnotation(df = data.frame(cancer=slice_stem$cancer),
                    col = list(cancer=c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                                        'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                                        'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                                        'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
                    ))
ha <- HeatmapAnnotation(df = data.frame(cancer = slice_stem$cancer),
                        col = list(cancer=c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                                            'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                                            'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                                            'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
                                   
                        ))
plot_data<-slice_stem[,1:3]%>%as.matrix()
plot_data<-t(plot_data)
range(plot_data)
col_fun <- circlize::colorRamp2(
  c(0,0.4,0.6,0.8),        
  #seq(from=0.0,to=0.85,length.out=3), 
  c("white",'#DCCECF',"#DC9DA3","#770015")
)

p_Heatmap<-Heatmap(plot_data,
                   col = col_fun,
                   top_annotation =ha,#####顶部注释
                   #left_annotation=la,#####左侧注释
                   show_row_names = T,#####不显示行名
                   show_column_names =F,#####不显示列名
                   cluster_rows =F,
                   cluster_columns = F,
                   # column_split = rep(c("M3_like", "other"),c(277,734)),
                   # row_split = rep(c("1CMP","2myeloid","imm_1","imm_2","imm_3","TME"),
                   #                 c(1,3,7,7,2,22)) ,
                   name="stemness",
                   use_raster=F,
                   row_names_gp = gpar(fontsize = 10))
print(p_Heatmap)

pdf('E:/Mirror/ST_analysis/pic/re/3/stemness_all_2.pdf',width = 12,height = 2)
print(p_Heatmap)
dev.off()


###挑选一些切片 每种癌症两个
select_stem<-slice_stem[c('brca11/slice1','brca08/slice1','cesc02/slice2','cesc02/slice3','crc05/slice2','crc06/slice2',
                          'cscc01/slice4','cscc02/slice3','gbm03/slice9','gbm05/slice14','gist01/slice1','gist01/slice2',
                          'hgsc01/slice2','hgsc01/slice4','hn-as02/slice3','hn-as03/slice2','ipmn01/slice10','ipmn01/slice11',
                          'lihc01/slice1','lihc02/slice4','luad02/slice1','luad03/slice1','mibc01/slice1','mibc01/slice3',
                          'oscc01/slice5','oscc01/slice7','ovca06/slice2','ovca06/slice5','pcnsl01/slice1','pcnsl01/slice3',
                          'pdac03/slice2','pdac03/slice3','prad02/slice1','prad07/slice1','rcc01/slice1','rcc01/slice5',
                          'skcm12/slice1','skcm13/slice1'),]

la <- rowAnnotation(df = data.frame(cancer=select_stem$cancer),
                    col = list(cancer=c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                                        'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                                        'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                                        'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
                    ))
plot_data<-select_stem[,1:3]%>%as.matrix()
range(plot_data)
col_fun <- circlize::colorRamp2(
  c(0,0.3,0.7,0.8),        
  #seq(from=0.4,to=0.85,length.out=3), 
  c("white",'#DCCECF',"#DC8890","#770015")
)

p_Heatmap<-Heatmap(plot_data,
                   col = col_fun,
                   #top_annotation =ha,#####顶部注释
                   left_annotation=la,#####左侧注释
                   show_row_names = F,#####不显示行名
                   show_column_names =T,#####不显示列名
                   cluster_rows =F,
                   cluster_columns = F,
                   width = 5,
                   height = 5,
                   # column_split = rep(c("M3_like", "other"),c(277,734)),
                   # row_split = rep(c("1CMP","2myeloid","imm_1","imm_2","imm_3","TME"),
                   #                 c(1,3,7,7,2,22)) ,
                   #border = 'white',
                   rect_gp = gpar(col = "white", lwd = 2),
                   name="stemness",
                   use_raster=F,
                   row_names_gp = gpar(fontsize = 10))
print(p_Heatmap)

pdf('E:/Mirror/ST_analysis/pic/re/3/stemness_select.pdf',width = 2,height = 12)
print(p_Heatmap)
dev.off()



?Heatmap


annotation_row = data.frame(
  cancer = select_stem$cancer
)
rownames(annotation_row) = rownames(select_stem)
head(annotation_row)

ann_colors = list(
  cancer=c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
           'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
           'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
           'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
)

p<-pheatmap(plot_data, 
            scale = "none",
            color=col_fun,
            cluster_rows = F,
            cluster_cols = F,
            show_rownames = F,
            show_colnames = T,
            border='white',
            cellwidth=10,
            cellheight=10,
            name = 'stemness',
            #breaks = bk,
            annotation_row = annotation_row,
            annotation_colors = ann_colors,
            use_raster=F
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/re/3/stemness_select2.pdf',width = 2,height = 12)
print(p_Heatmap)
dev.off()

?pheatmap











