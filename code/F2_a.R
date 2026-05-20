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

all_slice_stem = read.table('/stemness/all_slice_stem.txt',header = T,row.names = 1,sep = '\t');gc()
all_slice_spottype = read.table('/stemness/all_slice_spottype.txt',header = T,row.names = 1,sep = '\t');gc()
table(rownames(all_slice_spottype)%in%rownames(all_slice_stem))
table(all_slice_spottype$FinalLocalType)
all_slice_spottype[all_slice_spottype$FinalLocalType%in%c('Budding'),'FinalLocalType']='Dispersion'
all_slice_spottype[all_slice_spottype$FinalLocalType%in%c('Immune','Normal'),'FinalLocalType']='Stromal'
all_slice_stem = all_slice_stem[match(rownames(all_slice_spottype),rownames(all_slice_stem)),]
identical(rownames(all_slice_spottype),rownames(all_slice_stem))
all_slice_spottype$score = all_slice_stem$score
rm(all_slice_stem);gc()
all_slice_spottype = all_slice_spottype[-c(which(is.na(all_slice_spottype$score))),]

slice_stem<-aggregate(all_slice_spottype$score,by=list(all_slice_spottype$dataset_slice,all_slice_spottype$FinalLocalType),mean)
#slice_stem<-slice_stem[which(slice_stem$Group.2=='Core'|slice_stem$Group.2=='Boundary'|slice_stem$Group.2=='Dispersion'),]
slice_stem<-reshape2::dcast(slice_stem,fun.aggregate = sum,Group.1~Group.2)
rownames(slice_stem)<-slice_stem$Group.1
slice_stem<-slice_stem[,-1]
slice_stem<-slice_stem[,c('Core','Boundary','Dispersion','Stromal')]

cancer = all_slice_spottype %>% distinct(cancer,dataset_slice) %>% as.data.frame()
table(rownames(slice_stem)%in%cancer$dataset_slice)
slice_stem$cancer = cancer[match(rownames(slice_stem),cancer$dataset_slice),]$cancer
slice_stem<-slice_stem[order(slice_stem$Core,decreasing = T),]
slice_stem<-slice_stem[order(slice_stem$cancer),]

slice_stem = slice_stem[-c(which(slice_stem$cancer%in%c('GBM','BRCA','PDAC') & (slice_stem$Stromal>slice_stem$Core))),]

color<-read.delim('/35cancer_color.txt',stringsAsFactors = F,check.names = F)
cancer_color=color$color
names(cancer_color)=color$cancer
la <- rowAnnotation(df = data.frame(cancer=slice_stem$cancer),
                    col = list(cancer=cancer_color
                    ))
ha <- HeatmapAnnotation(df = data.frame(cancer = slice_stem$cancer),
                        col = list(cancer=cancer_color
                                   
                        ))
plot_data<-slice_stem[,1:4]%>%as.matrix()
plot_data<-t(plot_data)
range(plot_data)
col_fun <- circlize::colorRamp2(
  c(0,0.1,0.3,0.5,0.7),        
  #seq(from=0.4,to=0.85,length.out=3), 
  c("#04040b",'#5b2177',"#b93c6d","#eb7d60","#f6f0b7")
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

pdf('/stemness/stemness_all_heatmap.pdf',width = 12,height = 2)
print(p_Heatmap)
dev.off()
