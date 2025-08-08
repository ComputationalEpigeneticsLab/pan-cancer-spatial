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
library(rlang)
library(ggplot2)
library(ggraph)
library(ggpubr)
library(psych)
#########################################################################################################################
####按每种染色体的基因变异均值绘制
dir_chr<-'E:/Mirror/ST_analysis/data/10X Visium/CNV_stem/chr_stat/'
file_chr<-list.files(pattern = 'txt',path = dir_chr)
data_slice<-unlist(lapply(strsplit(file_chr,'_chr'),function(x)x[1]))


i=1
chr_data<-read.delim(paste0(dir_chr,file_chr[i]),stringsAsFactors = F,check.names = F)
chr_data$chr_gene<-paste0(chr_data$chr,'_of_',chr_data$gene)
colnames(chr_data)[1:2]<-paste0(data_slice[i],'_of_',colnames(chr_data)[1:2])

all_chr_data<-chr_data[,c(1,2,5)]
for(i in 2:length(file_chr)){
  #i=2
  chr_data<-read.delim(paste0(dir_chr,file_chr[i]),stringsAsFactors = F,check.names = F)
  chr_data$chr_gene<-paste0(chr_data$chr,'_of_',chr_data$gene)
  colnames(chr_data)[1:2]<-paste0(data_slice[i],'_of_',colnames(chr_data)[1:2])
  all_chr_data<-merge(all_chr_data,chr_data[,c(1,2,5)],by='chr_gene',all=T)
}

rownames(all_chr_data)<-all_chr_data$chr_gene
all_chr_data<-all_chr_data[,-1]
# all_chr_data<-as.matrix(all_chr_data)
# all_chr_data[is.na(all_chr_data)]<-1
table(unlist(lapply(strsplit(rownames(all_chr_data),'_of_'),function(x)x[1])))

chr_data_core<-all_chr_data[,grep('core',colnames(all_chr_data))]
chr_data_BdyBud<-all_chr_data[,grep('bdy',colnames(all_chr_data))]

select_slice1<-lapply(unique(chr_allMean$chr),function(x){#x='chr1'
  aa_chr<-chr_allMean[chr_allMean$chr%in%x,]
  aa_chr<-aa_chr[which(aa_chr$p<0.01&aa_chr$diff>0),]
  #aa_chr<-aa_chr[order(aa_chr$diff,decreasing = T),]
  return(aa_chr$slices)
})
names(select_slice1)<-unique(chr_allMean$chr)
# select_slice1<-as.matrix(t(do.call(cbind,select_slice1)))
# select_slice1<-reshape2::melt(select_slice1)
select_slice1[['chr2']]

color_data<-c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
              'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
              'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
              'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')

#pdf(paste0('E:/Mirror/ST_analysis/pic/clone/CNA_chr/','chr_heatmap_select2','.pdf'),width = 10,height = 6)
pdf(paste0('E:/Mirror/ST_analysis/pic/clone/CNA_chr/','chr_box_select4','.pdf'),width = 4,height = 4)
for(j in paste0('chr',1:22)){
  #j='chr1'
  s_slice<-select_slice1[[j]]
  chr_single<-cbind(chr_data_core[grep(paste0(j,'_of_'),rownames(chr_data_core)),paste0(s_slice,'_of_core_mean')],
                    chr_data_BdyBud[grep(paste0(j,'_of_'),rownames(chr_data_BdyBud)),paste0(s_slice,'_of_bud_bdy_mean')])
  rownames(chr_single)<-unlist(lapply(strsplit(rownames(chr_single),'_of_'),function(x)x[2]))
  chr_single<-chr_single[which(apply(chr_single,1,function(x) sum(is.na(x))/ncol(chr_single) )<0.5),]
  
  row_order<-data.frame(gene=rownames(chr_single),
                        order=apply(chr_single,1,function(x) sum(x,na.rm = T))
  )
  row_order<-row_order[order(row_order$order,decreasing = T),]
  chr_single<-chr_single[row_order$gene,]
  
  col_data<-data.frame(slice=unlist(lapply(strsplit(colnames(chr_single),'_of_'),function(x)x[1])),
                       spotType=rep(c('Core','BdyBud'),each=ncol(chr_single)/2),
                       order=apply(chr_single,2,function(x) sum(x,na.rm = T)),
                       row.names = colnames(chr_single))
  col_data$cancer<-unlist(lapply(strsplit(col_data$slice,'_'),function(x)x[1]))
  col_data$cancer<-substr(col_data$cancer,1,nchar(col_data$cancer)-2) %>% toupper()
  col_data<-col_data[order(col_data$order,decreasing = T),]
  col_data<-col_data[order(col_data$spotType,decreasing = T),]
  chr_single<-chr_single[,rownames(col_data)]
  
  
  # annotation_col = data.frame(
  #   spotType = col_data$spotType,
  #   cancer=col_data$cancer
  # )
  # rownames(annotation_col) = rownames(col_data)
  # head(annotation_col)
  # 
  # ann_colors = list(
  #   spotType=c('Core'='#990033','BdyBud'='#FF9900'),
  #   cancer=color_data[unique(col_data$cancer)]
  # )
  
  p_data<-chr_single
  range(p_data,na.rm = T)
  p_data<-p_data[which(apply(p_data,1,function(x) sum(is.na(x))/ncol(p_data) )<0.5),]
  p_data<-p_data[,which(apply(p_data,2,function(x) sum(is.na(x))/nrow(p_data) )<1)]
  
  col_data<-col_data[colnames(p_data),]
  col_data$chrGeneExp<-apply(p_data,2,function(x) mean(x,na.rm=T))
  #?mean
  #sum(is.na(c(1,3,9,NA,0,NA,8,NA,NA)))
  ylim1<-lapply(unique(col_data$spotType),function(x){#x='Core'
    xx<-col_data$chrGeneExp[which(col_data$spotType==x)]
    Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
    iqr <- IQR(xx)###计算四分位间距
    up <- Q[2]+1.5*iqr # Upper Range 上限
    low<- Q[1]-1.5*iqr # Lower Range 下限
    return(c(up,low))
  })
  ylim1<-do.call(rbind,ylim1)
  
  col_data$spotType<-factor(col_data$spotType, levels = c('Core','BdyBud'))
  p_box<-ggplot(col_data, aes(x = spotType, y = chrGeneExp,fill=spotType))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = spotType), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c('Core'='#d62d28','BdyBud'='#f6b86d'))+
    scale_color_manual(values = c('Core'='#d62d28','BdyBud'='#f6b86d'))+
    ggtitle(j)+
    ylab('chrCNVmean')+
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 1.4,label.y = max(ylim1[,1])*1.02)+
    #stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
    #theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  
  ppx<-p_box + coord_cartesian(ylim = c(min(ylim1[,2]),max(ylim1[,1])*1.04))
  print(ppx)
  
  
  #p_data<-t(apply(p_data,1,function(x) (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T))))
  # range(p_data,na.rm = T)
  # col_fun <- colorRamp2(
  #   c(0.9, 1, 1.1), 
  #   c("#377EB8","#F0F0F0","#E41A1C"))
  # #?seq
  # 
  # p<-pheatmap(as.matrix(p_data), 
  #             scale = "none",
  #             color=col_fun,
  #             #border_color = 'white',
  #             show_rownames = F,
  #             show_colnames = F,
  #             cluster_rows = F,
  #             cluster_cols = F,
  #             fontsize = 10,
  #             gaps_col = c(ncol(chr_single)/2),
  #             main = j,
  #             #breaks = bk,
  #             annotation_colors = ann_colors,
  #             annotation_col = annotation_col,
  #             name = 'CNV',
  #             use_raster=F
  # )
  # 
  # print(p)
  
  
  
}
dev.off()