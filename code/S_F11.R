###使用SpatialIDM对挑选的LR对进行验证
library(scatterpie)
library(ggpubr)
library(png)
library(jsonlite)
library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)



###COL1A2_SDC4
bdy_CAF_LR1_slice<-c('brca01/slice3','brca11/slice1','brca14/slice1','brca26/slice7','brca26/slice8','cesc01/slice1',
                     'crc02/slice1','crc04/slice1','gbm04/slice10','lihc02/slice1','lihc02/slice12','lihc02/slice13',
                     'lihc02/slice2','lihc02/slice3','lihc02/slice5','luad01/slice2','oscc01/slice4','prad02/slice1')
###COL1A1_SDC1
bdy_CAF_LR2_slice<-c('brca01/slice3','brca07/slice1','brca11/slice1','brca14/slice1','brca26/slice7','cesc01/slice1',
                     'crc02/slice1','lihc02/slice1','lihc02/slice12','lihc02/slice13','lihc02/slice2','lihc02/slice3',
                     'lihc02/slice5','lihc02/slice6','lihc03/slice1','luad01/slice2','oscc01/slice4','ovca07/slice1')

###SPP1_ITGB1
bdy_TAM_LR1_slice<-c('brca01/slice3','brca08/slice4','brca26/slice1','brca26/slice6','cesc02/slice2','lihc02/slice4',
                     'lihc02/slice6','lihc02/slice7','lihc02/slice9','ovca03/slice1','ovca04/slice1','ovca06/slice1',
                     'ovca06/slice4','ovca07/slice5','rcc01/slice1','rcc01/slice2','rcc03/slice1','pdac03/slice3')
###LAMB2_CD44
bdy_TAM_LR2_slice<-c('brca01/slice3','brca08/slice4','brca11/slice1','brca12/slice1','brca13/slice1','brca14/slice1',
                     'brca14/slice2','brca15/slice2','lihc02/slice7','lihc02/slice8','lihc02/slice9','luad01/slice1',
                     'luad01/slice2','luad03/slice1','ovca03/slice1','ovca06/slice1','ovca07/slice5','rcc03/slice1')

###FN1_CD44
matrixCAF_LR1_slice<-c('brca08/slice3','brca15/slice1','brca15/slice2','brca26/slice1','brca26/slice6','crc01/slice1',
                       'crc05/slice1','crc06/slice1','crc06/slice2','crc07/slice1','crc07/slice2','crc08/slice2',
                       'gbm05/slice3','lihc02/slice14','lihc02/slice4','luad01/slice1','luad02/slice2','skcm13/slice1')
###ICAM1_ITGB2
matrixCAF_LR2_slice<-c('brca26/slice8','crc07/slice1','crc07/slice2','crc08/slice2','lihc02/slice7','lihc02/slice14',
                       'lihc02/slice2','lihc02/slice3','lihc02/slice4','lihc02/slice5','lihc02/slice9','lihc03/slice7',
                       'lihc03/slice8','luad01/slice1','luad02/slice2','oscc01/slice4','pdac03/slice1','skcm12/slice1')

###CD99_PILRB
matrixTAM_LR1_slice<-c('brca08/slice4','brca11/slice1','brca14/slice2','brca15/slice1','brca15/slice2','brca26/slice6',
                       'hn-as02/slice2','hn-as03/slice1','lihc02/slice10','lihc02/slice18','lihc02/slice7','ovca03/slice1',
                       'ovca07/slice5','pcnsl01/slice3','pdac03/slice1','pdac03/slice2','rcc01/slice4','rcc03/slice1')
###THBS1_SDC1
matrixTAM_LR2_slice<-c('brca08/slice4','brca11/slice1','brca12/slice1','brca14/slice2','brca15/slice1','brca15/slice2',
                       'brca26/slice1','brca26/slice6','hn-as01/slice2','hn-as01/slice4','hn-as02/slice2','lihc02/slice18',
                       'lihc02/slice4','lihc02/slice7','luad02/slice2','luad03/slice1','pdac03/slice1','rcc03/slice1')

'-−'
all_slice<-c(bdy_CAF_LR1_slice,bdy_CAF_LR2_slice,bdy_TAM_LR1_slice,bdy_TAM_LR2_slice,
             matrixCAF_LR1_slice,matrixCAF_LR2_slice,matrixTAM_LR1_slice,matrixTAM_LR2_slice) %>% unique()

dir_pic<-'E:/Mirror/ST_analysis/pic/cellchat/LR_express/select_LR/IDM/'
#################################################################################################################
###bdy_CAF
select_slice<-bdy_CAF_LR1_slice ##bdy_CAF_LR1_slice  ###COL1A2 _ SDC4   ###COL1A1 _ SDC1
L_gene<-'COL1A2';R_gene<-'SDC4'##bdy_CAF_LR1_slice
select_slice<-bdy_CAF_LR2_slice
L_gene<-'COL1A1';R_gene<-'SDC1'##bdy_CAF_LR2_slice

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_z_score.txt')### .h5ad_local_I.txt
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')

pdf(paste0(dir_pic,'Zscore_bdyCAF_',L_gene,'_',R_gene,'2.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'
  
  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
  
  near_spot<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'CAF'])
  p_data$type[p_data$cell_name%in%near_spot]<-'step1_spot'
  p_data$type[p_data$type%in%NA]<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Zscore')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()



#################################################################################################################
###bdy_TAM
select_slice<-bdy_TAM_LR1_slice  ##bdy_TAM_LR1_slice  ######SPP1 _ ITGB1   ###LAMB2 _ CD44
L_gene<-'SPP1';R_gene<-'ITGB1'
select_slice<-bdy_TAM_LR2_slice 
L_gene<-'LAMB2';R_gene<-'CD44'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_z_score.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')


pdf(paste0(dir_pic,'Zscore_bdyTAM_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'
  
  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
  
  near_spot<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'TAM'])
  p_data$type[p_data$cell_name%in%near_spot]<-'step1_spot'
  p_data$type[p_data$type%in%NA]<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Zscore')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()



#################################################################################################################
###matrix_CAF
select_slice<-matrixCAF_LR1_slice  ###matrixCAF_LR1_slice  ######FN1 _ CD44   ###ICAM1 _ ITGB2
L_gene<-'FN1';R_gene<-'CD44'
select_slice<-matrixCAF_LR2_slice
L_gene<-'ICAM1';R_gene<-'ITGB2'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_z_score.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')



pdf(paste0(dir_pic,'Zscore_MatrixCAF_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other_spot'
  p_data$type[p_data$type%in%'TAM']<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Zscore')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()



#################################################################################################################
###matrix_TAM
select_slice<-matrixTAM_LR1_slice  ##matrixTAM_LR1_slice  ######CD99 _ PILRB   ###THBS1 _ SDC1
L_gene<-'CD99';R_gene<-'PILRB'
select_slice<-matrixTAM_LR2_slice
L_gene<-'THBS1';R_gene<-'SDC1'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_z_score.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')

pdf(paste0(dir_pic,'Zscore_MatrixTAM_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other_spot'
  p_data$type[p_data$type%in%'CAF']<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Zscore')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()



###每组挑两对LR 每对LR挑6个切片 先绘制箱式图看差异
library(ggpubr)
library(png)
library(jsonlite)
library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)



dir_pic<-'E:/Mirror/ST_analysis/pic/cellchat/LR_express/select_LR/IDM2/'
#################################################################################################################
###bdy_CAF
###COL1A1 _???_ ITGB1
bdy_CAF_LR1_slice<-c('brca07/slice1',
                     'prad02/slice1',
                     'brca26/slice8',
                     'brca01/slice3',
                     'crc02/slice1',
                     'brca14/slice1')

###COL1A2 _ SDC4
bdy_CAF_LR2_slice<-c('brca07/slice1',
                     'brca14/slice1',
                     'brca01/slice3',
                     'luad01/slice2',
                     'brca26/slice7',
                     'lihc02/slice13')

select_slice<-bdy_CAF_LR1_slice 
L_gene<-'COL1A1';R_gene<-'ITGB1'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_I.txt')### .h5ad_local_I.txt   .h5ad_local_z_score.txt
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/' 
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')

pdf(paste0(dir_pic,'Moran_bdyCAF_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  #Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'
  
  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
  
  near_spot<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'CAF'])
  p_data$type[p_data$cell_name%in%near_spot]<-'step1_spot'
  p_data$type[p_data$type%in%NA]<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        # if(class(Moran_score)=='data.frame'){
        #   Moran_score<-apply(Moran_score,1,mean)
        #   p_data$Moran<-Moran_score
        # }else{
        #   p_data$Moran<-Moran_score
        # }
        
        if(length(grep('ITGA10_ITGB1',colnames(Moran_score)))>0){
          p_data$Moran<-Moran_score[,grep('ITGA10_ITGB1',colnames(Moran_score))]
        }else{
          p_data$Moran<-0
        }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Moran')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()




#################################################################################################################
###bdy_TAM
###SPP1 _???_ ITGB1
bdy_TAM_LR2_slice<-c('brca08/slice4',
                     'brca14/slice1',
                     'lihc02/slice7',
                     'ovca06/slice4',
                     'rcc01/slice2',
                     'ovca06/slice1')
# ###ANGPTL4 _ ITGB1
# bdy_TAM_LR2_slice<-c('brca13/slice1',#
#                      'brca08/slice4',#
#                      'lihc02/slice7',#
#                      'brca12/slice1',#
#                      'brca15/slice2',#
#                      'brca14/slice2'
#                      )
###ANGPTL4 _ ITGA5
bdy_TAM_LR3_slice<-c('brca13/slice1',#
                     'brca08/slice4',#
                     'lihc02/slice8',#
                     'lihc02/slice7',#
                     'brca12/slice1',#
                     'brca15/slice2')

select_slice<-bdy_TAM_LR3_slice 
L_gene<-'ANGPTL4';R_gene<-'ITGA5'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_I.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')


pdf(paste0(dir_pic,'Moran_bdyTAM_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  #Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'
  
  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
  
  near_spot<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'TAM'])
  p_data$type[p_data$cell_name%in%near_spot]<-'step1_spot'
  p_data$type[p_data$type%in%NA]<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        # if(length(grep('ITGA8_ITGB1',colnames(Moran_score)))>0){
        #   p_data$Moran<-Moran_score[,grep('ITGA8_ITGB1',colnames(Moran_score))]
        # }else{
        #   p_data$Moran<-Moran_score[,1]
        # }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Moran')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()



#################################################################################################################
###matrix_CAF
###FN1 _ CD44
matrixCAF_LR1_slice<-c('crc07/slice1',#
                       'lihc02/slice14',#
                       'skcm13/slice1',#
                       'brca08/slice3',#
                       'luad01/slice1',#
                       'crc01/slice1')

##ICAM1 _???_ ITGB2
matrixCAF_LR3_slice<-c('skcm13/slice1',#
                       'lihc03/slice8',#
                       'skcm12/slice1',#
                       'brca26/slice6',#
                       'lihc03/slice7',#
                       'crc06/slice1')
select_slice<-matrixCAF_LR3_slice
L_gene<-'ICAM1';R_gene<-'ITGB2'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_I.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')


pdf(paste0(dir_pic,'Moran_MatrixCAF_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  #Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other_spot'
  p_data$type[p_data$type%in%'TAM']<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        # if(length(grep('ITGA8_ITGB1',colnames(Moran_score)))>0){
        #   p_data$Moran<-Moran_score[,grep('ITGA8_ITGB1',colnames(Moran_score))]
        # }else{
        #   p_data$Moran<-Moran_score[,1]
        # }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Moran')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()





#################################################################################################################
###matrix_TAM
###THBS1 _ CD47
matrixTAM_LR1_slice<-c('pdac03/slice1',#
                       'rcc03/slice1',#
                       'hn-as02/slice2',#
                       'luad03/slice1',#
                       'hn-as03/slice1',#
                       'rcc01/slice4')
###THBS1 _ SDC1
matrixTAM_LR2_slice<-c('pdac03/slice2',#
                       'brca11/slice1',#
                       'brca08/slice4',#
                       'rcc01/slice4',#
                       'brca15/slice1',#
                       'luad03/slice1')

select_slice<-matrixTAM_LR2_slice
L_gene<-'THBS1';R_gene<-'SDC1'

dir_IDM<-'E:/Mirror/ST_analysis/data/10X Visium/SpatailDM/re_noFileter/'
file_IDM<-paste0('adata_',gsub('/','_',select_slice),'.h5ad_local_I.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(select_slice,'_BdyNearCore.txt')
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
file_RCTD<-paste0(select_slice,'_RCTD_celltype.txt')


pdf(paste0(dir_pic,'Moran_MatrixTAM_',L_gene,'_',R_gene,'.pdf'),width = 3,height = 4)
for(i in 1:length(select_slice)){
  #i=1
  Moran_score<-read.delim(paste0(dir_IDM,file_IDM[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  #Moran_score<-as.data.frame(t(Moran_score))
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),]
  
  p_data<-st_near[,1:4]
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other_spot'
  p_data$type[p_data$type%in%'CAF']<-'other_spot'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  p_data$type2<-p_data$type
  p_data$type2[!p_data$type%in%'other_spot']<-'aim_spot'
  Moran_score<-Moran_score[rownames(p_data),]
  
  if('aim_spot'%in%p_data$type2&&'other_spot'%in%p_data$type2){
    if(length(grep(L_gene,colnames(Moran_score)))>0){
      Moran_score<-Moran_score[,grep(L_gene,colnames(Moran_score))]
      if(length(grep(R_gene,colnames(Moran_score)))>0){
        Moran_score<-Moran_score[,grep(R_gene,colnames(Moran_score))]
        
        if(class(Moran_score)=='data.frame'){
          Moran_score<-apply(Moran_score,1,mean)
          p_data$Moran<-Moran_score
        }else{
          p_data$Moran<-Moran_score
        }
        
        # if(length(grep('ITGA8_ITGB1',colnames(Moran_score)))>0){
        #   p_data$Moran<-Moran_score[,grep('ITGA8_ITGB1',colnames(Moran_score))]
        # }else{
        #   p_data$Moran<-Moran_score[,1]
        # }
        
        plot_data<-lapply(unique(p_data$type2),function(x){#x='aim_spot'
          xx<-p_data$Moran[which(p_data$type2==x)]
          Q <- quantile(xx, probs=c(.25, .75), na.rm = FALSE)###计算25%与75%的值
          iqr <- IQR(xx)###计算四分位间距
          up <- Q[2]+1.5*iqr # Upper Range 上限
          low<- Q[1]-1.5*iqr # Lower Range 下限
          eliminated<- xx[which(xx > (Q[1] - 1.5*iqr) & xx < (Q[2]+1.5*iqr))]
          eliminated<-eliminated-min(eliminated)
          eliminated<-eliminated/max(eliminated)
          xx_data<-data.frame(type2=x,Moran=eliminated)
          return(xx_data)
        })
        plot_data<-do.call(rbind,plot_data)
        
        # plot_data$Moran<-plot_data$Moran-min(plot_data$Moran)
        # plot_data$Moran<-plot_data$Moran/max(plot_data$Moran)
        
        p_box<-ggplot(plot_data, aes(x = type2, y = Moran,fill=type2))+ 
          # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
          stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
          geom_boxplot(aes(fill = type2), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
                       position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
          theme_classic(base_size = 12)+
          theme(axis.text = element_text(color = 'black'))+
          scale_fill_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          scale_color_manual(values = c("aim_spot"="#179192","other_spot"="#d2e7cf"))+
          ggtitle(paste0(select_slice[i],'_',L_gene,'_',R_gene))+
          theme(plot.title = element_text(hjust = 0))+
          theme(plot.title = element_text(size = 12))+
          ylab('scale_Moran')+
          #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
          stat_compare_means(aes(label = ..p.signif..),label.x = 1.5,label.y =1)+
          #theme_bw()+
          theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),legend.position = 'none',
                axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
        print(p_box)
        
      }
    }
  }
  
}
dev.off()










