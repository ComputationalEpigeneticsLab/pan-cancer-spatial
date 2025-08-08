###每种癌症的step1的CAF+TAM与其他免疫细胞加一起比较：FC值和秩和检验p值
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)


dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_near<-list.files(pattern = '_nearSpotStep1to3.txt',path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_near,'_'),function(x) x[1]))

dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD3/'
file_RCTD<-list.files(pattern = 'Deconvolution.txt',path = dir_RCTD,recursive = T)
table(unlist(lapply(strsplit(file_RCTD,'_'),function(x) x[1]))==dataset_slice)
cancer<-unlist(lapply(strsplit(file_near,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

all_FC<-c()
all_P<-c()
for(i in 1:length(cancer)){#i=1
  file_near_cancer<-file_near[grep(cancer[i],file_near)]
  file_RCTD_cancer<-file_RCTD[grep(cancer[i],file_RCTD)]
  dataset_slice_cancer<-unlist(lapply(strsplit(file_near_cancer,'_'),function(x) x[1]))
  
  cancer_RCTD<-c()
  for(j in 1:length(file_near_cancer)){#j=1
    st_near<-read.delim(paste0(dir_bdy,file_near_cancer[j]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD_cancer[j]),stringsAsFactors = F,check.names = F)
    
    spot_near<-unlist(strsplit(st_near[st_near$FinalLocalType%in%c('Boundary','Dispersion'),'step_1'],',')) %>% unique()
    spot_near<-st_near[spot_near,]
    spot_near<-spot_near$cell_name[spot_near$FinalLocalType%in%c('Immune','Normal')]
    
    bdy_dis_RCTD<-st_RCTD[spot_near,]
    cancer_RCTD<-rbind(cancer_RCTD,bdy_dis_RCTD)
  }
  cancer_RCTD<-cancer_RCTD[,setdiff(colnames(cancer_RCTD),c('Endothelial','Epithelial'))]
  CAF_TAM<-cancer_RCTD$CAF+cancer_RCTD$TAM
  if(length(intersect(c('B lymphocytes','T lymphocytes'),colnames(cancer_RCTD)))==2){
    Imm_other<-apply(cancer_RCTD[,intersect(c('B lymphocytes','T lymphocytes'),colnames(cancer_RCTD))],1,sum)
  }else{
    Imm_other<-cancer_RCTD[,intersect(c('B lymphocytes','T lymphocytes'),colnames(cancer_RCTD))]
  }
  
  
  FC<-mean(CAF_TAM)/mean(Imm_other)
  P_value<-wilcox.test(CAF_TAM,Imm_other,alternative ='greater')$p.value
  
  all_FC<-c(all_FC,FC)
  all_P<-c(all_P,P_value)
}

# ?wilcox.test
# P_value

plot_d<-data.frame(logFC=log2(all_FC),
                   logP=-log10(all_P),
                   P_value=all_P,
                   cancer=cancer,
                   x='logFC')
plot_d$logP[which(plot_d$logP>10)]<-10
plot_d$color<-'white'
plot_d$color[which(plot_d$P_value<0.05)]<-'#336633'####p值显著的画圈 也可以不画

p_dot<-ggplot(plot_d, aes(x=logFC, y=cancer,size=logFC,fill=logP)) +
  geom_point(shape = 21,  color = plot_d$color) + # 使用shape = 21画圈
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=2,col="black",lwd=0.6) +
  scale_fill_gradientn(colours = c(colorRampPalette(c("#899DA4","#FBEDD2"))(10),
                                   colorRampPalette(c("#FBEDD2","#C13710"))(90)) )+ #设置填充颜色
  scale_size_continuous(range = c(3, 10))+
  #geom_text(aes(label = num), vjust = 0.5,size = 5)+
  theme_minimal()
print(p_dot)
pdf('E:/Mirror/ST_analysis/pic/re/4/CAF_TAM_cluster/cancer_step1_celltypeCompare_T_B.pdf',height = 6,width = 6)
print(p_dot)
dev.off()