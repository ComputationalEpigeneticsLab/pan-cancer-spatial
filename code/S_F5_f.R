####各种基因集得分 对MP聚类的5个簇可视化热图
##使用恶性spot的均值，再对每个MP类的所有切片算均值
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gghalves)
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)


###转移相关基因集
geneset1<-read.delim('E:/Mirror/ST_analysis/data/geneset/genest_mp.txt',stringsAsFactors = F,check.names = F)
all_geneset<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(all_geneset)<-geneset1$pathway
####可以用AddModuleScore计算的
all_geneset_use<-all_geneset

dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dir_EMT<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'

for(i in 1:length(file_bdy)){
  #i=2
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  
  st_rds<-st_rds[,rownames(st_bdy)]
  
  all_geneset_use<-all_geneset
  inter_gene<-lapply(all_geneset,function(x) length(intersect(x,rownames(st_rds))))%>%unlist()
  if(length(which(inter_gene<=0))>0){
    for(j in 1:length(which(inter_gene<=0))){
      all_geneset_use[[which(inter_gene<=0)[j]]]<-rownames(st_rds)[1:3]
    }
  }
  
  st_rds<-AddModuleScore(st_rds,all_geneset_use,name = names(all_geneset_use))
  score_re<-st_rds@meta.data
  score_re<-score_re[,paste0(names(all_geneset_use),1:length(names(all_geneset_use)))]
  colnames(score_re)<-names(all_geneset_use)
  score_re[,which(inter_gene<=0)]<-NA######没有交集的基因集得分为NA
  
  score_re$LocalType<-st_bdy$FinalLocalType
  score_re$imagerow<-st_bdy$imagerow
  score_re$imagecol<-st_bdy$imagecol
  score_re$cell_name<-rownames(score_re)
  
  write.table(score_re,paste0(dir_EMT,dataset_slice[i],'_MP5_score.txt'),quote = F,sep = '\t')
  print(dataset_slice[i])
}


dir_EMT<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'
file_metas<-list.files(pattern = '_MP5_score.txt',path = dir_EMT,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_metas,'_'),function(x) x[1]))


slice_MP5<-data.frame(pathway=names(all_geneset))
for(i in 1:length(file_metas)){
  #i=1
  metas_data<-read.delim(paste0(dir_EMT,file_metas[i]),stringsAsFactors = F,check.names = F)
  #mal_spot<-metas_data$cell_name[metas_data$LocalType%in%c('Core','Boundary','Dispersion')]
  mal_spot<-metas_data$cell_name[metas_data$LocalType%in%c('Core')]
  
  mal_metas<-metas_data[mal_spot,]
  slice_MP5<-cbind(slice_MP5,data.frame(score=apply(mal_metas[1:18],2,mean)))
  colnames(slice_MP5)[i+1]<-dataset_slice[i]
}
slice_MP5<-slice_MP5[,-1]
slice_MP5<-t(slice_MP5) %>% as.data.frame()

cluster_MP<-read.csv('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/pam_pearson/pam_pearson.k=5.consensusClass.csv',
                     header = F,stringsAsFactors = F,check.names = F)
cluster_MP$V1<-gsub('hn.as','hn-as',cluster_MP$V1)
cluster_MP$V1<-gsub('\\.','/',cluster_MP$V1)
table(cluster_MP$V2)

slice_MP5<-slice_MP5[cluster_MP$V1,]
slice_MP5<-cbind(cluster_MP,slice_MP5)
slice_MP5$V2<-paste0('C',slice_MP5$V2)
#slice_MP5<-slice_MP5[,-c(3,4)]
write.table(slice_MP5,'E:/Mirror/ST_analysis/data/pic_data/MP_cluster_genesetScore_core.txt',quote = F,sep = '\t')
slice_MP5<-read.delim('E:/Mirror/ST_analysis/data/pic_data/MP_cluster_genesetScore_core.txt',stringsAsFactors = F,check.names = F)
table(slice_MP5$V2)

quantile_mean<-function(x){#x=cancer_data$Imm_score[which(cancer_data$LocalType=='Core')]
  quantile_a<-quantile(x,na.rm = T)
  mean_robust<-(0.5*quantile_a[3]+0.25*(quantile_a[2]+quantile_a[4]))
}

library(base)
MP5_scoreMean<-aggregate(slice_MP5[,3:20],by=list(slice_MP5$V2),function(x){mean(x,na.rm=T)})
# mean(c(1,2,3,NA,NA),na.rm=T)
# mean(c(1,2,3,NA,NA),na.rm=F)
rownames(MP5_scoreMean)<-MP5_scoreMean$Group.1
MP5_scoreMean<-MP5_scoreMean[,-1]
MP5_scoreMean<-MP5_scoreMean[,-c(2,14,17,18)]
MP5_scoreMean<-MP5_scoreMean[,c(1:3,8,4:7,9:ncol(MP5_scoreMean))]
#apply(MP5_scoreMean,2,function(x){wilcox.test()})
kw_pValue<-c()
for(i in 3:20){
  fit <- kruskal.test(slice_MP5[,i] ~ slice_MP5[,2], data = slice_MP5)
  kw_pValue<-c(kw_pValue,fit[["p.value"]])
}
kw_pValue<-data.frame(p_value=kw_pValue,pathway=colnames(slice_MP5)[3:20])


heat_data<-MP5_scoreMean
range(heat_data)
heat_data<-t(scale(heat_data))###按列即按通路归一化
#heat_data<-t(scale(t(heat_data)))###按行即按MP类别归一化
range(heat_data)
bk<-seq(-1,1,length.out=100)
color_pheatmap<-c(colorRampPalette(c("#0669AD",'#89BDD9'))(20),
                  colorRampPalette(c("#89BDD9",'white'))(30),
                  colorRampPalette(c("white",'#E9C1C6'))(10),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(40)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(heat_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      #show_colnames = F,
                      cluster_rows = T,
                      cluster_cols = F,
                      treeheight_row = F,
                      # treeheight_col = T,
                      #display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "black",
                      fontsize = 10,
                      # cellwidth=15,
                      # cellheight=15,
                      main = "MP_cluster_pathwayScore",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/NMF_module/MP_cluster_pathwayScore_core2.pdf',width = 5,height = 6)
print(p)
dev.off()







