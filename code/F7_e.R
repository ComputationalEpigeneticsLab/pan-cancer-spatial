####基因集计算得分
####对每种癌症 使用所有切片core区marker的RRA排序结果和各种基因集计算GSEA
library(org.Hs.eg.db) #人类注释数据
# BiocManager::install('GO.db')
# BiocManager::install('clusterProfiler')
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(dplyr)
library(enrichplot)
library(fgsea)
library(stats)
library(Seurat)
library(stringr)
library(plyr)
#library(estimate)
library(RobustRankAggreg)


###干性基因集
geneset1<-read.delim('E:/Mirror/ST_analysis/data/geneset/stem_geneset.txt',stringsAsFactors = F,check.names = F)
stemness_geneset<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(stemness_geneset)<-paste0('Stem_of_',geneset1$geneset_name)

####可以用GSEA计算的
all_geneset<-stemness_geneset
pathway_data<-data.frame(pathway=rep(names(all_geneset),unlist(lapply(all_geneset,length))),
                         gene=unlist(all_geneset))

####FC数据
dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))


dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
table(unlist(lapply(strsplit(file_bdy,'_'),function(x) x[1]))==dataset_slice)

GSEA_ESscore<-data.frame(pathway=names(all_geneset))
GSEA_FDR<-data.frame(pathway=names(all_geneset))

for(i in 1:4){
  #i=1
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_count<-as.matrix(st_rds@assays[["Spatial"]]@counts)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
  st_count<-st_count[,st_bdy$cell.names]
  core_site<-which(st_bdy$FinalLocalType=='Core')
  other_site<-which(st_bdy$FinalLocalType!='Core')
  FC_slice<-apply(st_count,1,function(x){#x=st_count[1,]
    aa<-mean(x[core_site])-mean(x[other_site])
  })
  FC_slice<-FC_slice[order(FC_slice,decreasing = T)]
  
  
  GSEA_res <- GSEA(FC_slice,TERM2GENE = pathway_data,minGSSize = 2, maxGSSize = 1000, pvalueCutoff=1)
  
  slice_ESscore<-data.frame(pathway=GSEA_res@result$Description,
                            ESscore=GSEA_res@result$enrichmentScore)#enrichmentScore
  colnames(slice_ESscore)[2]<-toupper(dataSet[i])
  slice_FDR<-data.frame(pathway=GSEA_res@result$Description,
                        FDR=GSEA_res@result$p.adjust)
  colnames(slice_FDR)[2]<-toupper(dataSet[i])
  
  GSEA_ESscore<-merge(GSEA_ESscore,slice_ESscore,by='pathway',all=T)
  GSEA_FDR<-merge(GSEA_FDR,slice_FDR,by='pathway',all=T)
  print(dataSet[i])
}
rownames(GSEA_FDR)<-GSEA_FDR$pathway
GSEA_FDR<-GSEA_FDR[,-1]
rownames(GSEA_ESscore)<-GSEA_ESscore$pathway
GSEA_ESscore<-GSEA_ESscore[,-1]
#pathway_use<-names(which(apply(GSEA_ESscore,1,function(x){length(which(is.na(x)==T))})<19))

# GSEA_FDR<-GSEA_FDR[pathway_use,]
# GSEA_ESscore<-GSEA_ESscore[pathway_use,]

write.table(GSEA_ESscore,'E:/Mirror/ST_analysis/data/ESCC/pic_data/GSEA_Stem&genesetScore_ESscore.txt',quote = F,sep = '\t')
write.table(GSEA_FDR,'E:/Mirror/ST_analysis/data/ESCC/pic_data/GSEA_Stem&genesetScore_FDR.txt',quote = F,sep = '\t')


GSEA_ESscore<-read.delim('E:/Mirror/ST_analysis/data/ESCC/pic_data/GSEA_Stem&genesetScore_ESscore.txt',
                         stringsAsFactors = F,check.names = F)
GSEA_FDR<-read.delim('E:/Mirror/ST_analysis/data/ESCC/pic_data/GSEA_Stem&genesetScore_FDR.txt',
                     stringsAsFactors = F,check.names = F)



##############################################################################################
####为每个spot计算
dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))
dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
table(unlist(lapply(strsplit(file_bdy,'_'),function(x) x[1]))==dataset_slice)


for(i in 1:4){
  #i=1
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  
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
  
  write.table(score_re,paste0(dir_bdy,dataset_slice[i],'_genesetScore.txt'),quote = F,sep = '\t')
  print(dataset_slice[i])
  
}



dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
file_stem<-list.files(pattern = '_genesetScore.txt',path = dir_bdy,recursive = T)
patient<-unlist(lapply(strsplit(file_bdy,'/'),function(x)x[1]))

all_stem<-c()
for(i in 1:2){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  colnames(st_bdy)
  st_stem<-read.delim(paste0(dir_bdy,file_stem[i]),stringsAsFactors = F,check.names = F)
  colnames(st_stem)<-unlist(lapply(strsplit(colnames(st_stem),'_of_'),function(x)x[2]))
  st_stem<-reshape2::melt(as.matrix(st_stem))
  colnames(st_stem)<-c('cell.names','stem','score')
  st_stem<-merge(st_stem,st_bdy[,c('cell.names','FinalLocalType')],by='cell.names',all=T)
  st_stem$patient<-patient[i]
  all_stem<-rbind(all_stem,st_stem)
}


pdf('E:/Mirror/ST_analysis/pic/ESCC/StemAddModuleScore_compare.pdf',width = 8,height = 4)
P_stem<-all_stem[all_stem$FinalLocalType%in%c('Core','Boundary','Dispersion'),]
p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=patient))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = patient), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
  scale_color_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
  ggtitle('all_mal')+
  #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =0.6)+
  stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
  #theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
        axis.title.x = element_text(size=12),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))
print(p_box)
for(i in c('Core','Boundary','Dispersion')){#i='Core'
  P_stem<-all_stem[all_stem$FinalLocalType%in%i,]
  #P_stem<-P_stem[P_stem$stem%in%c(),]
  p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=patient))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = patient), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
    scale_color_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
    ggtitle(i)+
    #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =0.6)+
    stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
    #theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size=10),
          axis.title.x = element_text(size=12),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10),
          panel.background = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(size = 12,colour = "black"),
          axis.title = element_text(size = 15))
  print(p_box)
}
dev.off()









