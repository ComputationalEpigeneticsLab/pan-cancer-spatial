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
library(ggpubr)

###干性基因集
geneset1<-read.delim('/data/10X_Visium/new_st_StemGSEA/stem_geneset.txt',stringsAsFactors = F,check.names = F)
stemness_geneset<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(stemness_geneset)<-paste0('Stem_of_',geneset1$geneset_name)

####可以用GSEA计算的
all_geneset<-stemness_geneset
pathway_data<-data.frame(pathway=rep(names(all_geneset),unlist(lapply(all_geneset,length))),
                         gene=unlist(all_geneset))


##############################################################################################
####为每个spot计算
dir_rds<-'/data/10X_Visium/new_st/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)

dir_bdy<-'/data/10X_Visium/new_st_copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)

g_r_noTIB=c("GBM_GSE235672_GSM7507327","GBM_GSE235672_GSM7507330",
            "LIHC_GSE238264_GSM7661255","LIHC_GSE238264_GSM7661256","LIHC_GSE238264_GSM7661257","LIHC_GSE238264_GSM7661258",'LIHC_lihc03_slice7' )

g_nr_TIB=c("GBM_GSE235672_GSM7507312","GBM_GSE235672_GSM7507328","LIHC_GSE238264_GSM7661260",'LIHC_lihc03_slice5')

g_nr_noTIB=c("GBM_GSE235672_GSM7507323","GBM_GSE235672_GSM7507329","LIHC_GSE238264_GSM7661259","LIHC_GSE238264_GSM7661261")

g_nr <-c(g_nr_TIB,g_nr_noTIB)
g <- c(g_nr,g_r_noTIB)


file_bdy_name <- gsub("/", "_", file_bdy)
file_bdy_name <- gsub("_BdyTumorCore.txt", "", file_bdy_name)

g_nr_positions <- which(file_bdy_name %in% g_nr)
g_r_positions <- which(file_bdy_name %in% g_r_noTIB)
g_nr_TIB_positions <- which(file_bdy_name %in% g_nr_TIB)
g_nr_noTIB_positions <- which(file_bdy_name %in% g_nr_noTIB)
g_positions <- which(file_bdy_name %in% g)

dir_out <- "/data/10X_Visium/plot/stemnessScore/"

for(i in g_positions){
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
  
  slice_name <- file_bdy_name[i]
  
  write.table(score_re,paste0(dir_out,slice_name,'_genesetScore.txt'),quote = F,sep = '\t')
  print(slice_name)
  
}


file_stem<-list.files(pattern = '_genesetScore.txt',path = dir_out,recursive = T)


all_stem<-c()
for(i in g_positions){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  colnames(st_bdy)
  slice_name <- file_bdy_name[i]
  
  st_stem<-read.delim(paste0(dir_out,slice_name,'_genesetScore.txt'),stringsAsFactors = F,check.names = F)
  colnames(st_stem)<-unlist(lapply(strsplit(colnames(st_stem),'_of_'),function(x)x[2]))
  st_stem<-reshape2::melt(as.matrix(st_stem))
  colnames(st_stem)<-c('cell_name','stem','score')
  st_stem<-merge(st_stem,st_bdy[,c('cell_name','FinalLocalType')],by='cell_name',all=T)
  st_stem$patient<-slice_name
  all_stem<-rbind(all_stem,st_stem)
  
}
write.table(all_stem,paste0(dir_out,'allslice_genesetScore.txt'),quote = F,sep = '\t')

all_stem <- read.table("/allslice_genesetScore.txt",header = T,sep = "\t")

all_stem$group <- "g"
all_stem[all_stem$patient %in% g_r_noTIB,]$group <- "g_r_noTIB"
all_stem[all_stem$patient %in% g_nr_TIB,]$group <- "g_nr_TIB"
all_stem[all_stem$patient %in% g_nr_noTIB,]$group <- "g_nr_noTIB"
all_stem$group2 <- all_stem$group 
all_stem[all_stem$group == "g_nr_TIB",]$group2 <- "g_nr"
all_stem[all_stem$group == "g_nr_noTIB",]$group2 <- "g_nr"

pdf('/StemAddModuleScore_compare.pdf',width = 8,height = 4)
P_stem<-all_stem[all_stem$FinalLocalType%in%c('Core','Boundary','Dispersion'),]
p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=group2))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = group2), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("g_nr"="#15629e","g_r_noTIB"="#ddb424"))+
  scale_color_manual(values = c("g_nr"="#15629e","g_r_noTIB"="#ddb424"))+
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

for(i in c('Core','Boundary','Budding')){#i='Core'
  P_stem<-all_stem[all_stem$FinalLocalType%in%i,]
  p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=group2))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = group2), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("g_nr"="#15629e","g_r_noTIB"="#ddb424"))+
    scale_color_manual(values = c("g_nr"="#15629e","g_r_noTIB"="#ddb424"))+
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


#nrtib_vs_r
pdf('7/StemAddModuleScore_compare_nrtib_vs_r.pdf',width = 8,height = 4)
P_stem<-all_stem[all_stem$FinalLocalType%in%c('Core','Boundary','Dispersion'),]
P_stem <- P_stem[P_stem$group%in%c("g_r_noTIB","g_nr_TIB"),]

p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=group))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = group), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("g_nr_TIB"="#15629e","g_r_noTIB"="#ddb424"))+
  scale_color_manual(values = c("g_nr_TIB"="#15629e","g_r_noTIB"="#ddb424"))+
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

for(i in c('Core','Boundary','Budding')){#i='Core'
  P_stem<-all_stem[all_stem$FinalLocalType%in%i,]
  P_stem <- P_stem[P_stem$group%in%c("g_r_noTIB","g_nr_TIB"),]
  p_box<-ggplot(P_stem, aes(x = stem, y = score,fill=group))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.3,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = group), color='black',width = 0.6,lwd=0.1,#fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("g_nr_TIB"="#15629e","g_r_noTIB"="#ddb424"))+
    scale_color_manual(values = c("g_nr_TIB"="#15629e","g_r_noTIB"="#ddb424"))+
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

