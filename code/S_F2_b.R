#####免疫受体表达
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(psych)


imm_check<-read.delim('E:/Mirror/ST_analysis/other_data/imm-check-gene.txt',stringsAsFactors = F,check.names = F)


dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)


for(i in 1:229){
  #i=1
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[imm_check$ligand,]
  
  imm_score<-data.frame(cell_name=st_bdy$cell_name,
                        LocalType=st_bdy$FinalLocalType,
                        dataSlice=dataset_slice[i],
                        Imm_score=apply(st_conut,2,function(x) mean(x,na.rm = T)))
  write.table(imm_score,paste0(dir_bdy,dataset_slice[i],'_ImmLigandScore.txt'),quote = F,sep = '\t')
  print(dataset_slice[i])
}

file_Imm<-list.files(pattern = '_ImmLigandScore.txt',path = dir_bdy,recursive = T)
all_slice_ImmScore<-c()
for(i in 1:229){
  imm_score<-read.delim(paste0(dir_bdy,file_Imm[i]),stringsAsFactors = F,check.names = F)
  all_slice_ImmScore<-rbind(all_slice_ImmScore,imm_score)
}
all_slice_ImmScore$cancer<-unlist(lapply(strsplit(all_slice_ImmScore$dataSlice,'/'),function(x) x[1]))
all_slice_ImmScore$cancer<-substr(all_slice_ImmScore$cancer,1,(nchar(all_slice_ImmScore$cancer)-2))%>%toupper()
write.table(all_slice_ImmScore,'E:/Mirror/ST_analysis/data/pic_data/all_slice_ImmLigandScore.txt',quote = F,sep = '\t')


all_slice_ImmScore<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_slice_ImmLigandScore.txt',stringsAsFactors = F,check.names = F)

quantile_mean<-function(x){#x=cancer_data$Imm_score[which(cancer_data$LocalType=='Core')]
  quantile_a<-quantile(x)
  mean_robust<-(0.5*quantile_a[3]+0.25*(quantile_a[2]+quantile_a[4]))
}

cancer<-unique(all_slice_ImmScore$cancer)
Imm_FC_CoreBdy<-c()
Imm_P_CoreBdy<-c()
Imm_FC_CoreDis<-c()
Imm_P_CoreDis<-c()
Imm_FC_BdyDis<-c()
Imm_P_BdyDis<-c()
Imm_FC_CoreVsBdyDis<-c()
Imm_P_CoreVsBdyDis<-c()


#quantile_a<-quantile(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')])
for(i in 1:length(cancer)){
  #i=1
  cancer_data<-all_slice_ImmScore[which(all_slice_ImmScore$cancer==cancer[i]),]
  table(cancer_data$LocalType)
  
  Imm_FC_CoreBdy<-c(Imm_FC_CoreBdy,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')])/
                      quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')]))
  Imm_P_CoreBdy<-c(Imm_P_CoreBdy,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')],
                                             cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')],
                                             alternative = 'greater')[["p.value"]])
  
  Imm_FC_CoreDis<-c(Imm_FC_CoreDis,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')])/
                      quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Dispersion')]))
  Imm_P_CoreDis<-c(Imm_P_CoreDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')],
                                             cancer_data$Imm_score[which(cancer_data$LocalType=='Dispersion')],
                                             alternative = 'greater')[["p.value"]])
  
  Imm_FC_BdyDis<-c(Imm_FC_BdyDis,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')])/
                     quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Dispersion')]))
  Imm_P_BdyDis<-c(Imm_P_BdyDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')],
                                           cancer_data$Imm_score[which(cancer_data$LocalType=='Dispersion')],
                                           alternative = 'greater')[["p.value"]])
  
  Imm_FC_CoreVsBdyDis<-c(Imm_FC_CoreVsBdyDis,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')])/
                           quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary'|cancer_data$LocalType=='Dispersion')]))
  Imm_P_CoreVsBdyDis<-c(Imm_P_CoreVsBdyDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')],
                                                       cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary'|cancer_data$LocalType=='Dispersion')],
                                                       alternative = 'greater')[["p.value"]])
}


compare_Imm_FC<-data.frame(CoreVsBdyDis=Imm_FC_CoreVsBdyDis,
                           CoreVsBdy=Imm_FC_CoreBdy,
                           CoreVsDis=Imm_FC_CoreDis,
                           BdyVsDis=Imm_FC_BdyDis,
                           row.names = cancer)
compare_Imm_P<-data.frame(CoreVsBdyDis=Imm_P_CoreVsBdyDis,
                          CoreVsBdy=Imm_P_CoreBdy,
                          CoreVsDis=Imm_P_CoreDis,
                          BdyVsDis=Imm_P_BdyDis,
                          row.names = cancer)

###统计每种癌症种各种大于的基因个数
num_CoreBdy<-c()
num_CoreDis<-c()
num_BdyDis<-c()
num_CoreBdyDis<-c()
all_immLigand<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_immLigand_CoreBdyDis.txt',stringsAsFactors = F,check.names = F)
for(i in 1:length(cancer)){
  #i=1
  cancer_data2<-all_immLigand[which(all_immLigand$cancer==cancer[i]),]
  
  j=0
  k=0
  l=0
  m=0
  for(nn in 5:17){##nn=5
    gene_data<-cancer_data2[,c(2,nn)]
    Core<-mean(gene_data[gene_data$LocalType%in%c('Core'),2])
    Bdy<-mean(gene_data[gene_data$LocalType%in%c('Boundary'),2])
    BdyDis<-mean(gene_data[gene_data$LocalType%in%c('Boundary','Dispersion'),2])
    Dis<-mean(gene_data[gene_data$LocalType%in%c('Dispersion'),2])
    #print(c(Core,Bdy,BdyDis,Dis))
    if(Core>Bdy) j=j+1
    if(Core>Dis) k=k+1
    if(Bdy>Dis) l=l+1
    if(Core>BdyDis) m=m+1
  }
  num_CoreBdy<-c(num_CoreBdy,j)
  num_CoreDis<-c(num_CoreDis,k)
  num_BdyDis<-c(num_BdyDis,l)
  num_CoreBdyDis<-c(num_CoreBdyDis,m)
  
}
compare_Imm_num<-data.frame(CoreVsBdyDis=num_CoreBdyDis,
                            CoreVsBdy=num_CoreBdy,
                            CoreVsDis=num_CoreDis,
                            BdyVsDis=num_BdyDis,
                            row.names = cancer)


compare_Imm_FC<-reshape2::melt(as.matrix(compare_Imm_FC))
compare_Imm_P<-reshape2::melt(as.matrix(compare_Imm_P))
compare_Imm_num<-reshape2::melt(as.matrix(compare_Imm_num))
compare_Imm<-compare_Imm_FC
colnames(compare_Imm)<-c('cancer','versus','FC')
compare_Imm$P<-compare_Imm_P$value
compare_Imm$logP<-(-log10(compare_Imm$P))
compare_Imm$logP[which(compare_Imm$logP>10)]<-10
compare_Imm$num<-compare_Imm_num$value
compare_Imm$color<-'white'
compare_Imm$color[which(compare_Imm$P<0.05)]<-'#336633'####p值显著的画圈 也可以不画
#compare_Imm<-compare_Imm[order(compare_Imm$cancer,decreasing = F),]
compare_Imm<-mutate(compare_Imm,cancer = factor(compare_Imm$cancer, levels = rev(c('PCNSL','GIST','CSCC','HGSC','MIBC','OSCC','GBM',
                                                                                   'RCC','HN-AS','IPMN','LUAD','BRCA','LIHC','CESC',
                                                                                   'CRC','SKCM','OVCA','PDAC','PRAD'))))
unique(compare_Imm$cancer)


p_dot<-ggplot(compare_Imm, aes(x=versus, y=cancer,size=FC,fill=logP)) +
  geom_point(shape = 21,  color = compare_Imm$color) + # 使用shape = 21画圈
  scale_fill_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(20),
                                   colorRampPalette(c("#F39C67","#B20A1C"))(80)) )+ #设置填充颜色
  scale_size_continuous(range = c(2, 8))+
  geom_text(aes(label = num), vjust = 0.5,size = 5)+
  theme_minimal()
print(p_dot)
pdf('E:/Mirror/ST_analysis/pic/re/3/ImmLigand_compare_new_mean2.pdf',width = 5,height = 6)
print(p_dot)
dev.off()

?geom_text



