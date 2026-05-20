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


dir_out<-'/Fs3/'

imm_check<-read.delim('/imm-check-gene.txt',stringsAsFactors = F,check.names = F)

all_slice_ImmScore<-read.delim(paste0(dir_out,'all_slice_ImmLigandScore.txt'),stringsAsFactors = F,check.names = F)

unique(all_slice_ImmScore$dataSlice)

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
                      quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Budding')]))
  Imm_P_CoreDis<-c(Imm_P_CoreDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')],
                                             cancer_data$Imm_score[which(cancer_data$LocalType=='Budding')],
                                             alternative = 'greater')[["p.value"]])
  
  Imm_FC_BdyDis<-c(Imm_FC_BdyDis,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')])/
                     quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Budding')]))
  Imm_P_BdyDis<-c(Imm_P_BdyDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary')],
                                           cancer_data$Imm_score[which(cancer_data$LocalType=='Budding')],
                                           alternative = 'greater')[["p.value"]])
  
  Imm_FC_CoreVsBdyDis<-c(Imm_FC_CoreVsBdyDis,quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')])/
                           quantile_mean(cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary'|cancer_data$LocalType=='Dispersion')]))
  Imm_P_CoreVsBdyDis<-c(Imm_P_CoreVsBdyDis,wilcox.test(cancer_data$Imm_score[which(cancer_data$LocalType=='Core')],
                                                       cancer_data$Imm_score[which(cancer_data$LocalType=='Boundary'|cancer_data$LocalType=='Dispersion')],
                                                       alternative = 'greater')[["p.value"]])
}


compare_Imm_FC<-data.frame(CoreVsBdyBud=Imm_FC_CoreVsBdyDis,
                           CoreVsBdy=Imm_FC_CoreBdy,
                           CoreVsBud=Imm_FC_CoreDis,
                           BdyVsBud=Imm_FC_BdyDis,
                           row.names = cancer)
compare_Imm_P<-data.frame(CoreVsBdyBud=Imm_P_CoreVsBdyDis,
                          CoreVsBdy=Imm_P_CoreBdy,
                          CoreVsBud=Imm_P_CoreDis,
                          BdyVsBud=Imm_P_BdyDis,
                          row.names = cancer)



compare_Imm_FC<-reshape2::melt(as.matrix(compare_Imm_FC))
compare_Imm_P<-reshape2::melt(as.matrix(compare_Imm_P))
compare_Imm<-compare_Imm_FC
colnames(compare_Imm)<-c('cancer','versus','FC')
compare_Imm$P<-compare_Imm_P$value
compare_Imm$logP<-(-log10(compare_Imm$P))
compare_Imm$logP[which(compare_Imm$logP>10)]<-10
compare_Imm$color<-'white'
compare_Imm$color[which(compare_Imm$P<0.05)]<-'#336633'
#compare_Imm<-compare_Imm[order(compare_Imm$cancer,decreasing = F),]
compare_Imm$cancer<-factor(compare_Imm$cancer, levels = rev(c('CRC','HNSCC','IPMN','CESC','CSCC','LIHC','DSRCT','EC','LUAD','LUSC','MIBC','GIST',
                                                               'HB','HGSC','PCNSL','PDAC','PN','PRAD','RCC','TC','OS','OSCC','OVCA','GBM',
                                                               'LGACC','LNC','HN-AS','NPC','GC','BRCA','DIPG','SKCM','LAM','PTCL','DLBCL')))
#unique(compare_Imm$cancer)


p_dot<-ggplot(compare_Imm, aes(x=versus, y=cancer,size=FC)) +
  geom_point(shape = 19,  aes(colour = logP)) + # 使用shape = 21画圈
  scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(20),
                                    colorRampPalette(c("#F39C67","#B20A1C"))(80)) )+ #设置填充颜色
  scale_size_continuous(range = c(2, 8))+
  theme_minimal()+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
    axis.text = element_text(size = 10,colour = "black"),
    strip.background = element_blank(),
    strip.text.x = element_text(size=9),
    strip.text.y = element_text(size=9, face="bold"))
print(p_dot)


pdf('/ImmLigand_compare_new_mean2.pdf',width = 5,height = 8.5)
print(p_dot)
dev.off()

