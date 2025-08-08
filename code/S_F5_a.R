####MP14基因交集热图
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




Program_713_org<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/program_top50gene.txt',
                            stringsAsFactors = F,check.names = F)
Program_713_org$cluster<-gsub("Cluster","MP",Program_713_org$cluster)
table(Program_713_org$cluster)
Program_713<-Program_713_org

order_MP<-c('MP_3','MP_2','MP_5','MP_4','MP_7','MP_6','MP_9','MP_11','MP_10','MP_8','MP_13','MP_12','MP_14')

Program_713<-Program_713[lapply(order_MP,function(x) which(Program_713$cluster==x))%>%unlist(),]
#lapply(order_MP,function(x) which(Program_713$cluster==x))%>%unlist()

Program_713_gene<-lapply(Program_713$gene,function(x) unlist(strsplit(x,',')))
Program_713_gene<-do.call(cbind,Program_713_gene)
colnames(Program_713_gene)<-Program_713$slice

Program_713$cancer<-unlist(lapply(strsplit(Program_713$slice,'_'),function(x) x[1]))
Program_713$cancer<-substr(Program_713$cancer,1,(nchar(Program_713$cancer)-2))%>%toupper()
Program_713$data_slice<-unlist(lapply(strsplit(Program_713$slice,'.RDS'),function(x) x[1]))
aaa<-as.data.frame.array(table(Program_713$data_slice,Program_713$cluster))
program_cancer<-as.data.frame.array(table(Program_713$cancer,Program_713$cluster))
program_cancer<-reshape2::melt(as.matrix(program_cancer))
MP_order<-c('MP_1',
            'MP_2',
            'MP_3','MP_6','MP_7','MP_14',
            'MP_4', 
            'MP_5','MP_12', 'MP_13',  
            'MP_8', 'MP_9', 
            'MP_10',  
            'MP_11')#%>%toupper()
program_cancer<-mutate(program_cancer,Legend = factor(program_cancer$Var2, levels = as.character(MP_order)))
table(program_cancer$Var1,program_cancer$Var2)

p_compare<-ggplot(program_cancer,aes(x=Legend,y=value,fill=Var1)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  scale_fill_manual(values = c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                               'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                               'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                               'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b'))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("MP")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('713program_cancer')
print(p_compare)
pdf('E:/Mirror/ST_analysis/pic/re/5/program_cancer.pdf',height = 6,width = 8)
print(p_compare)
dev.off()

