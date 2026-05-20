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

dir_sc<-'/data/SC_data/re/'
file_sc<-list.files(pattern = '_meta.txt',path = dir_sc,recursive = T)
file_sc
cancer<-unlist(lapply(strsplit(file_sc,'/'),function(x)x[1]))
cancer

res_cluster<-c()
for(i in 1:length(file_sc)){
  sc_meta<-read.delim(paste0(dir_sc,file_sc[i]),stringsAsFactors = F,check.names = F)
  sc_meta$subCAFTAM[grep('CAF',sc_meta$subCAFTAM)]<-'CAF'
  sc_meta$subCAFTAM[grep('TAM',sc_meta$subCAFTAM)]<-'TAM'
  
  aa<-data.frame(cancer=cancer[i],
                 cellSubType=sc_meta$subCAFTAM)
  res_cluster<-rbind(res_cluster,aa)
  
}
write.table(res_cluster,'/data/SC_data/re/all_subType.txt',
            quote = F,sep = '\t',row.names = F)

dir_pic<-'/Fs6/'
res_cluster<-read.delim('/Fs6/FigS6ab/all_subType.txt',
                        stringsAsFactors = F,check.names = F)
res_cluster<-res_cluster[which(res_cluster$cellSubType!='unknown'),]
sub_data<-as.data.frame(table(res_cluster$cancer,res_cluster$cellSubType))
# sub_data<-sub_data[!sub_data$Var2%in%c('Endothelial','Epithelial'),]
sub_data$Var2<-factor(sub_data$Var2,levels = c(setdiff(unique(sub_data$Var2),c('CAF','TAM')),
                                               c('CAF','TAM')))
unique(sub_data$Var2)
p_compare<-ggplot(sub_data,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  coord_flip()+
  scale_fill_manual(values = c("B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                               "Macrophage"="#CC3333","Myeloid cell"="#ED703F","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                               'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                               "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",
                               
                               "GC B cells in the DZ"='#CC9933',
                               "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',
                               "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
                               "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                               "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066' ))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 12,colour = "black"),
    axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('SC_SubCelltype')
print(p_compare)
pdf(paste0(dir_pic,'sc_SubCelltype_bar.pdf'),height = 7,width = 7)
pdf(paste0(dir_pic,'sc_SubCelltype_bar_outEN_EP.pdf'),height = 9,width = 7)
print(p_compare)
dev.off()



