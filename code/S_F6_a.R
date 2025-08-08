####单细胞重新注释
###与单个数据注释的相一致的结果
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(circlize)



umap_meta<-read.csv('E:/Mirror/ST_analysis/data/SC_data/integrate/meta/umap_meta.csv',
                    stringsAsFactors = F,check.names = F,row.names = 1)
table(umap_meta$cancer)
umap_coor<-read.csv('E:/Mirror/ST_analysis/data/SC_data/integrate/meta/umap.csv',
                    stringsAsFactors = F,check.names = F)
umap_meta<-cbind(umap_meta,umap_coor)
colnames(umap_meta)
table(umap_meta$celltype)
umap_meta<-umap_meta[which(umap_meta$celltype!='unknown'),]

sc_cluster1<-read.delim('E:/Mirror/ST_analysis/other_data/SC_anno/sc_merge_medata2.txt',stringsAsFactors = F,check.names = F)
sc_cluster2<-read.delim('E:/Mirror/ST_analysis/other_data/SC_anno/sc_merge_medata3.txt',stringsAsFactors = F,check.names = F)
sc_cluster3<-read.delim('E:/Mirror/ST_analysis/other_data/SC_anno/sc_merge_medata4.txt',stringsAsFactors = F,check.names = F)
sc_cluster4<-read.delim('E:/Mirror/ST_analysis/other_data/SC_anno/sc_merge_medata5.txt',stringsAsFactors = F,check.names = F)

table(sc_cluster1[,1]==sc_cluster4[,1])
table(sc_cluster1$celltype_res1==sc_cluster2$celltype_res1)

sc_cluster<-data.frame(cell_name=sc_cluster1[,1],
                       cancer=sc_cluster1$cancer,
                       celltype=sc_cluster1$celltype_res1,
                       leiden_res_0.50=sc_cluster1$leiden_res_0.50,
                       leiden_res_0.80=sc_cluster1$leiden_res_0.80,
                       leiden_res_1.00=sc_cluster1$leiden_res_1.00,
                       leiden_res_1.20=sc_cluster1$leiden_res_1.20,
                       leiden_res_2.00=sc_cluster1$leiden_res_2.00,
                       leiden_res_3.00=sc_cluster1$leiden_res_3.00,
                       leiden_res_4.00=sc_cluster1$leiden_res_4.00,
                       leiden_res_5.00=sc_cluster1$leiden_res_5.00,
                       leiden_res_6.00=sc_cluster4$leiden_res_6.00,
                       leiden_res_7.00=sc_cluster4$leiden_res_7.00,
                       leiden_res_8.00=sc_cluster4$leiden_res_8.00,
                       leiden_res_9.00=sc_cluster2$leiden_res_9.00,
                       leiden_res_10.00=sc_cluster2$leiden_res_10.00,
                       leiden_res_11.00=sc_cluster2$leiden_res_11.00,
                       leiden_res_12.00=sc_cluster3$leiden_res_12.00,
                       leiden_res_13.00=sc_cluster3$leiden_res_13.00,
                       leiden_res_14.00=sc_cluster3$leiden_res_14.00,
                       leiden_res_15.00=sc_cluster3$leiden_res_15.00
)


dir_scAnno<-'E:/Mirror/ST_analysis/other_data/SC_anno/'
file_csAnno<-list.files(pattern = 'maxCelltype.txt',path = dir_scAnno)
res_n<-unlist(lapply(strsplit(file_csAnno,'_org_ratio'),function(x)x[1]))

res_n<-intersect(colnames(sc_cluster),res_n)
dir_pic<-'E:/Mirror/ST_analysis/pic/re/sc_integrate/match/'

#####定为res=14  i=17
sc_anno<-read.delim(paste0(dir_scAnno,file_csAnno[grep(res_n[i],file_csAnno)]),stringsAsFactors = F,check.names = F)
#sc_anno<-sc_anno[which(sc_anno$celltype!='-'),]

res_cluster<-sc_cluster[,c('cell_name','cancer',res_n[i])]
colnames(res_cluster)[3]<-'cluster'
res_cluster<-merge(res_cluster,sc_anno,by='cluster')
rownames(res_cluster)<-res_cluster$cell_name
res_cluster<-res_cluster[rownames(umap_meta),]
res_cluster$Umap1<-umap_meta$Umap1
res_cluster$Umap2<-umap_meta$Umap2
res_cluster$cellSubType<-umap_meta$celltype
table(res_cluster$celltype)
res_cluster<-res_cluster[which(res_cluster$celltype!='-'),]
res_cluster$cellNew<-res_cluster$celltype
res_cluster$cellNew[which(res_cluster$celltype=='Fibroblasts'&res_cluster$cellSubType=='CAF')]<-'CAF'
res_cluster$cellNew[which(res_cluster$celltype=='Macrophage'&res_cluster$cellSubType=='TAM')]<-'TAM'
table(res_cluster$cellNew)


res_cluster$Major<-res_cluster$cellSubType
res_cluster$Major[which(res_cluster$cellSubType=='CAF')]<-'Fibroblasts'
res_cluster$Major[which(res_cluster$cellSubType=='TAM')]<-'Macrophage'
####条形图#####################################
major_data<-as.data.frame(table(res_cluster$cancer,res_cluster$Major))
major_data<-major_data[!major_data$Var2%in%c('Endothelial','Epithelial'),]
major_data$Var2<-factor(major_data$Var2,levels = c(setdiff(unique(major_data$Var2),c('Fibroblasts','Macrophage')),
                                                   c('Fibroblasts','Macrophage')))
p_compare<-ggplot(major_data,aes(x=Var1,y=Freq,fill=Var2)) +
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
  ggtitle('SC_MajorCelltype')
print(p_compare)
pdf(paste0(dir_pic,'sc_MajorCelltype_bar.pdf'),height = 7,width = 7)
pdf(paste0(dir_pic,'sc_MajorCelltype_bar_outEN_EP.pdf'),height = 7,width = 7)
print(p_compare)
dev.off()











