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

p1<-ggplot(data = res_cluster,aes(x = Umap1, y = Umap2,color=cellNew)) + 
  geom_point(size=0.1)+
  scale_color_manual(values = c("B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                                "Macrophage"="#CC3333","Myeloid cell"="#ED703F","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                                'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                                "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'unknown'='grey85',
                                
                                "GC B cells in the DZ"='#CC9933',
                                "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#CCCC99',
                                "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
                                "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                                "Naive"='#FF9900',"B cell Regulatory"='#990033',"Naive B cell"='#990066',
                                'Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d'))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggtitle(res_n[i])
pdf(paste0(dir_pic,'res14_subType.pdf'),width = 6.5,height = 5)
print(p1)
dev.off()

