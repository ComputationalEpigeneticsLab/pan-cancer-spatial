####对单细胞数据中的CAF细胞分亚型
###从文章中获取8个CAF亚型的marker
library(devtools)
#remove.packages(c("Seurat","SeuratObject"))
# Version("Matrix")
# packageDescription('Matrix')
# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
library(Seurat)
library(SeuratObject)
packageVersion("Seurat")

library(stringr)
library(dplyr)
library(coda)
library(rjags)
library(tidyverse)
library(infercnv)
library(tidyverse)
library(circlize)
library(hdf5r)
#BiocManager::install("rhdf5")
library(rhdf5)
library(Matrix)

source("E:/Mirror/ST_analysis/program/sctype_score.R")
marker<-read.csv('E:/Mirror/ST_analysis/other_data/subCAF_marker.csv',stringsAsFactors = F,check.names = F)
positive<-lapply(unique(marker$CAF_type),function(x) marker$marker[which(marker$CAF_type==x)])
names(positive)<-unique(marker$CAF_type)
negative<-lapply(positive, function(x) character())
marker_list<-list(positive=positive,
                  negative=negative)


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)
cancer<-unlist(lapply(strsplit(file_sc,'/'),function(x)x[1]))

dir_out<-'E:/Mirror/ST_analysis/data/SC_data/subCAF/'

for(i in 2:length(file_sc)){
  #i=2
  if(!dir.exists(paste0(dir_out,cancer[i]))) dir.create(paste0(dir_out,cancer[i]))
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  #intersect(rownames(sc_rds),'C15ORF48')
  sc_CAF<-subset(sc_rds,subset = seurat_clusters == "CAF")
  
  DataSeuratObject <- CreateSeuratObject(counts = sc_CAF@assays[["RNA"]]@counts, project = "cancer",min.cells = 3)
  DataSeuratObject <- NormalizeData(DataSeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  DataSeuratObject <- FindVariableFeatures(DataSeuratObject, selection.method = "vst", nfeatures = 2000)
  #Scaling the data
  all.genes <- rownames(DataSeuratObject)
  DataSeuratObject <- ScaleData(DataSeuratObject, features = all.genes)
  DataSeuratObject <- RunPCA(DataSeuratObject, features = VariableFeatures(object = DataSeuratObject))
  DataSeuratObject <- FindNeighbors(DataSeuratObject, dims = 1:25)
  DataSeuratObject <- FindClusters(DataSeuratObject, resolution = 1.5) 
  DataSeuratObject <- RunUMAP(DataSeuratObject, reduction = "pca", dims = 1:30)####顺便把单细胞的UMAP坐标输出
  
  #aa<-DataSeuratObject@assays[["RNA"]]@layers[["scale.data"]]
  ##########注意有几套单细胞数据######################################################################
  es.max = sctype_score(scRNAseqData = DataSeuratObject[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = marker_list$positive, gs2 = marker_list$negative)
  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(DataSeuratObject@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(DataSeuratObject@meta.data[DataSeuratObject@meta.data$seurat_clusters==cl, ])]), 
                     decreasing = !0)
    ##构建数据框
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, 
                    ncells = sum(DataSeuratObject@meta.data$seurat_clusters==cl),stringsAsFactors = F), 10)
  }))
  ####cluster re ####
  cluster_celltype = data.frame(cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores) ,stringsAsFactors = F )
  #cluster_celltype$type[as.numeric(as.character(cluster_celltype$scores)) < cluster_celltype$ncells/4] = "unknown"
  ####change seuratobject Idents & seurat_clusters####
  cluster_celltype_mtx <- as.matrix(cluster_celltype)
  #cluster_celltype_mtx[,2]=add_postfix(cluster_celltype_mtx[,2])
  rownames(cluster_celltype_mtx) <- cluster_celltype_mtx[,1]
  Idents_new=cluster_celltype_mtx[as.character(DataSeuratObject$seurat_clusters),2]
  names(Idents_new)=names(DataSeuratObject$seurat_clusters)
  Seurat::Idents(object = DataSeuratObject)=Idents_new
  DataSeuratObject$seurat_clusters=Idents_new
  ###提出细胞和细胞类型######
  cell_celltype<-data.frame(DataSeuratObject$seurat_clusters)
  cell_celltype$cellname=rownames(cell_celltype)
  colnames(cell_celltype)=c("celltype","cellname")
  #cell_celltype$source<-GSM[i]
  #cell_celltype$source<-as.vector(DataSeuratObject@meta.data[["orig.ident"]])
  table(cell_celltype$celltype)
  
  all_meta<-sc_rds@meta.data
  all_meta_CAF<-all_meta[setdiff(rownames(all_meta),rownames(cell_celltype)),]
  CAF_add<-data.frame(row.names = c(rownames(all_meta_CAF),rownames(cell_celltype)),
                      cellType=c(all_meta_CAF$seurat_clusters,cell_celltype$celltype),
                      lab='add')
  CAF_add<-CAF_add[colnames(sc_rds),]
  
  sc_rds<-AddMetaData(sc_rds,metadata = CAF_add$cellType,col.name = 'subCAF')
  table(sc_rds@meta.data[["seurat_clusters"]],sc_rds@meta.data[["subCAF"]])
  
  print(cancer[i])
  
  saveRDS(DataSeuratObject,paste0(dir_out,cancer[i],'/subCAF.rds'))
  saveRDS(sc_rds,paste0(dir_sc,file_sc[i]))
  
  rm(DataSeuratObject)
  rm(sc_rds)
  Sys.sleep(5)
  
}




library(tidyverse)
library(ggplot2)
library(ggpubr)


dir_Cyto<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_Cyto<-list.files(pattern = 'new_assigned_locations.csv',path = dir_Cyto,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Cyto,'/'),function(x)x[1]))
###"brca09/slice1"  "gbm02/slice1"   "ipmn01/slice12"
lapply(c("brca09_slice1","gbm02_slice1","ipmn01_slice12"),function(x) grep(x,dataset_slice)) %>% unlist()
file_Cyto<-file_Cyto[-c(10,61,128)]
dataset_slice<-dataset_slice[-c(10,61,128)]

cancer<-unlist(lapply(strsplit(dataset_slice,'_'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()%>% toupper()


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)

for(i in 5:length(cancer)){#i=1
  cancer_file<-file_Cyto[grep(cancer[i],file_Cyto,ignore.case = T)]
  cancer_slice<-unlist(lapply(strsplit(cancer_file,'/'),function(x)x[1]))
  
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  meta_data<-sc_rds@meta.data
  rm(sc_rds)
  
  for(j in 1:length(cancer_file)){#j=1
    st_cyto<-read.csv(paste0(dir_Cyto,cancer_file[j]),stringsAsFactors = F,check.names = F)
    
    st_cyto$subCAF_celltype<-lapply(st_cyto$OriginalCID,function(x){#x=st_cyto$OriginalCID[1]
      cell_name<-unlist(strsplit(x,','))
      cell_celltype<-paste0(meta_data[cell_name,'subCAF'],collapse = ',')
    }) %>% unlist()
    
    write.table(st_cyto,paste0(dir_Cyto,cancer_slice[j],'/new_assigned_subCAF.txt'),quote = F,sep = '\t')
    print(cancer_slice[j])
  }
  Sys.sleep(5)
  
}


all_sliceCAF$cancer<-unlist(lapply(strsplit(all_sliceCAF$slice,'_'),function(x)x[1]))
all_sliceCAF$cancer<-substr(all_sliceCAF$cancer,1,nchar(all_sliceCAF$cancer)-2) %>% toupper()
unique(all_sliceCAF$cellType)
all_sliceCAF$cellType<-factor(all_sliceCAF$cellType,levels = c('apCAF','iCAF',"ifnCAF","rCAF","mCAF","vCAF","dCAF","tCAF"))

p_compare<-ggplot(all_sliceCAF,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("apCAF"="#c26701",'iCAF'='#e2a503',"ifnCAF"="#ffc000","rCAF"="#bbc100",
                               "mCAF"="#007182","vCAF"="#00827c","dCAF"="#4a806c","tCAF"="#737f6b"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subCAF_step1')
print(p_compare)

pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subCAF_bar_step1.pdf',height = 4,width = 6.5)
print(p_compare)
dev.off()
