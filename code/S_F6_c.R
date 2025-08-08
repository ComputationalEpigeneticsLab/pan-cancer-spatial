####对单细胞数据中的TAM细胞分亚型
###从文章中获取5个TAM亚型的marker
#remove.packages(c("Seurat","SeuratObject"))
# Version("Matrix")
# packageDescription('Matrix')
# remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
library(Seurat)
library(SeuratObject)
library(stringr)
library(dplyr)
library(coda)
library(rjags)
library(tidyverse)
library(circlize)
library(hdf5r)
#BiocManager::install("rhdf5")
library(rhdf5)
library(Matrix)


source("E:/Mirror/ST_analysis/program/sctype_score.R")
marker<-read.csv('E:/Mirror/ST_analysis/other_data/subTAM_marker.csv',stringsAsFactors = F,check.names = F)
marker<-marker[!duplicated(marker[,1:2]),]
positive<-lapply(unique(marker$TAM_type),function(x) marker$marker[which(marker$TAM_type==x)])
names(positive)<-unique(marker$TAM_type)
negative<-lapply(positive, function(x) character())
marker_list<-list(positive=positive,
                  negative=negative)


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)
cancer<-unlist(lapply(strsplit(file_sc,'/'),function(x)x[1]))

dir_out<-'E:/Mirror/ST_analysis/data/SC_data/subTAM/'

for(i in 2:length(file_sc)){
  #i= 17 19
  if(!dir.exists(paste0(dir_out,cancer[i]))) dir.create(paste0(dir_out,cancer[i]))
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  #intersect(rownames(sc_rds),'C15ORF48')
  sc_TAM<-subset(sc_rds,subset = seurat_clusters == "TAM")
  
  rm(sc_TAM)
  DataSeuratObject <- CreateSeuratObject(counts = sc_TAM@assays[["RNA"]]@counts, project = "cancer",min.cells = 3)
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
  all_meta_TAM<-all_meta[setdiff(rownames(all_meta),rownames(cell_celltype)),]
  TAM_add<-data.frame(row.names = c(rownames(all_meta_TAM),rownames(cell_celltype)),
                      cellType=c(all_meta_TAM$seurat_clusters,cell_celltype$celltype),
                      lab='add')
  TAM_add<-TAM_add[colnames(sc_rds),]
  
  sc_rds<-AddMetaData(sc_rds,metadata = TAM_add$cellType,col.name = 'subTAM')
  table(sc_rds@meta.data[["seurat_clusters"]],sc_rds@meta.data[["subTAM"]])
  
  meta_data<-sc_rds@meta.data
  meta_data$subCAFTAM<-meta_data$subTAM
  meta_data$subCAFTAM[which(meta_data$subCAFTAM=='CAF')]<-meta_data$subCAF[which(meta_data$subCAFTAM=='CAF')]
  sc_rds<-AddMetaData(sc_rds,metadata = meta_data$subCAFTAM,col.name = 'subCAFTAM')
  print(cancer[i])
  
  saveRDS(DataSeuratObject,paste0(dir_out,cancer[i],'/subTAM.rds'))
  saveRDS(sc_rds,paste0(dir_sc,file_sc[i]))
  
  rm(DataSeuratObject)
  rm(sc_rds)
  Sys.sleep(5)
  
}





####映射5中TAM亚型到空转切片上
library(tidyverse)
library(ggplot2)
library(ggpubr)


dir_Cyto<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_Cyto<-list.files(pattern = 'new_assigned_subCAF.txt',path = dir_Cyto,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Cyto,'/'),function(x)x[1]))
###"brca09/slice1"  "gbm02/slice1"   "ipmn01/slice12"
lapply(c("brca09_slice1","gbm02_slice1","ipmn01_slice12"),function(x) grep(x,dataset_slice)) %>% unlist()
# file_Cyto<-file_Cyto[-c(10,61,128)]
# dataset_slice<-dataset_slice[-c(10,61,128)]

cancer<-unlist(lapply(strsplit(dataset_slice,'_'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()%>% toupper()


dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/re/'
file_sc<-list.files(pattern = 'rds',path = dir_sc,recursive = T)
meta_file<-gsub('.rds','_meta.txt',file_sc)

for(i in 3:length(cancer)){#i=2
  cancer_file<-file_Cyto[grep(cancer[i],file_Cyto,ignore.case = T)]
  cancer_slice<-unlist(lapply(strsplit(cancer_file,'/'),function(x)x[1]))
  
  sc_rds<-readRDS(paste0(dir_sc,file_sc[i]))
  meta_data<-sc_rds@meta.data
  rm(sc_rds)
  write.table(meta_data,paste0(dir_sc,meta_file[i]),quote = F,sep = '\t')
  
  for(j in 1:length(cancer_file)){#j=1
    st_cyto<-read.delim(paste0(dir_Cyto,cancer_file[j]),stringsAsFactors = F,check.names = F)
    st_cyto<-st_cyto[,-4]
    
    st_cyto$subCAFTAM<-lapply(st_cyto$OriginalCID,function(x){#x=st_cyto$OriginalCID[1]
      cell_name<-unlist(strsplit(x,','))
      cell_celltype<-paste0(meta_data[cell_name,'subCAFTAM'],collapse = ',')
    }) %>% unlist()
    
    write.table(st_cyto,paste0(dir_Cyto,cancer_slice[j],'/new_assigned_subCAFTAM.txt'),quote = F,sep = '\t')
    print(cancer_slice[j])
  }
  Sys.sleep(5)
  
}



######绘制比例条形图
dir_subTAM<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_sub<-list.files(pattern = 'new_assigned_subCAFTAM.txt',path=dir_subTAM,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_sub,'/'),function(x)x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-paste0(gsub('_','/',dataset_slice),'_BdyTumorCore.txt')

all_sliceTAM<-c()
for(i in 1:length(file_sub)){
  #i=1
  st_subTAM<-read.delim(paste0(dir_subTAM,file_sub[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[st_subTAM$SpotID,]
  st_bdy$subCAFTAM<-st_subTAM$subCAFTAM
  st_bdy<-st_bdy[st_bdy$FinalLocalType%in%c('Immune','Normal'),]
  
  cellType<-lapply(st_bdy$subCAFTAM,function(x){
    strsplit(x,',') %>% unlist()
  }) %>% unlist()
  cellType<-cellType[grep('TAM',cellType)]
  cellType[grep('FCN1_TAM',cellType)]<-'SPP1_TAM'
  table(cellType)
  
  cellType<-as.data.frame(table(cellType))
  cellType$ratio<-cellType$Freq/sum(cellType$Freq)
  cellType$slice<-dataset_slice[i]
  
  all_sliceTAM<-rbind(all_sliceTAM,cellType)
  print(dataset_slice[i])
}

all_sliceTAM$cancer<-unlist(lapply(strsplit(all_sliceTAM$slice,'_'),function(x)x[1]))
all_sliceTAM$cancer<-substr(all_sliceTAM$cancer,1,nchar(all_sliceTAM$cancer)-2) %>% toupper()
unique(all_sliceTAM$cellType)
write.table(all_sliceTAM,'E:/Mirror/ST_analysis/data/pic_data/subTAM_all.txt',quote = F,sep = '\t',row.names = F)
all_sliceTAM<-read.delim('E:/Mirror/ST_analysis/data/pic_data/subTAM_all.txt',stringsAsFactors = F,check.names = F)

all_sliceTAM$cellType<-factor(all_sliceTAM$cellType,levels = c("SPP1_TAM",'C1Q_TAM',"Microglial_TAM","blood_TAM"))
ppp_data<-aggregate(all_sliceTAM$Freq,by=list(all_sliceTAM$cellType,all_sliceTAM$cancer),sum)
site_cancer<-lapply(unique(ppp_data$Group.2), function(x){
  aaa<-ppp_data[ppp_data$Group.2%in%x,]
  return(aaa[which(aaa[,1]=='SPP1_TAM'),3]-aaa[which(aaa[,1]=='C1Q_TAM'),3])
}) %>% unlist()
names(site_cancer)<-unique(ppp_data$Group.2)
site_cancer<-site_cancer[order(abs(site_cancer),decreasing = F)]

all_sliceTAM$cancer<-factor(all_sliceTAM$cancer,levels = names(site_cancer))

p_compare<-ggplot(all_sliceTAM,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("C1Q_TAM"="#71a3c5",'FCN1_TAM'='#cc889f',"SPP1_TAM"="#b83231",
                               "Microglial_TAM"="#80c598","blood_TAM"="#1d804e"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subTAM')
print(p_compare)

pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subTAM_bar_SPP1&FN1.pdf',height = 4,width = 6.5)
print(p_compare)
dev.off()



#######bdy的step1的TAM亚型比例
######绘制比例条形图
dir_subTAM<-'E:/Mirror/ST_analysis/data/10X Visium/Cytospace/'
file_sub<-list.files(pattern = 'new_assigned_subCAFTAM.txt',path=dir_subTAM,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_sub,'/'),function(x)x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-paste0(gsub('_','/',dataset_slice),'_BdyTumorCore.txt')
dir_near<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
file_near<-paste0(gsub('_','/',dataset_slice),'_BdyNearCore.txt')

all_sliceTAM<-c()
for(i in 1:length(file_sub)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_subTAM<-read.delim(paste0(dir_subTAM,file_sub[i]),stringsAsFactors = F,check.names = F)
  rownames(st_subTAM)<-st_subTAM$SpotID
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  spot_step1<-intersect(spot_step1,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  if(length(spot_step1)>1){
    st_bdy<-st_bdy[st_subTAM$SpotID,]
    st_bdy$subCAFTAM<-st_subTAM$subCAFTAM
    st_bdy<-st_bdy[spot_step1,]
    
    cellType<-lapply(st_bdy$subCAFTAM,function(x){
      strsplit(x,',') %>% unlist()
    }) %>% unlist()
    cellType<-cellType[grep('TAM',cellType)]
    cellType[grep('FCN1_TAM',cellType)]<-'SPP1_TAM'
    table(cellType)
    
    if(length(cellType)>1){
      cellType<-as.data.frame(table(cellType))
      cellType$ratio<-cellType$Freq/sum(cellType$Freq)
      cellType$slice<-dataset_slice[i]
      
      all_sliceTAM<-rbind(all_sliceTAM,cellType)
    }
    
  }
  
  print(dataset_slice[i])
}

all_sliceTAM$cancer<-unlist(lapply(strsplit(all_sliceTAM$slice,'_'),function(x)x[1]))
all_sliceTAM$cancer<-substr(all_sliceTAM$cancer,1,nchar(all_sliceTAM$cancer)-2) %>% toupper()
unique(all_sliceTAM$cellType)

write.table(all_sliceTAM,'E:/Mirror/ST_analysis/data/pic_data/subTAM_step1.txt',quote = F,sep = '\t',row.names = F)
all_sliceTAM<-read.delim('E:/Mirror/ST_analysis/data/pic_data/subTAM_step1.txt',stringsAsFactors = F,check.names = F)

all_sliceTAM$cellType<-factor(all_sliceTAM$cellType,levels = c("SPP1_TAM",'C1Q_TAM',"Microglial_TAM","blood_TAM"))
all_sliceTAM$cancer<-factor(all_sliceTAM$cancer,levels = names(site_cancer))

p_compare<-ggplot(all_sliceTAM,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("C1Q_TAM"="#71a3c5",'FCN1_TAM'='#cc889f',"SPP1_TAM"="#b83231",
                               "Microglial_TAM"="#80c598","blood_TAM"="#1d804e"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subTAM_step1')
print(p_compare)

pdf('E:/Mirror/ST_analysis/pic/subCAFTAM/subTAM_bar_step1_SPP1&FN1.pdf',height = 4,width = 6.5)
print(p_compare)
dev.off()

