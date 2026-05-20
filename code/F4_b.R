####将所有单细胞数据整合
library(Seurat)


dir_sc<-'/data/zhouweiwei/ST_analysis/data/SC_data/re/SC_re/'
#dir_sc<-'E:/Mirror/ST_analysis/data/SC_data/'
file_sc_rds<-list.files(path = dir_sc,pattern = 'rds',recursive = T)
#file_sc_rds<-file_sc_rds[grep('process',file_sc_rds,invert = T)]
cancer<-unlist(lapply(strsplit(file_sc_rds,'/'),function(x) x[1]))

ifnb.list<-c()
for(i in 1:length(file_sc_rds)){
  #i=2
  st_rds<-readRDS(paste0(dir_sc,file_sc_rds[i]))
  st_rds@meta.data$dataSetID<-cancer[i]
  colnames(st_rds)<-paste0(cancer[i],'_',colnames(st_rds))
  ifnb.list<-as.list(c(ifnb.list,st_rds))
  print(cancer[i])
}
saveRDS(ifnb.list,file = paste0(dir_sc,"all_SC_ifnb.list.rds"))


features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# x<-ifnb.list[[1]]
# table(unlist(lapply(ifnb.list,function(x){length(intersect(features,rownames(x)))})))

ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "cca")
#?FindIntegrationAnchors
saveRDS(immune.anchors,file = paste0(dir_sc,"all_SC_immune.anchors_cca.rds"))
rm(ifnb.list)
#immune.anchors<-readRDS(file=paste0(dir_sc,"all_SC_immune.anchors_cca.rds"))

immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
#saveRDS(immune.combined.sct,file = paste0(dir_out_cancer,cancerType[j],"_immune.combined.sct_Integrate.rds"))

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, verbose = FALSE,resolution = 1)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
saveRDS(immune.combined.sct,file = paste0(dir_sc,"all_SC_immune.combined.sct_Integrate_cca.rds"))
rm(immune.combined.sct)

###rpca
ifnb.list<-readRDS(paste0(dir_sc,"all_SC_ifnb.list.rds"))
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# x<-ifnb.list[[1]]
# table(unlist(lapply(ifnb.list,function(x){length(intersect(features,rownames(x)))})))

ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 10)
saveRDS(immune.anchors,file = paste0(dir_sc,"all_SC_immune.anchors_rpca.rds"))
rm(ifnb.list)
#immune.anchors<-readRDS(file=paste0(dir_out_cancer,cancerType[j],"_immune.anchors_k.anchor_10.rds"))

immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)
#saveRDS(immune.combined.sct,file = paste0(dir_out_cancer,cancerType[j],"_immune.combined.sct_Integrate.rds"))

immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, verbose = FALSE,resolution = 1)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
saveRDS(immune.combined.sct,file = paste0(dir_sc,"all_SC_immune.combined.sct_Integrate_rpca.rds"))










