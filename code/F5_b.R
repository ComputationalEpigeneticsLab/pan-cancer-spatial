###是在整合的ST数据里 将每个类对应的切片 把每个切片的bdy dis及其对应的step1的spot拿出来
####bdy dis一个标签 step1的spot一个标签 算这两者之间的cellchat
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

###导入配体受体数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

###构建scale数据  所有切片都一样
#image <- readPNG(source = file.path("E:/Mirror/ST_analysis/data/tissue_hires_image.png"))
#scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/1-project/spatial_analysis/cellchat/scalefactors_json.json"))
#scale.factors = jsonlite::fromJSON(txt = file.path("E:/Mirror/ST_analysis/data/scalefactors_json.json"))
# scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))
# 
# spot.size = 65 # the theoretical spot size (um) in 10X Visium
# conversion.factor = spot.size/scale.factors$spot_diameter_fullres
# spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)



dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/'
dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'##### 细胞类型由 195_RCTD点的细胞类型.R 定义

# dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/cellchat_Bdy_CAFstep1/'
# dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
# dir_bdy<-'/data/zhouweiwei/1-project/spatial_analysis/BdyNearCore/'
# dir_RCTD<-'/data/zhouweiwei/1-project/spatial_analysis/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/cellchat_Bdy_CAFstep1/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/RCTD_celltype/'

# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'

file_near<-list.files(pattern = 'BdyNearCore.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_near,'_'),function(x)x[1]))
table(paste0(dataslice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1])))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
slice_rda<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))

for(i in 1:length(file_rds)){
  #i=1
  if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  tryCatch({
    st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    table(st_near$FinalLocalType)
    
    bdy_spot<-st_near$cell_name[st_near$FinalLocalType%in%'Boundary']
    
    bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
    bdy_spot_near<-st_near[bdy_spot_near,]
    table(bdy_spot_near$FinalLocalType)
    bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
    
    step1_near_celltype<-st_RCTD[bdy_spot_near,]
    #step1_near_celltype<-step1_near_celltype[which(step1_near_celltype$celltype=='CAF'),]
    
    if(length(bdy_spot)>5&&nrow(step1_near_celltype)>5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_rds_cellchat<-st_rds[,c(bdy_spot,rownames(step1_near_celltype))]
      
      data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
      
      MP_count<-CreateSeuratObject(data.input)
      MP_count<-AddMetaData(MP_count,metadata = c(rep('bdy',length(bdy_spot)),step1_near_celltype$celltype),col.name = 'labs')
      cellchat <- createCellChat(object=MP_count,group.by = "labs")
      
      # #meta信息
      # meta = data.frame(labels = c(rep('bdy',length(bdy_spot)),rep('CAF',nrow(step1_near_celltype))), #名字自定义
      #                   samples='sample1',
      #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # # meta = data.frame(labels = as.vector(st_rds_cellchat@meta.data$predicted.celltype), #名字自定义
      # #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # meta$samples<-as.factor(meta$samples)
      # unique(meta$labels)
      # # 空间图像信息
      # spatial.locs = Seurat::GetTissueCoordinates(st_rds_cellchat, scale = NULL, cols = c("imagerow", "imagecol"))
      # 
      # cellchat <- createCellChat(object = data.input, 
      #                            meta = meta, 
      #                            group.by = "labels", #前面的meta ，定义的名字是labels
      #                            datatype = "spatial", ###
      #                            coordinates = spatial.locs, 
      #                            spatial.factors = spatial.factors)
      #?createCellChat
      groupSize <- as.numeric(table(cellchat@idents))
      print(groupSize)
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat) 
      #cellchat <- projectData(cellchat, PPI.human) 
      cellchat <- smoothData(cellchat, adj=PPI.human) 
      #?projectData
      #?smoothData
      
      cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE,contact.range=100,contact.knn.k = 6) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 5)
      
      ##CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率
      ##推断信号通路水平的细胞通讯网络
      ##可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
      cellchat <- computeCommunProbPathway(cellchat)
      #saveRDS(cellchat,file = paste0(dir_out,dataslice[i],'_cellchat.rds'))
      
      save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
      #load('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/brca01/test.rdata')
      
      #cellchat_test<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/brca01/slice3_cellchat.rds')
      
      df.net <- subsetCommunication(cellchat,thresh=1)
      write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
      
      df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
      write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
      
    }
  },error = function(e){
    stopMessage_bbb<-"Unable to calculate cellchat"
    return(stopMessage_bbb)
  })
  
  print(dataslice[i])
}


###是在整合的ST数据里 将每个类对应的切片 把每个切片的bdy dis及其对应的step1的spot拿出来
####bdy dis一个标签 step1的spot一个标签 算这两者之间的cellchat
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

###导入配体受体数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

###构建scale数据  所有切片都一样
#image <- readPNG(source = file.path("E:/Mirror/ST_analysis/data/tissue_hires_image.png"))
#scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/1-project/spatial_analysis/cellchat/scalefactors_json.json"))
#scale.factors = jsonlite::fromJSON(txt = file.path("E:/Mirror/ST_analysis/data/scalefactors_json.json"))
scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))

spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scale.factors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)



dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_CAF/'
dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'##### 细胞类型由 195_RCTD点的细胞类型.R 定义

# dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/cellchat_matrix_CAF/'
# dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
# dir_bdy<-'/data/zhouweiwei/1-project/spatial_analysis/BdyNearCore/'
# dir_RCTD<-'/data/zhouweiwei/1-project/spatial_analysis/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/spatial/data/10X_Visium/cellchat_matrix_CAF/'
dir_rds<-'/data/zhouweiwei/spatial/data/10X_Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/spatial/data/10X_Visium/nearMax_spot/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/spatial/data/10X_Visium/RCTD_celltype/'

# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'

file_near<-list.files(pattern = 'BdyNearCore.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_near,'_'),function(x)x[1]))
table(paste0(dataslice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1])))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
slice_rda<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))

for(i in 1:length(file_rds)){
  #i=1
  if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  tryCatch({
    st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    table(st_near$FinalLocalType)
    st_near<-st_near[rownames(st_RCTD),]
    st_RCTD$FinalLocalType<-st_near$FinalLocalType
    st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal')],]
    table(st_RCTD$celltype)
    CAF_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='CAF')]
    Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='CAF'&st_RCTD$celltype!='TAM')]
    
    if(length(CAF_spot)>5&&length(Imm_spot)>5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_rds_cellchat<-st_rds[,c(CAF_spot,Imm_spot)]
      
      data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
      
      MP_count<-CreateSeuratObject(data.input)
      MP_count<-AddMetaData(MP_count,metadata = c(rep('CAF',length(CAF_spot)),st_RCTD[Imm_spot,'celltype']),col.name = 'labs')
      cellchat <- createCellChat(object=MP_count,group.by = "labs")
      
      # #meta信息
      # meta = data.frame(labels = c(rep('bdy',length(bdy_spot)),rep('CAF',nrow(step1_near_celltype))), #名字自定义
      #                   samples='sample1',
      #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # # meta = data.frame(labels = as.vector(st_rds_cellchat@meta.data$predicted.celltype), #名字自定义
      # #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # meta$samples<-as.factor(meta$samples)
      # unique(meta$labels)
      # # 空间图像信息
      # spatial.locs = Seurat::GetTissueCoordinates(st_rds_cellchat, scale = NULL, cols = c("imagerow", "imagecol"))
      # 
      # cellchat <- createCellChat(object = data.input, 
      #                            meta = meta, 
      #                            group.by = "labels", #前面的meta ，定义的名字是labels
      #                            datatype = "spatial", ###
      #                            coordinates = spatial.locs, 
      #                            spatial.factors = spatial.factors)
      #?createCellChat
      groupSize <- as.numeric(table(cellchat@idents))
      print(groupSize)
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat) 
      #cellchat <- projectData(cellchat, PPI.human) 
      cellchat <- smoothData(cellchat, adj=PPI.human) 
      #?projectData
      #?smoothData
      
      cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE,contact.range=100,contact.knn.k = 6) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 5)
      
      ##CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率
      ##推断信号通路水平的细胞通讯网络
      ##可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
      cellchat <- computeCommunProbPathway(cellchat)
      #saveRDS(cellchat,file = paste0(dir_out,dataslice[i],'_cellchat.rds'))
      
      save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
      #load('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/brca01/test.rdata')
      
      #cellchat_test<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/brca01/slice3_cellchat.rds')
      
      df.net <- subsetCommunication(cellchat,thresh=1)
      write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
      
      df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
      write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
      
    }
  },error = function(e){
    stopMessage_bbb<-"Unable to calculate cellchat"
    return(stopMessage_bbb)
  })
  
  print(dataslice[i])
}



###是在整合的ST数据里 将每个类对应的切片 把每个切片的bdy dis及其对应的step1的spot拿出来
####bdy dis一个标签 step1的spot一个标签 算这两者之间的cellchat
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

###导入配体受体数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

###构建scale数据  所有切片都一样
#image <- readPNG(source = file.path("E:/Mirror/ST_analysis/data/tissue_hires_image.png"))
#scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/1-project/spatial_analysis/cellchat/scalefactors_json.json"))
#scale.factors = jsonlite::fromJSON(txt = file.path("E:/Mirror/ST_analysis/data/scalefactors_json.json"))
scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))

spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scale.factors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)



dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/'
dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'##### 细胞类型由 195_RCTD点的细胞类型.R 定义

# dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/cellchat_Bdy_TAMstep1/'
# dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
# dir_bdy<-'/data/zhouweiwei/1-project/spatial_analysis/BdyNearCore/'
# dir_RCTD<-'/data/zhouweiwei/1-project/spatial_analysis/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/spatial/data/10X_Visium/cellchat_Bdy_TAMstep1/'
dir_rds<-'/data/zhouweiwei/spatial/data/10X_Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/spatial/data/10X_Visium/nearMax_spot/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/spatial/data/10X_Visium/RCTD_celltype/'

# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'

file_near<-list.files(pattern = 'BdyNearCore.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_near,'_'),function(x)x[1]))
table(paste0(dataslice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1])))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
slice_rda<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))

for(i in 1:length(file_rds)){
  #i=3
  if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  tryCatch({
    st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    table(st_near$FinalLocalType)
    
    bdy_spot<-st_near$cell_name[st_near$FinalLocalType%in%'Boundary']
    
    bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
    bdy_spot_near<-st_near[bdy_spot_near,]
    table(bdy_spot_near$FinalLocalType)
    bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
    
    step1_near_celltype<-st_RCTD[bdy_spot_near,]
    step1_near_celltype<-step1_near_celltype[which(step1_near_celltype$celltype=='TAM'),]
    
    if(length(bdy_spot)>5&&nrow(step1_near_celltype)>5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_rds_cellchat<-st_rds[,c(bdy_spot,rownames(step1_near_celltype))]
      
      data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
      
      MP_count<-CreateSeuratObject(data.input)
      MP_count<-AddMetaData(MP_count,metadata = c(rep('bdy',length(bdy_spot)),rep('TAM',nrow(step1_near_celltype))),col.name = 'labs')
      cellchat <- createCellChat(object=MP_count,group.by = "labs")
      
      # #meta信息
      # meta = data.frame(labels = c(rep('bdy',length(bdy_spot)),rep('TAM',nrow(step1_near_celltype))), #名字自定义
      #                   samples='sample1',
      #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # # meta = data.frame(labels = as.vector(st_rds_cellchat@meta.data$predicted.celltype), #名字自定义
      # #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # meta$samples<-as.factor(meta$samples)
      # unique(meta$labels)
      # # 空间图像信息
      # spatial.locs = Seurat::GetTissueCoordinates(st_rds_cellchat, scale = NULL, cols = c("imagerow", "imagecol"))
      # 
      # cellchat <- createCellChat(object = data.input, 
      #                            meta = meta, 
      #                            group.by = "labels", #前面的meta ，定义的名字是labels
      #                            datatype = "spatial", ###
      #                            coordinates = spatial.locs, 
      #                            spatial.factors = spatial.factors)
      #?createCellChat
      groupSize <- as.numeric(table(cellchat@idents))
      print(groupSize)
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat) 
      #cellchat <- projectData(cellchat, PPI.human) 
      cellchat <- smoothData(cellchat, adj=PPI.human) 
      #?projectData
      #?smoothData
      
      cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE,contact.range=100,contact.knn.k = 6) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 5)
      
      ##CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率
      ##推断信号通路水平的细胞通讯网络
      ##可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
      cellchat <- computeCommunProbPathway(cellchat)
      #saveRDS(cellchat,file = paste0(dir_out,dataslice[i],'_cellchat.rds'))
      
      save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
      #load('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/test.rdata')
      
      #cellchat_test<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/slice3_cellchat.rds')
      
      df.net <- subsetCommunication(cellchat,thresh=1)
      write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
      
      df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
      write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
      
    }
  },error = function(e){
    stopMessage_bbb<-"Unable to calculate cellchat"
    return(stopMessage_bbb)
  })
  
  print(dataslice[i])
}






###是在整合的ST数据里 将每个类对应的切片 把每个切片的bdy dis及其对应的step1的spot拿出来
####bdy dis一个标签 step1的spot一个标签 算这两者之间的cellchat
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

###导入配体受体数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

###构建scale数据  所有切片都一样
#image <- readPNG(source = file.path("E:/Mirror/ST_analysis/data/tissue_hires_image.png"))
#scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/1-project/spatial_analysis/cellchat/scalefactors_json.json"))
#scale.factors = jsonlite::fromJSON(txt = file.path("E:/Mirror/ST_analysis/data/scalefactors_json.json"))
# scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))
# scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))
# scale.factors$spot_diameter_fullres==14.16615
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/14.16615
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)



dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'##### 细胞类型由 195_RCTD点的细胞类型.R 定义

# dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/cellchat_matrix_TAM/'
# dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
# dir_bdy<-'/data/zhouweiwei/1-project/spatial_analysis/BdyNearCore/'
# dir_RCTD<-'/data/zhouweiwei/1-project/spatial_analysis/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/spatial/data/10X_Visium/cellchat_matrix_TAM/'
dir_rds<-'/data/zhouweiwei/spatial/data/10X_Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/spatial/data/10X_Visium/nearMax_spot/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/spatial/data/10X_Visium/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X Visium/RCTD_celltype/'

# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'

file_near<-list.files(pattern = 'BdyNearCore.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_near,'_'),function(x)x[1]))
table(paste0(dataslice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1])))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
slice_rda<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))

for(i in 1:length(file_rds)){
  #i=1
  if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  tryCatch({
    st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    table(st_near$FinalLocalType)
    st_near<-st_near[rownames(st_RCTD),]
    st_RCTD$FinalLocalType<-st_near$FinalLocalType
    st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal')],]
    table(st_RCTD$celltype)
    TAM_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='TAM')]
    Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='CAF'&st_RCTD$celltype!='TAM')]
    
    if(length(TAM_spot)>5&&length(Imm_spot)>5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_rds_cellchat<-st_rds[,c(TAM_spot,Imm_spot)]
      
      data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
      
      MP_count<-CreateSeuratObject(data.input)
      MP_count<-AddMetaData(MP_count,metadata = c(rep('TAM',length(TAM_spot)),st_RCTD[Imm_spot,'celltype']),col.name = 'labs')
      cellchat <- createCellChat(object=MP_count,group.by = "labs")
      
      # #meta信息
      # meta = data.frame(labels = c(rep('bdy',length(bdy_spot)),rep('TAM',nrow(step1_near_celltype))), #名字自定义
      #                   samples='sample1',
      #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # # meta = data.frame(labels = as.vector(st_rds_cellchat@meta.data$predicted.celltype), #名字自定义
      # #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # meta$samples<-as.factor(meta$samples)
      # unique(meta$labels)
      # # 空间图像信息
      # spatial.locs = Seurat::GetTissueCoordinates(st_rds_cellchat, scale = NULL, cols = c("imagerow", "imagecol"))
      # 
      # cellchat <- createCellChat(object = data.input, 
      #                            meta = meta, 
      #                            group.by = "labels", #前面的meta ，定义的名字是labels
      #                            datatype = "spatial", ###
      #                            coordinates = spatial.locs, 
      #                            spatial.factors = spatial.factors)
      #?createCellChat
      groupSize <- as.numeric(table(cellchat@idents))
      print(groupSize)
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat) 
      #cellchat <- projectData(cellchat, PPI.human) 
      cellchat <- smoothData(cellchat, adj=PPI.human) 
      #?projectData
      #?smoothData
      
      cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE,contact.range=100,contact.knn.k = 6) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 5)
      
      ##CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率
      ##推断信号通路水平的细胞通讯网络
      ##可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
      cellchat <- computeCommunProbPathway(cellchat)
      #saveRDS(cellchat,file = paste0(dir_out,dataslice[i],'_cellchat.rds'))
      
      save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
      #load('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/test.rdata')
      
      #cellchat_test<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/slice3_cellchat.rds')
      
      df.net <- subsetCommunication(cellchat,thresh=1)
      write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
      
      df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
      write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
      
    }
  },error = function(e){
    stopMessage_bbb<-"Unable to calculate cellchat"
    return(stopMessage_bbb)
  })
  
  print(dataslice[i])
}






###是在整合的ST数据里 将每个类对应的切片 把每个切片的bdy dis及其对应的step1的spot拿出来
####bdy dis一个标签 step1的spot一个标签 算这两者之间的cellchat
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

###导入配体受体数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

###构建scale数据  所有切片都一样
#image <- readPNG(source = file.path("E:/Mirror/ST_analysis/data/tissue_hires_image.png"))
#scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/1-project/spatial_analysis/cellchat/scalefactors_json.json"))
#scale.factors = jsonlite::fromJSON(txt = file.path("E:/Mirror/ST_analysis/data/scalefactors_json.json"))
# scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))
# scale.factors = jsonlite::fromJSON(txt = file.path("/data/zhouweiwei/spatial/data/scalefactors_json.json"))
# scale.factors$spot_diameter_fullres==14.16615
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/14.16615
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)



dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'##### 细胞类型由 195_RCTD点的细胞类型.R 定义

# dir_out<-'/data/zhouweiwei/1-project/spatial_analysis/cellchat_matrix_TAM/'
# dir_rds<-'/data/zhouweiwei/1-project/spatial_analysis/0_ST_data/'
# dir_bdy<-'/data/zhouweiwei/1-project/spatial_analysis/BdyNearCore/'
# dir_RCTD<-'/data/zhouweiwei/1-project/spatial_analysis/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/spatial/data/10X_Visium/cellchat_matrix_TAM/'
dir_rds<-'/data/zhouweiwei/spatial/data/10X_Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/spatial/data/10X_Visium/nearMax_spot/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/spatial/data/10X_Visium/RCTD_celltype/'

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X Visium/ST_expression/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X Visium/BdyNearCore/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X Visium/RCTD_celltype/'

# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'

file_near<-list.files(pattern = 'BdyNearCore.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_near,'_'),function(x)x[1]))
table(paste0(dataslice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
table(dataslice==unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1])))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))
slice_rda<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[2]))

for(i in 1:length(file_rds)){
  #i=1
  if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  tryCatch({
    st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    table(st_near$FinalLocalType)
    st_near<-st_near[rownames(st_RCTD),]
    st_RCTD$FinalLocalType<-st_near$FinalLocalType
    st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal')],]
    table(st_RCTD$celltype)
    TAM_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='TAM')]
    Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='CAF'&st_RCTD$celltype!='TAM')]
    
    if(length(TAM_spot)>5&&length(Imm_spot)>5){
      st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
      st_rds_cellchat<-st_rds[,c(TAM_spot,Imm_spot)]
      
      data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
      
      MP_count<-CreateSeuratObject(data.input)
      MP_count<-AddMetaData(MP_count,metadata = c(rep('TAM',length(TAM_spot)),st_RCTD[Imm_spot,'celltype']),col.name = 'labs')
      cellchat <- createCellChat(object=MP_count,group.by = "labs")
      
      # #meta信息
      # meta = data.frame(labels = c(rep('bdy',length(bdy_spot)),rep('TAM',nrow(step1_near_celltype))), #名字自定义
      #                   samples='sample1',
      #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # # meta = data.frame(labels = as.vector(st_rds_cellchat@meta.data$predicted.celltype), #名字自定义
      # #                   row.names = colnames(st_rds_cellchat)) # manually create a dataframe consisting of the cell labels
      # meta$samples<-as.factor(meta$samples)
      # unique(meta$labels)
      # # 空间图像信息
      # spatial.locs = Seurat::GetTissueCoordinates(st_rds_cellchat, scale = NULL, cols = c("imagerow", "imagecol"))
      # 
      # cellchat <- createCellChat(object = data.input, 
      #                            meta = meta, 
      #                            group.by = "labels", #前面的meta ，定义的名字是labels
      #                            datatype = "spatial", ###
      #                            coordinates = spatial.locs, 
      #                            spatial.factors = spatial.factors)
      #?createCellChat
      groupSize <- as.numeric(table(cellchat@idents))
      print(groupSize)
      cellchat@DB <- CellChatDB.use
      cellchat <- subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat) 
      #cellchat <- projectData(cellchat, PPI.human) 
      cellchat <- smoothData(cellchat, adj=PPI.human) 
      #?projectData
      #?smoothData
      
      cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE,contact.range=100,contact.knn.k = 6) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
      # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
      cellchat <- filterCommunication(cellchat, min.cells = 5)
      
      ##CellChat 通过总结与每个信号通路相关的所有配体-受体相互作用的通信概率，来计算信号通路级别上的通信概率
      ##推断信号通路水平的细胞通讯网络
      ##可以通过计算链路的数量或汇总通信概率来计算细胞间的聚合通信网络
      cellchat <- computeCommunProbPathway(cellchat)
      #saveRDS(cellchat,file = paste0(dir_out,dataslice[i],'_cellchat.rds'))
      
      save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
      #load('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/test.rdata')
      
      #cellchat_test<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/brca01/slice3_cellchat.rds')
      
      df.net <- subsetCommunication(cellchat,thresh=1)
      write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
      
      df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
      write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
      
    }
  },error = function(e){
    stopMessage_bbb<-"Unable to calculate cellchat"
    return(stopMessage_bbb)
  })
  
  print(dataslice[i])
}


###cellchat分开统计
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(reshape2)


####
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_CAF/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_cellchat,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for(i in 1:length(file_LR)){
  #i=1
  LR_data<-read.csv(paste0(dir_cellchat,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  LR_data$slice<-dataSlice[i]
  colnames(LR_data)
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  LR_data<-LR_data[,c('source','target','pathway_name','ligand','receptor','annotation','evidence','slice')]
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    other<-lapply(setdiff(unique(c(LR_data$source,LR_data$target)),'CAF'),function(y){
      return(y%in%a)
    }) %>% unlist()
    return('CAF'%in%a&&TRUE%in%other)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  if(nrow(LR_data)>0){
    all_path_LR<-rbind(all_path_LR,LR_data)
  }
  print(dataSlice[i])
}
all_path_LR$pathway_name<-as.vector(all_path_LR$pathway_name)

write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)

# bdy_step1CAF<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF.txt',stringsAsFactors = F,check.names = F)
# sum(bdy_step1CAF$LR_num)
# 
# class(all_path_LR$pathway_name)
# a_all_path_LR<-all_path_LR[all_path_LR$slice%in%'brca01/slice1',]
# b_bdy_step1CAF<-bdy_step1CAF[bdy_step1CAF$slice%in%'brca01/slice1',]
#############################################################################要注意 不同细胞对之间可能有相同的通路和LR对
#########################################################################不能用table计算频率 重复出现的LR会导致pathway少

bdy_step1CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF_diviLR.txt',stringsAsFactors=F,check.names=F)
bdy_step1TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM_diviLR.txt',stringsAsFactors=F,check.names=F)

matrix_CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_diviLR.txt',stringsAsFactors = F,check.names = F)
matrix_TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_diviLR.txt',stringsAsFactors = F,check.names = F)
#sum(matrix_CAF_divi$LR_num)

bdy_step1CAF_divi$L_R<-'source'
bdy_step1CAF_divi$L_R[which(bdy_step1CAF_divi$target=='CAF')]<-'target'
bdy_CAFpath<-as.data.frame.array(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$L_R))
bdy_CAFpath$sum<-bdy_CAFpath$source+bdy_CAFpath$target

bdy_step1TAM_divi$L_R<-'source'
bdy_step1TAM_divi$L_R[which(bdy_step1TAM_divi$target=='TAM')]<-'target'
bdy_TAMpath<-as.data.frame.array(table(bdy_step1TAM_divi$pathway_name,bdy_step1TAM_divi$L_R))
bdy_TAMpath$sum<-bdy_TAMpath$source+bdy_TAMpath$target

matrix_CAF_divi$L_R<-'source'
matrix_CAF_divi$L_R[which(matrix_CAF_divi$target=='CAF')]<-'target'
matrix_CAFpath<-as.data.frame.array(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$L_R))
matrix_CAFpath$sum<-matrix_CAFpath$source+matrix_CAFpath$target

matrix_TAM_divi$L_R<-'source'
matrix_TAM_divi$L_R[which(matrix_TAM_divi$target=='TAM')]<-'target'
matrix_TAMpath<-as.data.frame.array(table(matrix_TAM_divi$pathway_name,matrix_TAM_divi$L_R))
matrix_TAMpath$sum<-matrix_TAMpath$source+matrix_TAMpath$target












