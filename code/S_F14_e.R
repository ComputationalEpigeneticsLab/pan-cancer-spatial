####计算cellchat
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

dir_out<-'/data/zhouweiwei/ST_analysis/data/ESCC/1re_data/cellchat_Bdy_CAFstep1/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/ESCC/1re_data/ST_expression/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/ESCC/1re_data/copykat/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/ESCC/1re_data/copykat/'


dir_out<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_Bdy_CAFstep1/'
dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
dir_RCTD<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'

file_near<-list.files(pattern = 'nearSpotStep1to20.txt',path = dir_bdy,recursive = T)
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
sample_group<-'1_2_9_10_12_public'
file_RCTD<-list.files(pattern = paste0('Deconvolution_',sample_group,'.txt'),path = dir_RCTD,recursive = T)
dataslice<-unlist(lapply(strsplit(file_RCTD,'_'),function(x)x[1]))
dataset<-unlist(lapply(strsplit(dataslice,'/'),function(x)x[1]))

###############################################################################################
### CAF bdy
i=2
if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))

st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
table(st_near$FinalLocalType)

st_RCTD<-st_RCTD[,colnames(st_RCTD)%in%setdiff(colnames(st_RCTD),c('Epithelial','Endothelial'))]
st_RCTD$celltype<-apply(st_RCTD,1,function(x){
  names(x)<-colnames(st_RCTD)
  x<-x[order(x,decreasing = T)]
  return(names(x)[1])
})

bdy_spot<-st_near$cell.names[st_near$FinalLocalType%in%'Boundary']

bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
bdy_spot_near<-st_near[bdy_spot_near,]
table(bdy_spot_near$FinalLocalType)
bdy_spot_near<-bdy_spot_near$cell.names[bdy_spot_near$FinalLocalType%in%c('Immune','Normal','not.defined')]

step1_near_celltype<-st_RCTD[bdy_spot_near,]
step1_near_celltype<-step1_near_celltype[which(step1_near_celltype$celltype=='CAF'),]

if(length(bdy_spot)>5&&nrow(step1_near_celltype)>1){
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_rds_cellchat<-st_rds[,c(bdy_spot,rownames(step1_near_celltype))]
  
  data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
  
  MP_count<-CreateSeuratObject(data.input)
  MP_count<-AddMetaData(MP_count,metadata = c(rep('bdy',length(bdy_spot)),step1_near_celltype$celltype),col.name = 'labs')
  cellchat <- createCellChat(object=MP_count,group.by = "labs")
  
  groupSize <- as.numeric(table(cellchat@idents))
  print(groupSize)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.human) 
  #cellchat <- smoothData(cellchat, adj=PPI.human) 
  
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
  
  df.net <- subsetCommunication(cellchat,thresh=1)
  write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
  
  df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
  write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
  
}



######################################################################################################
### TAM bdy
dir_out<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_Bdy_TAMstep1/'


i=1
if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))

st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
table(st_near$FinalLocalType)

st_RCTD<-st_RCTD[,colnames(st_RCTD)%in%setdiff(colnames(st_RCTD),c('Epithelial','Endothelial'))]
st_RCTD$celltype<-apply(st_RCTD,1,function(x){
  names(x)<-colnames(st_RCTD)
  x<-x[order(x,decreasing = T)]
  return(names(x)[1])
})

bdy_spot<-st_near$cell.names[st_near$FinalLocalType%in%'Boundary']

bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'step_1'],',')) %>% unique()
bdy_spot_near<-st_near[bdy_spot_near,]
table(bdy_spot_near$FinalLocalType)
bdy_spot_near<-bdy_spot_near$cell.names[bdy_spot_near$FinalLocalType%in%c('Immune','Normal','not.defined')]

step1_near_celltype<-st_RCTD[bdy_spot_near,]
step1_near_celltype<-step1_near_celltype[which(step1_near_celltype$celltype=='TAM'),]

if(length(bdy_spot)>5&&nrow(step1_near_celltype)>1){
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_rds_cellchat<-st_rds[,c(bdy_spot,rownames(step1_near_celltype))]
  
  data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
  
  MP_count<-CreateSeuratObject(data.input)
  MP_count<-AddMetaData(MP_count,metadata = c(rep('bdy',length(bdy_spot)),step1_near_celltype$celltype),col.name = 'labs')
  cellchat <- createCellChat(object=MP_count,group.by = "labs")
  
  groupSize <- as.numeric(table(cellchat@idents))
  print(groupSize)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.human) 
  #cellchat <- smoothData(cellchat, adj=PPI.human) 
  
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
  
  df.net <- subsetCommunication(cellchat,thresh=1)
  write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
  
  df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
  write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
  
}




###############################################################################################
### CAF matrix
dir_out<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_matrix_CAF/'

i=1
if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))

st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
table(st_near$FinalLocalType)

st_RCTD<-st_RCTD[,colnames(st_RCTD)%in%setdiff(colnames(st_RCTD),c('Epithelial','Endothelial'))]
st_RCTD$celltype<-apply(st_RCTD,1,function(x){
  names(x)<-colnames(st_RCTD)
  x<-x[order(x,decreasing = T)]
  return(names(x)[1])
})

st_near<-st_near[rownames(st_RCTD),]
st_RCTD$FinalLocalType<-st_near$FinalLocalType
st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal','not.defined')],]
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
  
  groupSize <- as.numeric(table(cellchat@idents))
  print(groupSize)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.human) 
  #cellchat <- smoothData(cellchat, adj=PPI.human) 
  
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
  
  df.net <- subsetCommunication(cellchat,thresh=1)
  write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
  
  df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
  write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
  
}



######################################################################################################
### TAM matrix
dir_out<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_matrix_TAM/'

i=2
if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))

st_near<-read.delim(paste0(dir_bdy,file_near[i]),stringsAsFactors = F,check.names = F)
st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
table(st_near$FinalLocalType)

st_RCTD<-st_RCTD[,colnames(st_RCTD)%in%setdiff(colnames(st_RCTD),c('Epithelial','Endothelial'))]
st_RCTD$celltype<-apply(st_RCTD,1,function(x){
  names(x)<-colnames(st_RCTD)
  x<-x[order(x,decreasing = T)]
  return(names(x)[1])
})

st_near<-st_near[rownames(st_RCTD),]
st_RCTD$FinalLocalType<-st_near$FinalLocalType
st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal','not.defined')],]
table(st_RCTD$celltype)
CAF_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='TAM')]
Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='CAF'&st_RCTD$celltype!='TAM')]

if(length(CAF_spot)>5&&length(Imm_spot)>5){
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_rds_cellchat<-st_rds[,c(CAF_spot,Imm_spot)]
  
  data.input = Seurat::GetAssayData(st_rds_cellchat, slot = "count", assay = "Spatial") 
  
  MP_count<-CreateSeuratObject(data.input)
  MP_count<-AddMetaData(MP_count,metadata = c(rep('TAM',length(CAF_spot)),st_RCTD[Imm_spot,'celltype']),col.name = 'labs')
  cellchat <- createCellChat(object=MP_count,group.by = "labs")
  
  groupSize <- as.numeric(table(cellchat@idents))
  print(groupSize)
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat) 
  cellchat <- projectData(cellchat, PPI.human) 
  #cellchat <- smoothData(cellchat, adj=PPI.human) 
  
  cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  
  save(cellchat,file = paste0(dir_out,dataslice[i],'.rdata'))
  
  df.net <- subsetCommunication(cellchat,thresh=1)
  write.csv(df.net, paste0(dir_out,dataslice[i],"_net_lr.csv"))
  
  df.netp <- subsetCommunication(cellchat, slot.name = "netP",thresh=1)
  write.csv(df.netp, paste0(dir_out,dataslice[i],"_net_pathway.csv"))
  
}




######################################################################################
###结果统计
dir_cellchat<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_Bdy_CAFstep1/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/cellchat_matrix_CAF/'
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

write.table(all_path_LR,'E:/Mirror/ST_analysis/data/ESCC/pic_data/cellchat/bdy_step1CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/ESCC/pic_data/cellchat/matrix_CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)

dir_cell<-'E:/Mirror/ST_analysis/data/ESCC/pic_data/cellchat/'
bdy_step1CAF_divi<-read.delim(paste0(dir_cell,'bdy_step1CAF_diviLR.txt'),stringsAsFactors=F,check.names=F)
matrix_CAF_divi<-read.delim(paste0(dir_cell,'matrix_CAF_diviLR.txt'),stringsAsFactors = F,check.names = F)


bdy_step1CAF_divi$L_R<-'source'
bdy_step1CAF_divi$L_R[which(bdy_step1CAF_divi$target=='CAF')]<-'target'
bdy_CAFpath<-as.data.frame.array(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$slice))
bdy_CAFpath$sum<-bdy_CAFpath$`BGM/slice1`+bdy_CAFpath$`XYZ/slice1`

matrix_CAF_divi$L_R<-'source'
matrix_CAF_divi$L_R[which(matrix_CAF_divi$target=='CAF')]<-'target'
matrix_CAFpath<-as.data.frame.array(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$slice))
matrix_CAFpath$sum<-matrix_CAFpath$`BGM/slice1`+matrix_CAFpath$`XYZ/slice1`




#######################################################################################################################################
#######################################################################################################################################
###绘制网络图
library(ggraph)
library(tidygraph)
library(igraph)
library(dplyr)

library(svglite)
##wnt通路
net_data<-matrix_CAF_divi[which(matrix_CAF_divi$pathway_name=='WNT'&matrix_CAF_divi$slice=='BGM/slice1'),]
##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)
dir_pic<-'E:/Mirror/ST_analysis/pic/ESCC/cellchat/net/'
#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = 1),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1, 3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=2.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('matrix_step1CAF_WNT')
p_net
svglite(paste0(dir_pic,'BGM_net_matrix_step1CAF_WNT.svg'),width = 4.5,height = 3.5)
print(p_net)
dev.off()



##COLLAGEN通路
net_data<-matrix_CAF_divi[which(matrix_CAF_divi$pathway_name=='COLLAGEN'&matrix_CAF_divi$slice=='BGM/slice1'),]###XYZ  BGM
##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)
dir_pic<-'E:/Mirror/ST_analysis/pic/ESCC/cellchat/net/'
#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = 1),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1, 3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=2.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('matrix_step1CAF_COLLAGEN')
p_net
svglite(paste0(dir_pic,'BGM_net_matrix_step1CAF_COLLAGEN.svg'),width = 4.5,height = 3.5)
print(p_net)
dev.off()