library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(png)

package.version('CellChat')

aaa<-readRDS('/10X_Visium/new_st_cellchat_ST/BRCA/brca01/slice1.rds')

col_data<-c('Core'='#d62d28','Boundary'='#f6b86d','Budding'='#ee762d','Immune'='#a9d38a',
            "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
            "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
            'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
            "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',"GC B cells in the DZ"='#CC9933',
            "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',"Treg"='#FFCC33',
            "Cytotoxic"='#FF6666',"TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
            "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066')             

# Spatial plot
# aaa@netP$pathways
pathways.show <- c("COLLAGEN")
# levels(aaa@idents) 
# 
par(mfrow=c(1,1))
netVisual_aggregate(aaa, signaling = pathways.show, layout = "spatial",
                    edge.width.max = 3, vertex.size.max = 2,
                    alpha.image = 0.6, vertex.label.cex = 3.5,
                    color.use = col_data[levels(aaa@idents) ])


# 设置输出路径
dir_cellchat<-'/10X_Visium/new_st_cellchat_ST/'
file_rds<-list.files(pattern = 'rds',path = dir_cellchat,recursive = T)
dataslice<-gsub('.rds','',file_rds)
dataslice <- gsub('/','_',dataslice)



#COLLAGEN
##透明度1
output_path <- "/10X_Visium/plot/cellchat_new/COLLAGEN/"

for(i in 1:length(dataslice)){
  pdf(paste0(output_path, dataslice[i], ".pdf"), width = 10, height = 8)
  
  aaa <- readRDS(paste0(dir_cellchat, file_rds[i]))
  
  col_data <- c('Core'='#d62d28','Boundary'='#f6b86d','Budding'='#ee762d','Immune'='#a9d38a',
                "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',"GC B cells in the DZ"='#CC9933',
                "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',"Treg"='#FFCC33',
                "Cytotoxic"='#FF6666',"TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066')
  
  pathways.show <- c("COLLAGEN")
  par(mfrow = c(1, 1))
  
  # 执行绘图
  p <- netVisual_aggregate(aaa, signaling = pathways.show, layout = "spatial",
                           edge.width.max = 3, vertex.size.max = 2,
                           alpha.image = 1, vertex.label.cex = 3.5,
                           color.use = col_data[levels(aaa@idents)])
  
  # 显示图形
  print(p)
  
  # 强制刷新
  dev.flush()
  
  # 短暂延迟（如果需要）
  Sys.sleep(0.5)
  
  cat("完成:", dataslice[i], "\n")
  dev.off()
}



#LAMININ
##透明度1
output_path <- "/10X_Visium/plot/cellchat_new/LAMININ/"
for(i in 1:length(dataslice)){
  pdf(paste0(output_path, dataslice[i], ".pdf"), width = 10, height = 8)
  
  aaa <- readRDS(paste0(dir_cellchat, file_rds[i]))
  
  col_data <- c('Core'='#d62d28','Boundary'='#f6b86d','Budding'='#ee762d','Immune'='#a9d38a',
                "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',"GC B cells in the DZ"='#CC9933',
                "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',"Treg"='#FFCC33',
                "Cytotoxic"='#FF6666',"TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066')
  
  pathways.show <- c("LAMININ")
  par(mfrow = c(1, 1))
  
  # 执行绘图
  p <- netVisual_aggregate(aaa, signaling = pathways.show, layout = "spatial",
                           edge.width.max = 3, vertex.size.max = 2,
                           alpha.image = 1, vertex.label.cex = 3.5,
                           color.use = col_data[levels(aaa@idents)])
  
  # 显示图形
  print(p)
  
  # 强制刷新
  dev.flush()
  
  # 短暂延迟（如果需要）
  Sys.sleep(0.5)
  
  cat("完成:", dataslice[i], "\n")
  dev.off()
}



