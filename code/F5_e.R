#####挑选配体受体画切片表达
library(scatterpie)
library(ggpubr)
library(png)
library(jsonlite)
library(spacexr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)



convert_to_processed_rds <- function(paths) {
  sapply(paths, function(path) {
    # 分割路径
    parts <- strsplit(path, "/")[[1]]
    
    if (length(parts) == 3) {
      cancer_type <- parts[1]
      dataset <- parts[2]
      sample <- parts[3]
      
      # 检查样本名是否已经包含"processed_"
      if (!grepl("^processed_", sample)) {
        sample <- paste0("processed_", sample)
      }
      
      # 添加.rds后缀
      if (!grepl("\\.rds$", sample)) {
        sample <- paste0(sample, ".rds")
      }
      
      # 重新组合路径
      return(paste(cancer_type, dataset, sample, sep = "/"))
    } else {
      # 如果路径格式不是预期的，返回原路径
      warning(sprintf("路径格式异常: %s", path))
      return(path)
    }
  }, USE.NAMES = FALSE)
}


dir_pic<-'/10X_Visium/plot/cellchat_exp/'
####绘制所有切片结果#########################################################################

dir_near<-'/10X Visium/new_st_copykat/'
file_near<-list.files(pattern = '_BdyCoreBud.txt',path = dir_near,recursive = T)
dir_RCTD<-'/10X Visium/RCTD3/'
file_RCTD<-list.files(pattern = '_Deconvolution.txt',path = dir_RCTD,recursive = T)
dir_rds<-'/10X Visium/new_st/'
file_rds<-list.files(pattern = 'processed_',path = dir_rds,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_RCTD,'_Deconvolution.txt'),function(x)x[1]))


####bdy_step1CAF########################
select_slice_zhutu<-c('BRCA/brca07/slice1','BRCA/brca15/slice1','CRC/crc02/slice1',
                      'CRC/GSE225857/GSM7058757','GBM/GSE237183/GSM7596597','HGSC/GSE263920/GSM8207499',
                      'HGSC/GSE288483/GSM8768682','HNSCC/hnscc01/slice2','LUAD/GSE307534/GSM9226182',
                      'LUAD/luad02/slice2','LUSC/Genomics/10xgenomicsLUSC','OSCC/GSE220978/GSM6833485',
                      'PDAC/GSE235315/GSM7498815','PDAC/GSE274103/GSM8443451','PRAD/GSE278936/GSM8557986')
select_slice_futu<-c('CESC/cesc02/slice2','CRC/crc08/slice2','HB/GSE261958/GSM8155175',
                     'LUAD/GSE307534/GSM9226177','OSCC/GSE220978/GSM6833485','PDAC/GSE294669/GSM8915498',
                     'PDAC/pdac03/slice1','PRAD/GSE278936/GSM8557999','PRAD/GSE278936/GSM8558004',
                     'PRAD/GSE278936/GSM8558009','RCC/GSE250163/GSM7974887','TC/GSE250521/GSM7980860')
select_slice<-select_slice_futu   ####  select_slice_futu 主图/附图
file_near<-paste0(select_slice,"_BdyCoreBud.txt")
#dataSlice<-dataSlice[match(select_slice,dataSlice)]
file_RCTD<-paste0(select_slice,"_Deconvolution.txt")
file_rds <- convert_to_processed_rds(select_slice)


pdf(paste0(dir_pic,'bdy_step1CAF_point.pdf'),width = 5, height = 4)
for(i in 1:length(select_slice)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),c("cell_name", "imagerow","imagecol","near_spot","FinalLocalType")]
  colnames(st_near)

  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
  st_RCTD$celltype<-apply(st_RCTD,1,function(x){
    return(names(x)[which.max(x)])
  })
  table(st_RCTD$celltype)


  p_data<-st_near
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'

  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'near_spot'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  table(bdy_spot_near$FinalLocalType)
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]

  near_CAF<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'CAF'])
  p_data$type[p_data$cell_name%in%near_CAF]<-'step1_CAF'
  p_data$type[p_data$type%in%NA]<-'other'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]

  p_spot2<-ggplot(data = p_data,aes(x = imagerow, y = imagecol,color=type)) +
    geom_point(size=0.7)+
    scale_color_manual(values =  c('Boundary'='#F6AC52','other'='grey95',
                                   'TAM'='#7B3257','step1_TAM'='#C64242',
                                   'CAF'='#A78E41','step1_CAF'='#669999'))+
    theme_classic()+
    ggtitle(select_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=3)))
  print(p_spot2)

  print(select_slice[i])
}
dev.off()


select_LR<-c('ITGB1','COL1A1','COL1A2','SDC4')
####配体受体染色
for(j in 1:length(select_LR)){#j=1
  pdf(paste0(dir_pic,'bdy_step1CAF_',select_LR[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(select_slice)){#i=2
    st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
    st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
    st_RCTD$celltype<-apply(st_RCTD,1,function(x){
      return(names(x)[which.max(x)])
    })
    table(st_RCTD$celltype)

    bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'near_spot'],',')) %>% unique()
    bdy_spot_near<-st_near[bdy_spot_near,]
    table(bdy_spot_near$FinalLocalType)
    bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
    near_CAF<-intersect(bdy_spot_near,rownames(st_RCTD)[st_RCTD$celltype%in%'CAF'])

    if(select_LR[j]%in%rownames(st_rds)){
      p_data<-data.frame(imagerow=st_near$imagerow,imagecol=st_near$imagecol,
                         local=st_near$FinalLocalType,
                         gene=st_rds@assays$Spatial@counts[select_LR[j],])
      p_data<-p_data[which(p_data$local!='not.defined'),]
      p_data$local[rownames(p_data)%in%near_CAF]<-'step1'
      p_data$gene[p_data$local%in%setdiff(unique(p_data$local),c('Boundary','step1'))]<--0.001

      p_data_a<-p_data[which(p_data$gene>=0),]
      Q <- quantile(p_data_a$gene, probs=0.99, na.rm = FALSE)
      p_data_a$gene[which(p_data_a$gene>=Q)]<-round(Q)
      p_data$gene[match(rownames(p_data_a),rownames(p_data))]<-p_data_a$gene

      p<-ggplot(p_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=gene),size=0.7) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(1),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(20),
                                          colorRampPalette(c("#3b9283","#6ab062"))(30),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(50))
        )+
        theme_classic()+
        labs(title = paste0(select_slice[i],'_of_',select_LR[j]),x = "",y = "",col='score')
      print(p)
    }
  }
  dev.off()

}




####################################################################################################
#####bdy_step1TAM#########################
select_slice_zhutu<-c('BRCA/brca12/slice1','BRCA/brca13/slice1','BRCA/GSE242311/GSM7757984',
                      'GC/GSE246011/GSM7853984','HGSC/GSE263920/GSM8207505','HNSCC/GSE281978/GSM8633895',
                      'LIHC/GSE281759/GSM8627611','LIHC/lihc02/slice18','NPC/GSE206245/GSM6248653')
select_slice_futu<-c('LUAD/GSE307534/GSM9226170','LIHC/lihc02/slice8','LIHC/lihc02/slice7',
                     'LIHC/GSE281759/GSM8627606','HNSCC/GSE281978/GSM8633896','HNSCC/GSE281978/GSM8633895',
                     'HGSC/GSE274657/GSM8454235','GBM/gbm05/slice3','BRCA/GSE243275/GSM7782696')

select_slice<-select_slice_futu   ####  select_slice_futu 主图/附图
file_near<-paste0(select_slice,"_BdyCoreBud.txt")
#dataSlice<-dataSlice[match(select_slice,dataSlice)]
file_RCTD<-paste0(select_slice,"_Deconvolution.txt")
file_rds <- convert_to_processed_rds(select_slice)


pdf(paste0(dir_pic,'bdy_step1TAM_point.pdf'),width = 5, height = 4)
for(i in 1:length(select_slice)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),c("cell_name", "imagerow","imagecol","near_spot","FinalLocalType")]
  colnames(st_near)

  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
  st_RCTD$celltype<-apply(st_RCTD,1,function(x){
    return(names(x)[which.max(x)])
  })
  table(st_RCTD$celltype)

  p_data<-st_near
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-NA
  p_data$type[p_data$FinalLocalType%in%c('Boundary')]<-'Boundary'

  bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'near_spot'],',')) %>% unique()
  bdy_spot_near<-st_near[bdy_spot_near,]
  table(bdy_spot_near$FinalLocalType)
  bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]

  near_CAF<-intersect(bdy_spot_near,p_data$cell_name[p_data$cellType%in%'TAM'])
  p_data$type[p_data$cell_name%in%near_CAF]<-'step1_TAM'
  p_data$type[p_data$type%in%NA]<-'other'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]

  p_spot2<-ggplot(data = p_data,aes(x = imagerow, y = imagecol,color=type)) +
    geom_point(size=0.7)+
    scale_color_manual(values =  c('Boundary'='#F6AC52','other'='grey95',
                                   'TAM'='#7B3257','step1_TAM'='#C64242',
                                   'CAF'='#A78E41','step1_CAF'='#669999'))+
    theme_classic()+
    ggtitle(select_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=3)))
  print(p_spot2)

  print(select_slice[i])
}
dev.off()


select_LR<-c('ITGB1','SPP1','ITGA5','ANGPTL4')
####配体受体染色
for(j in 1:length(select_LR)){#j=2
  pdf(paste0(dir_pic,'bdy_step1TAM_',select_LR[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(select_slice)){#i=6
    st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
    st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
    st_RCTD$celltype<-apply(st_RCTD,1,function(x){
      return(names(x)[which.max(x)])
    })
    table(st_RCTD$celltype)
    
    bdy_spot_near<-unlist(strsplit(st_near[which(st_near$FinalLocalType=='Boundary'),'near_spot'],',')) %>% unique()
    bdy_spot_near<-st_near[bdy_spot_near,]
    table(bdy_spot_near$FinalLocalType)
    bdy_spot_near<-bdy_spot_near$cell_name[bdy_spot_near$FinalLocalType%in%c('Immune','Normal')]
    near_TAM<-intersect(bdy_spot_near,rownames(st_RCTD)[st_RCTD$celltype%in%'TAM'])
    
    if(select_LR[j]%in%rownames(st_rds)){
      p_data<-data.frame(imagerow=st_near$imagerow,imagecol=st_near$imagecol,
                         local=st_near$FinalLocalType,
                         gene=st_rds@assays$Spatial@counts[select_LR[j],])
      p_data<-p_data[which(p_data$local!='not.defined'),]
      p_data$local[rownames(p_data)%in%near_TAM]<-'step1'
      p_data$gene[p_data$local%in%setdiff(unique(p_data$local),c('Boundary','step1'))]<--0.001
      
      p_data_a<-p_data[which(p_data$gene>=0),]
      Q <- quantile(p_data_a$gene, probs=0.99, na.rm = FALSE)
      p_data_a$gene[which(p_data_a$gene>=Q)]<-round(Q)
      p_data$gene[match(rownames(p_data_a),rownames(p_data))]<-p_data_a$gene
      
      p<-ggplot(p_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=gene),size=0.7) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(1),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(20),
                                          colorRampPalette(c("#3b9283","#6ab062"))(30),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(50))
        )+
        theme_classic()+
        labs(title = paste0(select_slice[i],'_of_',select_LR[j]),x = "",y = "",col='score')
      print(p)
    }
  }
  dev.off()
  
}



####################################################################################################
#####matrix_CAF#########################
select_slice_zhutu<-c('BRCA/brca01/slice1','BRCA/brca08/slice4','BRCA/brca11/slice1',
                      'CESC/cesc01/slice1','CRC/crc01/slice1','CRC/crc06/slice1',
                      'DSRCT/GSE263523/GSM8279108','DSRCT/GSE263523/GSM8279112','GBM/GSE235672/GSM7507311',
                      'OSCC/oscc01/slice4')
select_slice_futu<-c('TC/GSE250521/GSM7980867','TC/GSE250521/GSM7980865','SKCM/skcm13/slice1',
                     'SKCM/skcm12/slice1','RCC/rcc03/slice1','RCC/GSE175540/GSM5924046',
                     'PRAD/GSE278936/GSM8558018','PRAD/GSE278936/GSM8558006','PRAD/GSE278936/GSM8557995',
                     'HGSC/GSE288483/GSM8768681','GC/GSE251950/GSM7990476')

select_slice<-select_slice_zhutu   ####  select_slice_futu 主图/附图
file_near<-paste0(select_slice,"_BdyCoreBud.txt")
#dataSlice<-dataSlice[match(select_slice,dataSlice)]
file_RCTD<-paste0(select_slice,"_Deconvolution.txt")
file_rds <- convert_to_processed_rds(select_slice)


pdf(paste0(dir_pic,'matrixCAF_point.pdf'),width = 5, height = 4)
for(i in 1:length(select_slice)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),c("cell_name", "imagerow","imagecol","near_spot","FinalLocalType")]
  colnames(st_near)
  
  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
  st_RCTD$celltype<-apply(st_RCTD,1,function(x){
    return(names(x)[which.max(x)])
  })
  table(st_RCTD$celltype)
  
  p_data<-st_near
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other'
  p_data$type[p_data$type%in%'TAM']<-'other'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  
  p_spot2<-ggplot(data = p_data,aes(x = imagerow, y = imagecol,color=type)) + 
    geom_point(size=0.7)+
    scale_color_manual(values =  c('Boundary'='#F6AC52','other'='grey95',
                                   'TAM'='#7B3257','step1_TAM'='#C64242',
                                   'CAF'='#A78E41','step1_CAF'='#669999',
                                   "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                                   "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                                   'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                                   "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',
                                   "GC B cells in the DZ"='#CC9933',
                                   "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',
                                   "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
                                   "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                                   "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066' ))+
    theme_classic()+
    ggtitle(select_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=3)))
  print(p_spot2)
  
  print(select_slice[i])
}
dev.off()



select_LR<-c('FN1','CD44','ITGB2','ICAM1')
####配体受体染色
for(j in 1:length(select_LR)){#j=1
  pdf(paste0(dir_pic,'matrixCAF_',select_LR[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(select_slice)){#i=1
    st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
    st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
    st_RCTD$celltype<-apply(st_RCTD,1,function(x){
      return(names(x)[which.max(x)])
    })
    table(st_RCTD$celltype)
    
    st_RCTD$FinalLocalType<-st_near$FinalLocalType
    st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal')],]
    table(st_RCTD$celltype)
    CAF_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='CAF')]
    Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='CAF'&st_RCTD$celltype!='TAM')]
    
    if(select_LR[j]%in%rownames(st_rds)){
      p_data<-data.frame(imagerow=st_near$imagerow,imagecol=st_near$imagecol,row.names = rownames(st_near),
                         local=st_near$FinalLocalType,
                         gene=st_rds@assays$Spatial@counts[select_LR[j],])
      p_data<-p_data[which(p_data$local!='not.defined'),]
      p_data$gene[!rownames(p_data)%in%c(CAF_spot,Imm_spot)]<--0.001
      
      p_data_a<-p_data[which(p_data$gene>=0),]
      Q <- quantile(p_data_a$gene, probs=0.99, na.rm = FALSE)
      p_data_a$gene[which(p_data_a$gene>=Q)]<-round(Q)
      p_data$gene[match(rownames(p_data_a),rownames(p_data))]<-p_data_a$gene
      
      p<-ggplot(p_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=gene),size=0.7) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(1),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(20),
                                          colorRampPalette(c("#3b9283","#6ab062"))(30),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(50))
        )+
        theme_classic()+
        labs(title = paste0(select_slice[i],'_of_',select_LR[j]),x = "",y = "",col='score')
      print(p)
    }
  }
  dev.off()
  
}



####################################################################################################
#####matrix_TAM#########################
select_slice_zhutu<-c('BRCA/brca08/slice4','CESC/cesc02/slice2','PDAC/GSE274557/GSM8452877',
                      'GBM/GSE235672/GSM7507311','LUAD/luad01/slice2','OVCA/ovca07/slice1')
select_slice_futu<-c('BRCA/GSE242311/GSM7757974','HGSC/GSE263920/GSM8207504','LIHC/lihc02/slice14',
                     'LUAD/luad03/slice1','NPC/GSE206245/GSM6248646','NPC/GSE206245/GSM8361993',
                     'OS/GSE299025/GSM9030939','PDAC/GSE278694/GSM8552951','PDAC/GSE294669/GSM8915494',
                     'PDAC/pdac03/slice1','RCC/rcc01/slice4','TC/GSE250521/GSM7980873')

select_slice<-select_slice_zhutu   ####  select_slice_futu 主图/附图
file_near<-paste0(select_slice,"_BdyCoreBud.txt")
#dataSlice<-dataSlice[match(select_slice,dataSlice)]
file_RCTD<-paste0(select_slice,"_Deconvolution.txt")
file_rds <- convert_to_processed_rds(select_slice)


pdf(paste0(dir_pic,'matrixTAM_point.pdf'),width = 5, height = 4)
for(i in 1:length(select_slice)){
  #i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  st_near<-st_near[rownames(st_RCTD),c("cell_name", "imagerow","imagecol","near_spot","FinalLocalType")]
  colnames(st_near)
  
  st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
  st_RCTD$celltype<-apply(st_RCTD,1,function(x){
    return(names(x)[which.max(x)])
  })
  table(st_RCTD$celltype)
  
  p_data<-st_near
  p_data$cellType<-st_RCTD$celltype
  p_data$type<-p_data$cellType
  p_data$type[p_data$FinalLocalType%in%setdiff(unique(p_data$FinalLocalType),c('Immune','Normal'))]<-'other'
  p_data$type[p_data$type%in%'CAF']<-'other'
  p_data<-p_data[which(p_data$FinalLocalType!='not.defined'),]
  
  p_spot2<-ggplot(data = p_data,aes(x = imagerow, y = imagecol,color=type)) + 
    geom_point(size=0.7)+
    scale_color_manual(values =  c('Boundary'='#F6AC52','other'='grey95',
                                   'TAM'='#7B3257','step1_TAM'='#C64242',
                                   'CAF'='#A78E41','step1_CAF'='#669999',
                                   "B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
                                   "Macrophage"="#FF6666","Myeloid cell"="#FF6600","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
                                   'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
                                   "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'no'='grey80',
                                   "GC B cells in the DZ"='#CC9933',
                                   "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',
                                   "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
                                   "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
                                   "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066' ))+
    theme_classic()+
    ggtitle(select_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=3)))
  print(p_spot2)
  
  print(select_slice[i])
}
dev.off()


select_LR<-c('CD47','THBS1','SDC1')
####配体受体染色
for(j in 1:length(select_LR)){#j=2
  pdf(paste0(dir_pic,'matrixTAM_',select_LR[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(select_slice)){#i=1
    st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
    st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
    st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
    st_RCTD$celltype<-apply(st_RCTD,1,function(x){
      return(names(x)[which.max(x)])
    })
    table(st_RCTD$celltype)
    
    st_RCTD$FinalLocalType<-st_near$FinalLocalType
    st_RCTD<-st_RCTD[rownames(st_RCTD)[st_RCTD$FinalLocalType%in%c('Immune','Normal')],]
    table(st_RCTD$celltype)
    CAF_spot<-rownames(st_RCTD)[which(st_RCTD$celltype=='TAM')]
    Imm_spot<-rownames(st_RCTD)[which(st_RCTD$celltype!='TAM'&st_RCTD$celltype!='CAF')]
    
    if(select_LR[j]%in%rownames(st_rds)){
      p_data<-data.frame(imagerow=st_near$imagerow,imagecol=st_near$imagecol,row.names = rownames(st_near),
                         local=st_near$FinalLocalType,
                         gene=st_rds@assays$Spatial@counts[select_LR[j],])
      p_data<-p_data[which(p_data$local!='not.defined'),]
      p_data$gene[!rownames(p_data)%in%c(CAF_spot,Imm_spot)]<--0.001
      
      p_data_a<-p_data[which(p_data$gene>=0),]
      Q <- quantile(p_data_a$gene, probs=0.99, na.rm = FALSE)
      p_data_a$gene[which(p_data_a$gene>=Q)]<-round(Q)
      p_data$gene[match(rownames(p_data_a),rownames(p_data))]<-p_data_a$gene
      
      p<-ggplot(p_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=gene),size=0.7) +
        # scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
        #                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(1),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(20),
                                          colorRampPalette(c("#3b9283","#6ab062"))(30),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(50))
        )+
        theme_classic()+
        labs(title = paste0(select_slice[i],'_of_',select_LR[j]),x = "",y = "",col='score')
      print(p)
    }
  }
  dev.off()
  
}







