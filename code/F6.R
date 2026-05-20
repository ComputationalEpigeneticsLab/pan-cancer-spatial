####slingshot后续分析挑选与轨迹进化相关的基因
library(slingshot)
#BiocManager::install("BUSpaRse")
library(BUSpaRse)
library(tidyverse)
#install.packages('tidymodels')
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
#library(Matrix)
library(tradeSeq)
library(SingleCellExperiment)
library(RColorBrewer)
library(DelayedMatrixStats)
library(ggplot2)
library(ggpubr)

col_all<-c("B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
           "Macrophage"="#CC3333","Myeloid cell"="#ED703F","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
           'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
           "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00","GC B cells in the DZ"='#CC9933',
           "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#A6B864',
           "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
           "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
           "Naive"='#E0B8B6',"B cell Regulatory"='#990033',"Naive B cell"='#990066' ,
           'Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d')
colors<-data.frame(celltype=names(col_all),col=col_all)

dir_sling<-'/data/10X_Visium/plot/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice_CAF<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice_CAF,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice_CAF==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
dir_pic<-'/data/10X_Visium/plot/plot_slingshot/'
cancer<-unlist(lapply(strsplit(dataSlice_CAF,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

pdf(paste0(dir_pic,'lineage_celltype_CAF.pdf'),width = 7,height = 6)
for(i in 1:length(file_CAF_slingshot)){
  #i=1
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[i]))
  colnames(colData(sim))
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')#[1]
  plot_lab<-paste0("lineage",1:length(SlingshotDataSet(sim)@lineages))
  plot_lab[CAF_site]<-'CAF_lineage'

  
  colors3<-colors[colors$celltype%in%intersect(colors$celltype,colData(sim)$celltype),]
  sim_meta<-data.frame(spot=rownames(colData(sim)),celltype=colData(sim)$celltype)
  plotcol_2<-merge(sim_meta,colors3,by='celltype',all=T)
  rownames(plotcol_2)<-plotcol_2$spot
  coor<-reducedDims(sim)$UMAP
  plotcol_2<-plotcol_2[rownames(coor),]
  
  plot(reducedDims(sim)$UMAP, col = plotcol_2$col, pch=16, asp = 1,cex=0.6,main=dataSlice_CAF[i])
  lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(length(SlingshotDataSet(sim)@lineages),"Set1"))
  legend("left",
         legend = plot_lab,
         col = unique(brewer.pal(length(SlingshotDataSet(sim)@lineages),"Set1")),
         inset=0.8,
         pch = 16)
  
  print(dataSlice_CAF[i])
}
dev.off()

file_TAM_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_TAM_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_TAM_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice_TAM<-unlist(lapply(strsplit(file_TAM_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice_TAM,unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
table(dataSlice_TAM==unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))

cancer<-unlist(lapply(strsplit(dataSlice_TAM,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

pdf(paste0(dir_pic,'lineage_celltype_TAM.pdf'),width = 7,height = 6)
for(i in 1:length(file_TAM_slingshot)){
  #i=1
  sim<-readRDS(paste0(dir_sling,file_TAM_slingshot[i]))
  colnames(colData(sim))
  
  lineages<-SlingshotDataSet(sim)@lineages
  TAM_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='TAM') rr<-'Y'
    return(rr)
  }) %>% unlist()
  TAM_site<-which(TAM_site=='Y')#[1]
  plot_lab<-paste0("lineage",1:length(SlingshotDataSet(sim)@lineages))
  plot_lab[TAM_site]<-'TAM_lineage'

  
  colors3<-colors[colors$celltype%in%intersect(colors$celltype,colData(sim)$celltype),]
  sim_meta<-data.frame(spot=rownames(colData(sim)),celltype=colData(sim)$celltype)
  plotcol_2<-merge(sim_meta,colors3,by='celltype',all=T)
  rownames(plotcol_2)<-plotcol_2$spot
  coor<-reducedDims(sim)$UMAP
  plotcol_2<-plotcol_2[rownames(coor),]
  
  plot(reducedDims(sim)$UMAP, col = plotcol_2$col, pch=16, asp = 1,cex=0.6,main=dataSlice_TAM[i])
  lines(SlingshotDataSet(sim), lwd=2, col=brewer.pal(length(SlingshotDataSet(sim)@lineages),"Set1"))
  legend("left",
         legend = plot_lab,
         col = unique(brewer.pal(length(SlingshotDataSet(sim)@lineages),"Set1")),
         inset=0.8,
         pch = 16)
  
  print(dataSlice_TAM[i])
}
dev.off()


#########################################################################################################################
#####按切片挑选对应的上升或下降的基因
dir_sling<-'/F6/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))


all_assoGene<-c()
for(i in 1:length(dataSlice)){
  #i=3
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,file_CAF_fitGAM[i]))
  startRes<-read.delim(paste0(dir_sling,file_CAF_startRes[i]))
  #startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='TAM') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')
  
  slice_assoGene<-data.frame(slice=dataSlice[i],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  
  all_assoGene<-rbind(all_assoGene,slice_assoGene)
  print(dataSlice[i])
}


write.table(all_assoGene,'/F6/TAM_assoGene.txt',quote = F,sep = '\t',row.names = F)

dir_sling<-'/F6/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))


all_assoGene<-c()
for(i in 1:length(dataSlice)){
  #i=3
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,file_CAF_fitGAM[i]))
  startRes<-read.delim(paste0(dir_sling,file_CAF_startRes[i]))
  #startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')
  
  slice_assoGene<-data.frame(slice=dataSlice[i],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  
  all_assoGene<-rbind(all_assoGene,slice_assoGene)
  print(dataSlice[i])
}


write.table(all_assoGene,'/F6/CAF_assoGene.txt',quote = F,sep = '\t',row.names = F)

########################################################################################################
###上下调基因判断
dir_out<-'/F6/'###FC>1.5
#dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot2/'###FC>1

CAF_assoGene<-read.delim('/F6/CAF_assoGene.txt',stringsAsFactors = F,check.names = F)
CAF_assoGene_FC<-reshape2::acast(CAF_assoGene[,c('gene','slice','logFC')],gene~slice)
CAF_assoGene_P<-reshape2::acast(CAF_assoGene[,c('gene','slice','pvalue')],gene~slice)
up<-c()
down<-c()
stable<-c()
slice_num<-c()
for(nn in 1:nrow(CAF_assoGene_FC)){#nn=1
  aa<-data.frame(FC=CAF_assoGene_FC[nn,],
                 Pvalue=CAF_assoGene_P[nn,])
  aa<-aa[!is.na(aa$FC),]
  slice_num<-c(slice_num,nrow(aa))
  up<-c(up,length(which(aa$Pvalue<0.05&aa$FC>log2(1))))
  down<-c(down,length(which(aa$Pvalue<0.05&aa$FC<(-log2(1)))))
  stable<-c(stable,length(which(aa$Pvalue>=0.05|abs(aa$FC)<log2(1))))
  # geneType = ifelse(aa$Pvalue< 0.05 & abs(aa$FC) >= log2(1.5), 
  #                                ifelse(aa$FC> log2(1.5) ,'Up','Down'),'Stable')
}
CAF_assoGeneType<-data.frame(gene=rownames(CAF_assoGene_FC),
                             up=up,
                             stable=stable,
                             down=down,
                             slice_num=slice_num)
CAF_assoGeneType$slice_num<-apply(CAF_assoGeneType[,2:4],1,sum)
write.table(CAF_assoGeneType,paste0(dir_out,'CAF_assoGeneType.txt'),quote = F,sep = '\t',row.names = F)
CAF_assoGeneType<-read.delim(paste0(dir_out,'CAF_assoGeneType.txt'),stringsAsFactors = F,check.names = F)


CAF_assoGene<-read.delim('/F6/CAF_assoGene.txt',stringsAsFactors = F,check.names = F)
cancer<-unlist(lapply(strsplit(CAF_assoGene$slice,'/'),function(x)x[1]))
cancer<-unique(cancer)

for(mm in 1:length(cancer)){#mm=1
  CAF_assoGene_cancer<-CAF_assoGene[grep(cancer[mm],CAF_assoGene$slice),]
  CAF_assoGene_FC<-reshape2::acast(CAF_assoGene_cancer[,c('gene','slice','logFC')],gene~slice)
  CAF_assoGene_P<-reshape2::acast(CAF_assoGene_cancer[,c('gene','slice','pvalue')],gene~slice)
  up<-c()
  down<-c()
  stable<-c()
  slice_num<-c()
  for(nn in 1:nrow(CAF_assoGene_FC)){#nn=1
    aa<-data.frame(FC=CAF_assoGene_FC[nn,],
                   Pvalue=CAF_assoGene_P[nn,])
    aa<-aa[!is.na(aa$FC),]
    slice_num<-c(slice_num,nrow(aa))
    up<-c(up,length(which(aa$Pvalue<0.05&aa$FC>log2(1))))
    down<-c(down,length(which(aa$Pvalue<0.05&aa$FC<(-log2(1)))))
    stable<-c(stable,length(which(aa$Pvalue>=0.05|abs(aa$FC)<log2(1))))
    # geneType = ifelse(aa$Pvalue< 0.05 & abs(aa$FC) >= log2(1.5), 
    #                                ifelse(aa$FC> log2(1.5) ,'Up','Down'),'Stable')
  }
  CAF_assoGeneType<-data.frame(gene=rownames(CAF_assoGene_FC),
                               up=up,
                               stable=stable,
                               down=down,
                               slice_num=slice_num)
  CAF_assoGeneType$slice_num<-apply(CAF_assoGeneType[,2:4],1,sum)
  write.table(CAF_assoGeneType,paste0(dir_out,cancer[mm],'_CAF_assoGeneType.txt'),quote = F,sep = '\t',row.names = F)
}


TAM_assoGene<-read.delim('/F6/TAM_assoGene.txt',stringsAsFactors = F,check.names = F)
TAM_assoGene_FC<-reshape2::acast(TAM_assoGene[,c('gene','slice','logFC')],gene~slice)
TAM_assoGene_P<-reshape2::acast(TAM_assoGene[,c('gene','slice','pvalue')],gene~slice)
up<-c()
down<-c()
stable<-c()
slice_num<-c()
for(nn in 1:nrow(TAM_assoGene_FC)){#nn=1
  aa<-data.frame(FC=TAM_assoGene_FC[nn,],
                 Pvalue=TAM_assoGene_P[nn,])
  aa<-aa[!is.na(aa$FC),]
  slice_num<-c(slice_num,nrow(aa))
  up<-c(up,length(which(aa$Pvalue<0.05&aa$FC>log2(1))))
  down<-c(down,length(which(aa$Pvalue<0.05&aa$FC<(-log2(1)))))
  stable<-c(stable,length(which(aa$Pvalue>=0.05|abs(aa$FC)<log2(1))))
  # geneType = ifelse(aa$Pvalue< 0.05 & abs(aa$FC) >= log2(1.5), 
  #                                ifelse(aa$FC> log2(1.5) ,'Up','Down'),'Stable')
}
TAM_assoGeneType<-data.frame(gene=rownames(TAM_assoGene_FC),
                             up=up,
                             stable=stable,
                             down=down,
                             slice_num=slice_num)
#TAM_assoGeneType$slice_num<-apply(TAM_assoGeneType[,2:4],1,sum)
write.table(TAM_assoGeneType,paste0(dir_out,'TAM_assoGeneType.txt'),quote = F,sep = '\t',row.names = F)

TAM_assoGene<-read.delim('/F6/TAM_assoGene.txt',stringsAsFactors = F,check.names = F)
cancer<-unlist(lapply(strsplit(TAM_assoGene$slice,'/'),function(x)x[1]))
cancer<-unique(cancer)

for(mm in 1:length(cancer)){#mm=1
  TAM_assoGene_cancer<-TAM_assoGene[grep(cancer[mm],TAM_assoGene$slice),]
  TAM_assoGene_FC<-reshape2::acast(TAM_assoGene_cancer[,c('gene','slice','logFC')],gene~slice)
  TAM_assoGene_P<-reshape2::acast(TAM_assoGene_cancer[,c('gene','slice','pvalue')],gene~slice)
  up<-c()
  down<-c()
  stable<-c()
  slice_num<-c()
  for(nn in 1:nrow(TAM_assoGene_FC)){#nn=1
    aa<-data.frame(FC=TAM_assoGene_FC[nn,],
                   Pvalue=TAM_assoGene_P[nn,])
    aa<-aa[!is.na(aa$FC),]
    slice_num<-c(slice_num,nrow(aa))
    up<-c(up,length(which(aa$Pvalue<0.05&aa$FC>log2(1))))
    down<-c(down,length(which(aa$Pvalue<0.05&aa$FC<(-log2(1)))))
    stable<-c(stable,length(which(aa$Pvalue>=0.05|abs(aa$FC)<log2(1))))
    # geneType = ifelse(aa$Pvalue< 0.05 & abs(aa$FC) >= log2(1.5), 
    #                                ifelse(aa$FC> log2(1.5) ,'Up','Down'),'Stable')
  }
  TAM_assoGeneType<-data.frame(gene=rownames(TAM_assoGene_FC),
                               up=up,
                               stable=stable,
                               down=down,
                               slice_num=slice_num)
  TAM_assoGeneType$slice_num<-apply(TAM_assoGeneType[,2:4],1,sum)
  write.table(TAM_assoGeneType,paste0(dir_out,cancer[mm],'_TAM_assoGeneType.txt'),quote = F,sep = '\t',row.names = F)
}



#####slingshot挑选的基因进行超几何功能富集
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggpubr)
library(dplyr)
library(circlize)
library(reshape2)
library(scales)
library(clusterProfiler)

#?read.gmt
dir_geneset<-'/F6/geneset/'
GOBP<-read.gmt(paste0(dir_geneset,'c5.go.bp.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
KEGG<-read.gmt(paste0(dir_geneset,'c2.cp.kegg_legacy.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
HALL<-read.gmt(paste0(dir_geneset,'h.all.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
MP17<-read.delim('/F6/geneset/intra23_inter23MP_list(20-300).txt',
                 stringsAsFactors = F,check.names = F)

MP17<-reshape2::melt(as.matrix(MP17)) %>% as.data.frame()
MP17<-MP17[,-1]
colnames(MP17)<-c('term','gene')
MP17$term<-paste0('MP_',MP17$term)

all_geneset<-rbind(GOBP,rbind(KEGG,rbind(HALL,MP17)))
all_geneset$term<-as.vector(all_geneset$term)
all_geneset$type<-unlist(lapply(strsplit(all_geneset$term,'_'),function(x)x[1]))
length(unique(all_geneset$term))
write.table(all_geneset,'/F6/geneset/GOBP_KEGG_HALL_MP17.txt',quote = F,sep = '\t',row.names = F)
all_geneset<-read.delim('/F6/geneset/GOBP_KEGG_HALL_MP17.txt',stringsAsFactors = F,check.names = F)



####阈值选的是1.5倍差异
#/F6/slingshot_stat
df <- CAF_assoGeneType
df$pattern <- with(df, {
  ifelse(up > stable & up > 1.5*down, "up_dominant",
         ifelse(down > stable & down > 1.5*up, "down_dominant",
                ifelse(stable >= up & stable >= down, "stable_dominant", "mixed")))
})


pan_caf_up <- df[df$pattern == "up_dominant",1]
pan_caf_down <- df[df$pattern == "down_dominant",1]

panCAF <- data.frame(
  "pan_caf_up","pan_caf_down"
)
panCAF[2,1] <- paste(pan_caf_up, collapse = ",")
panCAF[2,2] <- paste(pan_caf_down, collapse = ",")

panCAF <- t(panCAF)
rownames(panCAF) <- NULL
write.table(panCAF,'/F6/slingshot_stat/panCAF_greatStable&1.5differ.txt',quote = F,sep = '\t',row.names = F,col.names = F)


df <- CAF_assoGeneType
df$pattern <- with(df, {
  ifelse(up > 1.5*down, "up_dominant",
         ifelse(down > 1.5*up, "down_dominant",
                ifelse(stable >= up & stable >= down, "stable_dominant", "mixed")))
})


pan_caf_up <- df[df$pattern == "up_dominant",1]
pan_caf_down <- df[df$pattern == "down_dominant",1]

panCAF <- data.frame(
  "pan_caf_up","pan_caf_down"
)
panCAF[2,1] <- paste(pan_caf_up, collapse = ",")
panCAF[2,2] <- paste(pan_caf_down, collapse = ",")

panCAF <- t(panCAF)
rownames(panCAF) <- NULL
write.table(panCAF,'/F6/slingshot_stat/panCAF_1.5differ.txt',quote = F,sep = '\t',row.names = F,col.names = F)



df <- TAM_assoGeneType
df$pattern <- with(df, {
  ifelse(up > stable & up > 1.5*down, "up_dominant",
         ifelse(down > stable & down > 1.5*up, "down_dominant",
                ifelse(stable >= up & stable >= down, "stable_dominant", "mixed")))
})


pan_tam_up <- df[df$pattern == "up_dominant",1]
pan_tam_down <- df[df$pattern == "down_dominant",1]

panTAM <- data.frame(
  "pan_tam_up","pan_tam_down"
)
panTAM[2,1] <- paste(pan_tam_up, collapse = ",")
panTAM[2,2] <- paste(pan_tam_down, collapse = ",")

panTAM <- t(panTAM)
rownames(panTAM) <- NULL
write.table(panTAM,'/F6/slingshot_stat/panTAM_greatStable&1.5differ.txt',quote = F,sep = '\t',row.names = F,col.names = F)

df <- TAM_assoGeneType
df$pattern <- with(df, {
  ifelse(up > 1.5*down, "up_dominant",
         ifelse(down > 1.5*up, "down_dominant",
                ifelse(stable >= up & stable >= down, "stable_dominant", "mixed")))
})


pan_tam_up <- df[df$pattern == "up_dominant",1]
pan_tam_down <- df[df$pattern == "down_dominant",1]

panTAM <- data.frame(
  "pan_tam_up","pan_tam_down"
)
panTAM[2,1] <- paste(pan_tam_up, collapse = ",")
panTAM[2,2] <- paste(pan_tam_down, collapse = ",")

panTAM <- t(panTAM)
rownames(panTAM) <- NULL
write.table(panTAM,'/F6/slingshot_stat/panTAM_1.5differ.txt',quote = F,sep = '\t',row.names = F,col.names = F)


dir_gene<-'/F6/slingshot_stat/'
file_gene<-list.files(pattern = 'txt',path = dir_gene,recursive = T)
txt_name<-unlist(lapply(strsplit(file_gene,'.txt'),function(x)x[1]))
dir_out<-'/F6/slingshot_stat/FunEnrich/'

for(i in 1:length(file_gene)){
  #i=1
  slingshotGene<-read.delim(paste0(dir_gene,file_gene[i]),stringsAsFactors = F,check.names = F,header = F)
  
  enrich_re<-lapply(1:5,function(x){#x=1
    x_gene<-unlist(strsplit(slingshotGene[x,2],','))
    
    ph_re<-lapply(unique(all_geneset$term),function(y){#y=unique(all_geneset$term)[1]
      fun_gene<-all_geneset$gene[all_geneset$term%in%y]
      N<-19955###总的背景基因
      M<-length(x_gene)###待测基因/DEG
      n<-length(fun_gene)####目标基因/功能通路基因
      m<-length(intersect(x_gene,fun_gene))###交集
      enriched_P<-phyper(m-1, M, (N-M), n,lower.tail = F)###富集
      
      return(c(enriched_P,m,M,n,N))
    })
    
    ph_re<-do.call(rbind,ph_re) %>% as.data.frame()
    colnames(ph_re)<-c('Pvalue','InterNum','QueryNum','FunNum','BgNum')
    ph_re$geneID<-lapply(unique(all_geneset$term),function(y){#y=unique(all_geneset$term)[1]
      fun_gene<-all_geneset$gene[all_geneset$term%in%y]
      return(c(paste0(intersect(x_gene,fun_gene),collapse = ',')))
    }) %>% unlist()
    ph_re$Padjust<-p.adjust(ph_re$Pvalue, method = 'BH', n = length(ph_re$Pvalue))
    ph_re$Function<-unique(all_geneset$term)
    ph_re$slingshotGene<-slingshotGene$V1[x]
    return(ph_re)
  })
  enrich_re<-do.call(rbind,enrich_re) %>% as.data.frame()
  write.table(enrich_re,paste0(dir_out,txt_name[i],'_FunEnrich.txt'),quote = F,sep = '\t',row.names = F)
}




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
      
      warning(sprintf("路径格式异常: %s", path))
      return(path)
    }
  }, USE.NAMES = FALSE)
}



###染色
####20_slice
###########slice_gene_exp###########
#CAF_gene <- c("ELF3","CTSD","CD74","CXCL14","FASN","LGALS3")
#TAM_gene <- c("FTH1","APOE","SPARC","CFD","FOS","NUPR1")
#dataset_slice <- c('HB/GSE261958/GSM8155172','PRAD/GSE278936/GSM8557980','PN/GSE232766/GSM7373504')


CAF_gene <- c("TIMP1","FCGBP","CD74","SQSTM1","MMP2","TPM1","ID1",
              "SFRP2","NUPR1","CXCL14","MT-ND4","CIRBP","BGN","TM4SF1","CST3")

TAM_gene <- c("FTL","APOE","LAMB3","TSPAN1","IGHM","FOS",
              "SQOR","TMSB10","FN1","MMP2","FAM83H","IL32","PABPC1","CD24"
              )
dataset_slice <- c('OSCC/GSE220978/GSM6833484','GBM/GSE235672/GSM7507330','GC/GSE251950/GSM7990473',
                   'LIHC/lihc02/slice3','HGSC/GSE274657/GSM8454234','DSRCT/GSE263523/GSM8279108',
                   'MIBC/GSE246011/GSM7853988','PDAC/GSE254829/GSM8058244',"NPC/GSE206245/GSM6248650","HNSCC/GSE281978/GSM8633895")

dir_rds<-'/data/10X_Visium/new_st/'
#file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dir_copykat<-'/data/10X_Visium/new_st_copykat/'
#file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_copykat,recursive = T)

file_bdy<-paste0(dataset_slice,"_BdyCoreBud.txt")
#dataSlice<-dataSlice[match(select_slice,dataSlice)]
#file_RCTD<-paste0(select_slice,"_Deconvolution.txt")
file_rds <- convert_to_processed_rds(dataset_slice)



dir_exp <- "/data/10X_Visium/plot/gene_slingshot_plot/"
######CAF#######
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-dataset_slice[i]
  dataset_name<-str_split(dataset_he_slice,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[CAF_gene,]
  
  CAF_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                             LocalType=st_bdy$FinalLocalType,
                             dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                             imagerow=st_bdy$imagerow,  
                             imagecol=st_bdy$imagecol)
  CAF_gene_exp <- cbind(CAF_gene_exp,t(st_conut))
  write.table(CAF_gene_exp,paste0(dir_exp,slice_name,'_CAFgene.txt'),quote = F,sep = '\t')
  print(slice_name)
}
#####TAM######
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-dataset_slice[i]
  dataset_name<-str_split(dataset_he_slice,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[TAM_gene,]
  
  immune_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                               LocalType=st_bdy$FinalLocalType,
                               dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                               imagerow=st_bdy$imagerow,
                               imagecol=st_bdy$imagecol)
  immune_gene_exp <- cbind(immune_gene_exp,t(st_conut))
  write.table(immune_gene_exp,paste0(dir_exp,slice_name,'_TAMgene.txt'),quote = F,sep = '\t')
  print(slice_name)
}


######CAF####
dir_CAF_plot <- "/data/10X_Visium/plot/gene_slingshot_plot/CAF"
for(i in 1:length(dataset_slice)){
  #i <- 760
  #i <- 1
  dataset_he_slice<-dataset_slice[i]
  dataset_name<-str_split(dataset_he_slice,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  core_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_CAFgene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  core_gene_exp <- removeColsAllNa(core_gene_exp)
  if(dim(core_gene_exp)[2]>6){
    gene_data <- core_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    
    for(j in 1:length(gene_name)){
      gene_name0 <- gene_name[j]
      plot_data<-data.frame(imagerow=core_gene_exp$imagerow,imagecol=core_gene_exp$imagecol,core_gene_exp=core_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","core_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=core_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(25),
                                          colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                          colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                          colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = dataset_he_slice[i],x = "",y = "")
      
      pdf(paste0(dir_CAF_plot,"/",slice_name,"_CAF_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(slice_name)
}

#######TAM######## 
dir_TAM_plot <- "/data/10X_Visium/plot/gene_slingshot_plot/TAM"
for(i in 1:length(dataset_slice)){
  dataset_he_slice<-dataset_slice[i]
  dataset_name<-str_split(dataset_he_slice,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  immune_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_TAMgene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  immune_gene_exp <- removeColsAllNa(immune_gene_exp)
  if(dim(immune_gene_exp)[2]>6){
    gene_data <- immune_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    for (g in 1:length(gene_name)) {
      #g <- 1
      gene_name0 <- gene_name[g]
      plot_data<-data.frame(imagerow=immune_gene_exp$imagerow,imagecol=immune_gene_exp$imagecol,immune_gene_exp=immune_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","immune_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=immune_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(25),
                                          colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                          colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                          colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = slice_name,x = "",y = "")
      
      pdf(paste0(dir_TAM_plot,"/",slice_name,"_TAM_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(slice_name)
}  


