#####计算slingshot
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


# dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/slingshot/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/ST_expression/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/RCTD_celltype/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/copykat/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_RCTD,'_RCTD'),function(x)x[1]))
table(dataSlice==unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x)x[1])))
table(paste0(dataSlice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
dataSet<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[1]))

for(i in 4:130){
  #i=4
  if(!dir.exists(paste0(dir_out,dataSet[i]))) dir.create(paste0(dir_out,dataSet[i]))
  
  st_exp<-readRDS(paste0(dir_rds,file_rds[i]))
  RCTD_cell<-read.delim(paste0(dir_RCTD,file_RCTD[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]))
  RCTD_cell$localType<-st_bdy[rownames(RCTD_cell),'FinalLocalType']
  RCTD_cell$celltype[RCTD_cell$localType%in%c('Core','Boundary','Dispersion')]<-RCTD_cell$localType[RCTD_cell$localType%in%c('Core','Boundary','Dispersion')]
  RCTD_cell<-RCTD_cell[which(RCTD_cell$localType!='not.defined'),]
  
  tryCatch({
    if(length(which(RCTD_cell$celltype=='Core'))>10&&length(which(RCTD_cell$celltype=='CAF'))>10){
      st_exp<-st_exp[,rownames(RCTD_cell)]
      
      counts <- st_exp@assays[["SCT"]]@counts
      # 将表达矩阵转换为SingleCellExperiment对象
      sim <- SingleCellExperiment(assays = List(counts = counts))
      umap = st_exp@images[["slice1"]]@coordinates[,c('imagerow','imagecol')] %>% as.matrix()######必须变成矩阵
      colnames(umap) = c('UMAP_1', 'UMAP_2')
      class(umap)
      class(umap[1,1])
      # 将降维的结果添加到SingleCellExperiment对象中
      reducedDims(sim) = SimpleList(UMAP = umap)
      # metadata
      meta = st_exp@meta.data
      # colData(sim)相当于meta.data，但他不是data.frame格式
      # 所以需要分步赋予
      colData(sim)$orig.ident = meta$orig.ident
      colData(sim)$celltype = RCTD_cell$celltype
      
      #----------------------------------step3 使用slingshot构建分化谱系并进行拟时推断
      sim <- slingshot(sim, 
                       clusterLabels = 'celltype',  # 选择colData中细胞注释的列名
                       reducedDim = 'UMAP',  
                       start.clus= "Core",  # 选择起点
                       end.clus = 'CAF'     # 选择终点
      )
      saveRDS(sim,paste0(dir_out,dataSlice[i],'_slingshot_CAF.rds'))
      #aa<-readRDS(paste0(dir_out,dataSlice[i],'_slingshot.rds'))
      
      #-----------------------------------step5 tradeSeq 下游分析
      # Fit negative binomial model
      counts <- sim@assays@data$counts
      crv <- SlingshotDataSet(sim)
      #拟合负二项式模型 需要决定结的数量
      # set.seed(111)
      # icMat <- evaluateK(counts = counts, 
      #                    sds = crv, 
      #                    k = 3:10,    # no more than 12
      #                    nGenes = 500,
      #                    verbose = T)
      # saveRDS(icMat,paste0(dir_out,dataSlice[i],'_evaluateK.rds'))
      
      set.seed(111)
      pseudotime <- slingPseudotime(crv, na = FALSE)
      cellWeights <- slingCurveWeights(crv)
      # fit negative binomial GAM
      # 2k cells ~13 min
      # system.time()这个函数可以计算运行时间
      sce <- fitGAM(counts = counts, 
                    pseudotime = pseudotime, 
                    cellWeights = cellWeights,
                    nknots = 6, 
                    verbose = FALSE)
      saveRDS(sce,paste0(dir_out,dataSlice[i],'_fitGAM_CAF.rds'))
      #aa<-readRDS(paste0(dir_out,dataSlice[i],'_fitGAM.rds'))
      #探索基因表达与拟时序的相关性
      assoRes <- associationTest(sce)
      #head(assoRes)
      #寻找与起止点相关性最高的基因
      startRes <- startVsEndTest(sce)
      write.table(startRes,paste0(dir_out,dataSlice[i],'_startRes_CAF.txt'),quote = F,sep = '\t')
    }
  },error = function(e){
    stopMessage<-"something wrong"
  })
  
  
  print(dataSlice[i])
  
}





#####计算slingshot
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


# dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
# dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
# dir_RCTD<-'E:/Mirror/ST_analysis/data/10X Visium/RCTD_celltype/'
# dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'

dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/slingshot/'
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/ST_expression/'
dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/RCTD_celltype/'
dir_bdy<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/copykat/'

# dir_out<-'/data/zhouweiwei/spatial/data/10X_Visium/slingshot/'
# dir_rds<-'/data/zhouweiwei/spatial/data/10X_Visium/ST_expression/'
# dir_RCTD<-'/data/zhouweiwei/spatial/data/10X_Visium/RCTD_celltype/'
# dir_bdy<-'/data/zhouweiwei/spatial/data/10X_Visium/copykat/'

file_out<-list.files(pattern = '_startRes_TAM.txt',path = dir_out,recursive = T)
slice_out<-unlist(lapply(strsplit(file_out,'_star'),function(x)x[1]))

file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
file_RCTD<-list.files(pattern = 'txt',path = dir_RCTD,recursive = T)
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_RCTD,'_RCTD'),function(x)x[1]))
table(dataSlice==unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x)x[1])))
table(paste0(dataSlice,'.rds')==gsub('/ST/expression_position/','/',file_rds))
dataSet<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[1]))

for(i in 1:130){
  #i=1
  if(!dir.exists(paste0(dir_out,dataSet[i]))) dir.create(paste0(dir_out,dataSet[i]))
  
  if(!dataSlice[i]%in%slice_out){
    st_exp<-readRDS(paste0(dir_rds,file_rds[i]))
    RCTD_cell<-read.delim(paste0(dir_RCTD,file_RCTD[i]))
    st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]))
    RCTD_cell$localType<-st_bdy[rownames(RCTD_cell),'FinalLocalType']
    RCTD_cell$celltype[RCTD_cell$localType%in%c('Core','Boundary','Dispersion')]<-RCTD_cell$localType[RCTD_cell$localType%in%c('Core','Boundary','Dispersion')]
    RCTD_cell<-RCTD_cell[which(RCTD_cell$localType!='not.defined'),]
    
    tryCatch({
      if(length(which(RCTD_cell$celltype=='Core'))>10&&length(which(RCTD_cell$celltype=='TAM'))>10){
        st_exp<-st_exp[,rownames(RCTD_cell)]
        
        counts <- st_exp@assays[["SCT"]]@counts
        # 将表达矩阵转换为SingleCellExperiment对象
        sim <- SingleCellExperiment(assays = List(counts = counts))
        umap = st_exp@images[["slice1"]]@coordinates[,c('imagerow','imagecol')] %>% as.matrix()######必须变成矩阵
        colnames(umap) = c('UMAP_1', 'UMAP_2')
        # 将降维的结果添加到SingleCellExperiment对象中
        reducedDims(sim) = SimpleList(UMAP = umap)
        # metadata
        meta = st_exp@meta.data
        # colData(sim)相当于meta.data，但他不是data.frame格式
        # 所以需要分步赋予
        colData(sim)$orig.ident = meta$orig.ident
        colData(sim)$celltype = RCTD_cell$celltype
        
        #----------------------------------step3 使用slingshot构建分化谱系并进行拟时推断
        sim <- slingshot(sim, 
                         clusterLabels = 'celltype',  # 选择colData中细胞注释的列名
                         reducedDim = 'UMAP',  
                         start.clus= "Core",  # 选择起点
                         end.clus = 'TAM'     # 选择终点
        )
        saveRDS(sim,paste0(dir_out,dataSlice[i],'_slingshot_TAM.rds'))
        #aa<-readRDS(paste0(dir_out,dataSlice[i],'_slingshot.rds'))
        
        #-----------------------------------step5 tradeSeq 下游分析
        # Fit negative binomial model
        counts <- sim@assays@data$counts
        crv <- SlingshotDataSet(sim)
        #拟合负二项式模型 需要决定结的数量
        # set.seed(111)
        # icMat <- evaluateK(counts = counts, 
        #                    sds = crv, 
        #                    k = 3:10,    # no more than 12
        #                    nGenes = 500,
        #                    verbose = T)
        # saveRDS(icMat,paste0(dir_out,dataSlice[i],'_evaluateK.rds'))
        
        set.seed(111)
        pseudotime <- slingPseudotime(crv, na = FALSE)
        cellWeights <- slingCurveWeights(crv)
        # fit negative binomial GAM
        # 2k cells ~13 min
        # system.time()这个函数可以计算运行时间
        sce <- fitGAM(counts = counts, 
                      pseudotime = pseudotime, 
                      cellWeights = cellWeights,
                      nknots = 6, 
                      verbose = FALSE)
        saveRDS(sce,paste0(dir_out,dataSlice[i],'_fitGAM_TAM.rds'))
        #aa<-readRDS(paste0(dir_out,dataSlice[i],'_fitGAM.rds'))
        #探索基因表达与拟时序的相关性
        assoRes <- associationTest(sce)
        #head(assoRes)
        #寻找与起止点相关性最高的基因
        startRes <- startVsEndTest(sce)
        write.table(startRes,paste0(dir_out,dataSlice[i],'_startRes_TAM.txt'),quote = F,sep = '\t')
      }
    },error = function(e){
      stopMessage<-"something wrong"
    })
  }
  
  print(dataSlice[i])
  
}


###功能富集
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
dir_geneset<-'E:/Mirror/ST_analysis/data/geneset/'
GOBP<-read.gmt(paste0(dir_geneset,'c5.go.bp.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
KEGG<-read.gmt(paste0(dir_geneset,'c2.cp.kegg_legacy.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
HALL<-read.gmt(paste0(dir_geneset,'h.all.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
MP13<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/MP_list_intersect_initia25_intersect_cluster24.txt',
                 stringsAsFactors = F,check.names = F)
MP13<-MP13[,-1]
MP13<-reshape2::melt(as.matrix(MP13)) %>% as.data.frame()
MP13<-MP13[,-1]
colnames(MP13)<-c('term','gene')
MP13$term<-paste0('MP_',MP13$term)

all_geneset<-rbind(GOBP,rbind(KEGG,rbind(HALL,MP13)))
all_geneset$term<-as.vector(all_geneset$term)
all_geneset$type<-unlist(lapply(strsplit(all_geneset$term,'_'),function(x)x[1]))
length(unique(all_geneset$term))
write.table(all_geneset,'E:/Mirror/ST_analysis/data/pic_data/GOBP_KEGG_HALL_MP13.txt',quote = F,sep = '\t',row.names = F)
all_geneset<-read.delim('E:/Mirror/ST_analysis/data/pic_data/GOBP_KEGG_HALL_MP13.txt',stringsAsFactors = F,check.names = F)

dir_gene<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_stat/'
file_gene<-list.files(pattern = 'txt',path = dir_gene,recursive = T)
txt_name<-unlist(lapply(strsplit(file_gene,'.txt'),function(x)x[1]))
dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_stat/FunEnrich/'

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




#####对富集后的结果进行筛选
dir_enrich<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_stat/FunEnrich/'###FC>1.5
dir_enrich<-'E:/Mirror/ST_analysis/data/pic_data/slingshot2/stat/FunEnrich/'###FC>1
file_enrich<-list.files(pattern = '_FunEnrich.txt',path = dir_enrich)
enrich_name<-unlist(lapply(strsplit(file_enrich,'_FunEnrich.txt'),function(x)x[1]))

for(i in 1:length(file_enrich)){
  #i=8
  if(!dir.exists(paste0(dir_enrich,'SigRe/',enrich_name[i]))) dir.create(paste0(dir_enrich,'SigRe/',enrich_name[i]))
  enrich_re<-read.delim(paste0(dir_enrich,file_enrich[i]),stringsAsFactors = F,check.names = F)
  # enrich_re<-enrich_re[,c(9,8,7,1:6)]
  # write.table(enrich_re,paste0(dir_enrich,file_enrich[i]),quote = F,sep = '\t',row.names = F)
  unique(enrich_re$slingshotGene)
  
  for(j in unique(enrich_re$slingshotGene)){#j=unique(enrich_re$slingshotGene)[1]
    aa<-enrich_re[enrich_re$slingshotGene%in%j,]
    aa<-aa[which(aa$Padjust<0.05),]
    aa$Description<-gsub('_',' ',aa$Function)
    aa<-aa[,c(1,2,10,3:9)]
    write.table(aa,paste0(dir_enrich,'SigRe/',enrich_name[i],'/',j,'_sigResult.txt'),quote = F,sep = '\t',row.names = F)
  }
}





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

dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice_CAF<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice_CAF,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice_CAF==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'
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
  # colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
  # plotcol <- colors2[cut(sim$slingPseudotime_2, breaks=100)] # 这里我们用cut函数把 lineage2 分割成100个区间，同一区间的细胞给与同一个颜色
  # plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
  # plotcol
  # sim@colData$plotcol<-plotcol
  
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



line_col<-c('1'='#1e6cb6','2'='#b797cb','3'='#89d8b6','4'='#dab02a','5'='#ff7f85','6'='#fe8b15','7'='#e5b89e','8'='#00a257')
pdf(paste0(dir_pic,'associationGene_CAF.pdf'),width = 5,height = 4)
aa<-c()
for(i in 1:length(dataSlice)){
  #i=1
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,file_CAF_fitGAM[i]))
  startRes<-read.delim(paste0(dir_sling,file_CAF_startRes[i]))
  startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')#[1]
  
  aa<-c(aa,length(CAF_site))
  if(length(CAF_site)>=1){
    counts <- sim@assays@data$counts
    startRes<-startRes[order(startRes$waldStat,decreasing = T),]
    startRes<-startRes[which(abs(startRes[,(CAF_site+3)])>=1),]
    # 挑选相关性最强的基因，并可视化
    sigGeneStart <- rownames(startRes)[1]
    #sigGeneStart <- 'IL7R'
    col_use<-line_col[1:length(SlingshotDataSet(sim)@lineages)]
    names(col_use)<-1:length(col_use)
    # labels<-names(lineages)
    # labels[CAF_site]<-paste0(labels[CAF_site],'_CAF')
    
    pp<-plotSmoothers(sce, counts, gene = sigGeneStart,ylab=sigGeneStart,lwd=1,size=0.5,lineagesToPlot=CAF_site,plotLineages=T,
                      curvesCols=col_use)+
      scale_color_manual(values =  col_use)+
      # scale_color_manual(values =  c('#1ab4a0','#2375ae','#ea609e',
      #                                '#8792c8','#cca002'))+
      # scale_fill_manual(values =  c('1'='#d62d28','2'='#f6b86d','3'='#ee762d',
      #                                '4'='#a9d38a','Normal'='#2375ae'))+
      ggtitle(paste0(dataSlice[i],'_line',paste0(CAF_site,collapse = '_')))+
      #scale_fill_discrete(breaks=c("1", "2", "3",'4'),labels = labels)+
      guides(colour = guide_legend(override.aes = list(size=2)))
    print(pp)
  }
  print(dataSlice[i])
  
}
dev.off()
# ?plotSmoothers
# trace('plotSmoothers',edit = F)
# rownames(counts)[1]
# ?scale_fill_manual
# ?scale_fill_discrete
# ?guides
models <- sce; counts <- sim@assays@data$counts; gene = sigGeneStart

dm <- colData(models)$tradeSeq$dm # design matrix
# slingshotColData <- colData(models)$crv
# pseudotime <- slingshotColData[,grep(x = colnames(slingshotColData),
#                                      pattern = "pseudotime")]
# if (is.null(dim(pseudotime))) pseudotime <- matrix(pseudotime, ncol = 1)
nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
# betaMat <- rowData(models)$tradeSeq$beta[[1]]
# beta <- betaMat[id,]

#construct time variable based on cell assignments.
lcol <- timeAll <- rep(0, nrow(dm))
for (jj in seq_len(nCurves)) {
  for (ii in seq_len(nrow(dm))) {
    if (dm[ii, paste0("l", jj)] == 1) {
      timeAll[ii] <- dm[ii, paste0("t", jj)]
      lcol[ii] <- jj
    } else {
      next
    }
  }
}

# plot raw data
df <- data.frame("time" = timeAll,
                 "gene_count1" = unname(counts[names(models),]['YAP1',]),
                 "gene_count2" = unname(counts[names(models),]['CST2',]),
                 "lineage" = as.character(lcol))
colnames(df)[2:3]<-c('YAP1','CST2')
df <- df[df$lineage %in% 1,]
df<-data.frame(time=rep(df$time,2),
               gene=rep(c('YAP1','CST2'),each=nrow(df)),
               value=c(df$YAP1,df$CST2))
# rep(c(1,2,3),2)
# rep(c(1,2,3),each=3)

# p <- ggplot(df, aes(x = time, y = log1p(gene_count))) +
#   labs(x = xlab, y = ylab) +
#   theme_classic()+
#   geom_point(size = size, aes(col = lineage)) +
#   scale_color_viridis_d(alpha = alpha)
# p

formula <- y ~ poly(x, 4, raw = TRUE)
p <- ggplot(df, aes(x = time, y = log1p(value),color=gene)) +
  geom_point(size = size, aes(col = gene)) +
  geom_smooth(aes(fill = gene), method = "lm", formula = formula,se = F) +
  theme_classic()+
  scale_fill_manual(values = c("#E7B800","#00AFBB"))+
  scale_color_manual(values = c("#E7B800","#00AFBB"))
p
?geom_smooth






#########################################################################################################################
#####按切片挑选对应的上升或下降的基因
dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'

all_assoGene<-c()
for(i in 1:length(dataSlice)){
  #i=1
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

write.table(all_assoGene,'E:/Mirror/ST_analysis/data/pic_data/slingshot/CAF_assoGene.txt',quote = F,sep = '\t',row.names = F)
write.table(all_assoGene,'E:/Mirror/ST_analysis/data/pic_data/slingshot/TAM_assoGene.txt',quote = F,sep = '\t',row.names = F)


########################################################################################################
###上下调基因判断
dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot/'###FC>1.5
dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot2/'###FC>1

CAF_assoGene<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/CAF_assoGene.txt',stringsAsFactors = F,check.names = F)
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


CAF_assoGene<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/CAF_assoGene.txt',stringsAsFactors = F,check.names = F)
cancer<-unlist(lapply(strsplit(CAF_assoGene$slice,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

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




TAM_assoGene<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/TAM_assoGene.txt',stringsAsFactors = F,check.names = F)
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

TAM_assoGene<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/TAM_assoGene.txt',stringsAsFactors = F,check.names = F)
cancer<-unlist(lapply(strsplit(TAM_assoGene$slice,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

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



CAF_assoGeneType<-read.delim(paste0(dir_out,'CAF_assoGeneType.txt'),stringsAsFactors = F,check.names = F)
TAM_assoGeneType<-read.delim(paste0(dir_out,'TAM_assoGeneType.txt'),stringsAsFactors = F,check.names = F)

CAF_assoGene<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/CAF_assoGene.txt',stringsAsFactors = F,check.names = F)
cancer<-unlist(lapply(strsplit(CAF_assoGene$slice,'/'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))

j=1
cancer_CAF_assoGene<-CAF_assoGene[grep(cancer[j],CAF_assoGene$slice),]
cancer_CAF_assoGene<-cancer_CAF_assoGene[which(cancer_CAF_assoGene$pvalue<0.05),]
cancer_CAF_assoGene<-reshape2::acast(cancer_CAF_assoGene[,c('gene','slice','logFC')],gene~slice)


line_col<-c('1'='#1e6cb6','2'='#b797cb','3'='#89d8b6','4'='#dab02a','5'='#ff7f85','6'='#fe8b15','7'='#e5b89e','8'='#00a257')

UpDown<-apply(cancer_CAF_assoGene,1,function(x){
  up_asso<-length(which(x>0))
  down_asso<-length(which(x<0))
  return(c(up_asso,down_asso))
}) %>% t()
UpDown<-UpDown[which(apply(UpDown,1,sum)>0),] %>% as.data.frame()
colnames(UpDown)<-c('Up_asso','Down_asso')

cancer_file_CAF_slingshot<-file_CAF_slingshot[grep(cancer[j],file_CAF_slingshot)]
cancer_file_CAF_fitGAM<-file_CAF_fitGAM[grep(cancer[j],file_CAF_fitGAM)]
cancer_file_CAF_startRes<-file_CAF_startRes[grep(cancer[j],file_CAF_startRes)]
cancer_dataSlice<-dataSlice[grep(cancer[j],dataSlice)]

up_gene<-rownames(UpDown)[order(UpDown$Up_asso,decreasing = T)][1:2]
down_gene<-rownames(UpDown)[order(UpDown$Down_asso,decreasing = T)][1:2]

pdf(paste0(dir_pic,cancer[j],'_associationGene_CAF.pdf'),width = 5,height = 4)
for(i in 1:length(cancer_dataSlice)){
  #i=1
  sim<-readRDS(paste0(dir_sling,cancer_file_CAF_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,cancer_file_CAF_fitGAM[i]))
  #table(sce@colData@rownames==sim@colData@rownames)
  #table(sce@colData@listData[["crv"]]@rownames==sce@colData@rownames)
  startRes<-read.delim(paste0(dir_sling,cancer_file_CAF_startRes[i]))
  startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')#[1]
  
  if(length(CAF_site)>0&&length(intersect(c(up_gene,down_gene),rownames(startRes)))>0){
    models <- sce; counts <- sim@assays@data$counts
    
    up_gene<-up_gene[up_gene%in%rownames(startRes)]
    down_gene<-down_gene[down_gene%in%rownames(startRes)]
    gene_count<-data.frame(lab=1:ncol(counts))
    if(length(up_gene)>=1){
      gene_count_up<-unname(counts[names(models),][up_gene,]) %>% as.matrix() %>% as.data.frame()
      if(nrow(gene_count_up)==2) gene_count_up<-t(gene_count_up) %>% as.data.frame()
      colnames(gene_count_up)<-paste0('Up_',up_gene)
      gene_count<-cbind(gene_count,gene_count_up)
    }
    if(length(down_gene)>=1){
      gene_count_down<-unname(counts[names(models),][down_gene,]) %>% as.matrix() %>% as.data.frame()
      if(nrow(gene_count_down)==2) gene_count_down<-t(gene_count_down) %>% as.data.frame()
      colnames(gene_count_down)<-paste0('Down_',down_gene)
      gene_count<-cbind(gene_count,gene_count_down)
    }
    
    if(ncol(gene_count)>1){
      dm <- colData(models)$tradeSeq$dm %>% as.data.frame() # design matrix
      dm$spot<-as.vector(models@colData@rownames)
      
      nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
      
      lcol <- timeAll <- rep(0, nrow(dm))
      for (jj in seq_len(nCurves)) {
        for (ii in seq_len(nrow(dm))) {
          if (dm[ii, paste0("l", jj)] == 1) {
            timeAll[ii] <- dm[ii, paste0("t", jj)]
            lcol[ii] <- jj
          } else {
            next
          }
        }
      }
      
      df <- data.frame("time" = timeAll,
                       "lineage" = as.character(lcol))
      df<-cbind(df,gene_count) %>% as.data.frame()
      df<-df[,-3]
      df <- df[df$lineage %in% CAF_site,]
      df<-df[,-2]
      gene_value<-reshape2::melt(as.matrix(df[,2:ncol(df)]))
      df<-data.frame(time=rep(df$time,2),
                     gene=gene_value$Var2,
                     value=gene_value$value)
      
      formula <- y ~ poly(x, 4, raw = TRUE)
      p <- ggplot(df, aes(x = time, y = log1p(value),color=gene)) +
        #geom_point(size = 0.6, aes(col = gene)) +
        geom_smooth(aes(fill = gene), method = "lm", formula = formula,se = F) +
        theme_classic()+
        scale_fill_manual(values = c('#1e6cb6','#b797cb','#89d8b6','#dab02a'))+
        scale_color_manual(values = c('#1e6cb6','#b797cb','#89d8b6','#dab02a'))+
        ggtitle(cancer_dataSlice[i])
      p
    }
    
  }
  print(cancer_dataSlice[i])
  
}
dev.off()





dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_TAM_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_TAM_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_TAM_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice_TAM<-unlist(lapply(strsplit(file_TAM_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice_TAM,unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
table(dataSlice_TAM==unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'
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
  # colors2 <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100) # 我们把这些颜色变成100个梯度，模拟渐变色
  # plotcol <- colors2[cut(sim$slingPseudotime_2, breaks=100)] # 这里我们用cut函数把 lineage2 分割成100个区间，同一区间的细胞给与同一个颜色
  # plotcol[is.na(plotcol)] <- "lightgrey" # 不属于 lineage3 的点为NA，我们把他们变成灰色
  # plotcol
  # sim@colData$plotcol<-plotcol
  
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





####slingshot
###选择感兴趣基因或通路进行可视化
####slingshot相关结果可视化
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

line_col<-c('1'='#1e6cb6','2'='#b797cb','3'='#89d8b6','4'='#dab02a','5'='#ff7f85','6'='#fe8b15','7'='#e5b89e','8'='#00a257')
line_col<-c('#1e6cb6','#b797cb','#89d8b6','#dab02a','#ff7f85','#AB231C','#546748','#fe8b15','#00a257','#b05a28')

CAF_assoGeneType<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/CAF_assoGeneType.txt',stringsAsFactors = F,check.names = F)
TAM_assoGeneType<-read.delim('E:/Mirror/ST_analysis/data/pic_data/slingshot/TAM_assoGeneType.txt',stringsAsFactors = F,check.names = F)


interest_gene_CAF_up<-c('ACTA2','TGFB1','POSTN')
interest_gene_TAM_up<-c('CD163','IL10','PD-L1','CD274')
interest_gene_CAFTAM_up<-c('MMP2','MMP9')

interest_gene_CAF_down<-c('DCN')
interest_gene_TAM_down<-c('COL4A1','LAMB1')
interest_gene_CAFTAM_down<-c('AEBP1','MGP','SCD','COL3A1','COL1A1')

interest_gene<-data.frame(gene=c(interest_gene_CAF_up,interest_gene_TAM_up,interest_gene_TAM_down,interest_gene_CAFTAM_up),
                          type=c(rep('CAF_up',3),rep('TAM_up',4),rep('TAM_down',2),rep('CAFTAM_up',2)))

interest_gene<-data.frame(gene=c(interest_gene_CAF_down,interest_gene_TAM_down,interest_gene_CAFTAM_down),
                          type=c(rep('CAF_down',1),rep('TAM_down',2),rep('CAFTAM_down',5)))

dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'
# j=27
# # cancer_CAF_assoGene<-CAF_assoGene[grep(cancer[j],CAF_assoGene$slice),]
# # cancer_CAF_assoGene<-cancer_CAF_assoGene[which(cancer_CAF_assoGene$pvalue<0.05),]
# # cancer_CAF_assoGene<-reshape2::acast(cancer_CAF_assoGene[,c('gene','slice','logFC')],gene~slice)
# 
# cancer_file_CAF_slingshot<-file_CAF_slingshot[grep(cancer[j],file_CAF_slingshot)]
# cancer_file_CAF_fitGAM<-file_CAF_fitGAM[grep(cancer[j],file_CAF_fitGAM)]
# cancer_file_CAF_startRes<-file_CAF_startRes[grep(cancer[j],file_CAF_startRes)]
# cancer_dataSlice<-dataSlice[grep(cancer[j],dataSlice)]
dataSlice_select<-c('gist01/slice1','hn-as02/slice3','lihc03/slice6','pdac03/slice2','skcm12/slice1')
select_site<-lapply(dataSlice_select,function(x)grep(x,dataSlice)) %>% unlist()

pdf(paste0(dir_pic,'interestGene_CAFdown.pdf'),width = 5,height = 4)
for(i in 1:length(dataSlice)){
  #i=19
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,file_CAF_fitGAM[i]))
  #table(sce@colData@rownames==sim@colData@rownames)
  #table(sce@colData@listData[["crv"]]@rownames==sce@colData@rownames)
  startRes<-read.delim(paste0(dir_sling,file_CAF_startRes[i]))
  startRes<-startRes[which(startRes$pvalue<0.05),]
  
  inter_gene<-intersect(interest_gene$gene,rownames(startRes))
  if(length(inter_gene)>1){
    startRes_gene<-startRes[inter_gene,]
    interest_gene_ex<-interest_gene[match(rownames(startRes_gene),interest_gene$gene),]
    
    lineages<-SlingshotDataSet(sim)@lineages
    CAF_site<-lapply(lineages,function(x){
      rr<-'N'
      if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
      return(rr)
    }) %>% unlist()
    CAF_site<-which(CAF_site=='Y')#[1]
    
    if(length(CAF_site)==1){
      models <- sce; counts <- sim@assays@data$counts
      
      gene_count<-counts[names(models),][rownames(startRes_gene),] %>% as.matrix() %>% t() %>% as.data.frame()
      colnames(gene_count)<-paste0(interest_gene_ex$type,'_',interest_gene_ex$gene)
      
      dm <- colData(models)$tradeSeq$dm %>% as.data.frame() # design matrix
      dm$spot<-as.vector(models@colData@rownames)
      gene_count<-gene_count[dm$spot,]
      
      nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
      
      lcol <- timeAll <- rep(0, nrow(dm))
      for (jj in seq_len(nCurves)) {
        for (ii in seq_len(nrow(dm))) {
          if (dm[ii, paste0("l", jj)] == 1) {
            timeAll[ii] <- dm[ii, paste0("t", jj)]
            lcol[ii] <- jj
          } else {
            next
          }
        }
      }
      
      df <- data.frame("time" = timeAll,
                       "lineage" = as.character(lcol))
      df<-cbind(df,gene_count) %>% as.data.frame()
      df <- df[df$lineage %in% CAF_site,]
      df<-df[,-2]
      gene_value<-reshape2::melt(as.matrix(df[,2:ncol(df)]))
      df<-data.frame(time=rep(df$time,(ncol(df)-1)),
                     gene=gene_value$Var2,
                     value=gene_value$value)
      
      # formula <- y ~ poly(x, 4, raw = TRUE)
      p <- ggplot(df, aes(x = time, y = log1p(value),color=gene)) +
        #geom_point(size = 0.6, aes(col = gene)) +
        #geom_smooth(aes(fill = gene), method = "lm", formula = formula,se = T) +
        geom_smooth(aes(fill = gene),size = 1.5, span = 0.2, method = "loess", formula = y ~ x,se = F)+
        theme_classic()+
        scale_fill_manual(values = line_col[1:length(unique(df$gene))])+
        scale_color_manual(values = line_col[1:length(unique(df$gene))])+
        ggtitle(paste0(dataSlice[i],'_Core_CAF'))
      print(p)
      
    }
    
  }
  
  print(dataSlice[i])
}
dev.off()


# smooth_method = "loess"
# smooth_span = 0.2
# smooth_se = TRUE
# clrp = "milo"
# 
# ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = time,
#                                                          y = log1p(value),
#                                                          color = gene)) +
#   #trajectory_part_add_on +
#   ggplot2::geom_smooth(size = 1.5, span = smooth_span, method = smooth_method, formula = y ~ x,se = smooth_se) +
#   ggplot2::scale_y_continuous(breaks = base::seq(0 , 1, 0.2), labels = base::seq(0 , 1, 0.2)) +
#   confuns::scale_color_add_on(variable = "discrete", clrp = clrp) +
#   ggplot2::theme_classic() +
#   ggplot2::theme(
#     axis.text.x = ggplot2::element_blank(),
#     axis.ticks.x = ggplot2::element_blank(),
#     axis.line.x = ggplot2::element_line(arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "inches"))),
#     axis.line.y = ggplot2::element_line()
#   ) +
#   ggplot2::labs(x = "Trajectory Direction", y = y_title, color = "Genes")



################################################################################################################
###TAM
dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_TAM_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_TAM_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_TAM_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_TAM_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))

dataSlice_select<-c('brca08/slice2','crc05/slice2','lihc02/slice14','lihc03/slice8','pdac03/slice2')
select_site<-lapply(dataSlice_select,function(x)grep(x,dataSlice)) %>% unlist()

pdf(paste0(dir_pic,'interestGene_TAMdown.pdf'),width = 5,height = 4)
for(i in 1:length(dataSlice)){
  #i=51
  sim<-readRDS(paste0(dir_sling,file_TAM_slingshot[i]))
  sce<-readRDS(paste0(dir_sling,file_TAM_fitGAM[i]))
  #table(sce@colData@rownames==sim@colData@rownames)
  #table(sce@colData@listData[["crv"]]@rownames==sce@colData@rownames)
  startRes<-read.delim(paste0(dir_sling,file_TAM_startRes[i]))
  startRes<-startRes[which(startRes$pvalue<0.05),]
  
  inter_gene<-intersect(interest_gene$gene,rownames(startRes))
  if(length(inter_gene)>1){
    startRes_gene<-startRes[inter_gene,]
    interest_gene_ex<-interest_gene[match(rownames(startRes_gene),interest_gene$gene),]
    
    lineages<-SlingshotDataSet(sim)@lineages
    TAM_site<-lapply(lineages,function(x){
      rr<-'N'
      if(x[1]=='Core'&&x[length(x)]=='TAM') rr<-'Y'
      return(rr)
    }) %>% unlist()
    TAM_site<-which(TAM_site=='Y')#[1]
    
    if(length(TAM_site)==1){
      models <- sce; counts <- sim@assays@data$counts
      
      gene_count<-counts[names(models),][rownames(startRes_gene),] %>% as.matrix() %>% t() %>% as.data.frame()
      colnames(gene_count)<-paste0(interest_gene_ex$type,'_',interest_gene_ex$gene)
      
      dm <- colData(models)$tradeSeq$dm %>% as.data.frame() # design matrix
      dm$spot<-as.vector(models@colData@rownames)
      gene_count<-gene_count[dm$spot,]
      
      nCurves <- length(grep(x = colnames(dm), pattern = "t[1-9]"))
      
      lcol <- timeAll <- rep(0, nrow(dm))
      for (jj in seq_len(nCurves)) {
        for (ii in seq_len(nrow(dm))) {
          if (dm[ii, paste0("l", jj)] == 1) {
            timeAll[ii] <- dm[ii, paste0("t", jj)]
            lcol[ii] <- jj
          } else {
            next
          }
        }
      }
      
      df <- data.frame("time" = timeAll,
                       "lineage" = as.character(lcol))
      df<-cbind(df,gene_count) %>% as.data.frame()
      df <- df[df$lineage %in% TAM_site,]
      df<-df[,-2]
      gene_value<-reshape2::melt(as.matrix(df[,2:ncol(df)]))
      df<-data.frame(time=rep(df$time,(ncol(df)-1)),
                     gene=gene_value$Var2,
                     value=gene_value$value)
      
      # formula <- y ~ poly(x, 4, raw = TRUE)
      p <- ggplot(df, aes(x = time, y = log1p(value),color=gene)) +
        #geom_point(size = 0.6, aes(col = gene)) +
        #geom_smooth(aes(fill = gene), method = "lm", formula = formula,se = T) +
        geom_smooth(aes(fill = gene),size = 1.5, span = 0.2, method = "loess", formula = y ~ x,se = F)+
        theme_classic()+
        scale_fill_manual(values = line_col[1:length(unique(df$gene))])+
        scale_color_manual(values = line_col[1:length(unique(df$gene))])+
        ggtitle(paste0(dataSlice[i],'_Core_TAM'))
      print(p)
      
    }
    
  }
  
  print(dataSlice[i])
}
dev.off()


















