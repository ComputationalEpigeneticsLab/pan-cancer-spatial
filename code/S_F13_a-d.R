####根据轨迹变化对基因进行排序并进行GO富集
##对富集结果进行网络可视化
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
library(DOSE)
library(clusterProfiler)
library(dplyr)
library(GOSemSim)
library(ggridges)
library(patchwork)
library(org.Hs.eg.db)
??get_GO_data

# 加载GO语义数据
# hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)
source("E:/Mirror/ST_analysis/program/getGoTerm.R")
GO_DATA <- get_GO_data("org.Hs.eg.db", "BP", "SYMBOL")
class(GO_DATA)
table(GO_DATA$GO2ONT)

dir_geneset<-'E:/Mirror/ST_analysis/data/geneset/'
KEGG<-read.gmt(paste0(dir_geneset,'c2.cp.kegg_legacy.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
HALL<-read.gmt(paste0(dir_geneset,'h.all.v2024.1.Hs.symbols.txt')) %>% as.data.frame()
MP13<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/MP_list_intersect_initia25_intersect_cluster24.txt',
                 stringsAsFactors = F,check.names = F)
MP13<-MP13[,-1]
colnames(MP13)<-paste0('MP_',1:13)
MP13<-reshape2::melt(as.matrix(MP13)) %>% as.data.frame()
MP13<-MP13[,-1]
colnames(MP13)<-c('term','gene')
MP13$term<-paste0('MP_',MP13$term)
add_geneset<-rbind(MP13,KEGG,HALL)
table(add_geneset$term)
add_geneset$type<-unlist(lapply(strsplit(add_geneset$term,'_'),function(x)x[1]))

EXTID2PATHID<-GO_DATA$EXTID2PATHID
EXTID2PATHID_add<-EXTID2PATHID
addGene<-setdiff(add_geneset$gene,names(EXTID2PATHID))
addGene_set<-lapply(addGene,function(x){
  aa<-add_geneset[add_geneset$gene%in%x,]
  return(c(aa$term))
})
names(addGene_set)<-addGene
EXTID2PATHID_add<-lapply(names(EXTID2PATHID),function(x){#x=names(EXTID2PATHID)[1]
  aa<-add_geneset[add_geneset$gene%in%x,]
  return(c(EXTID2PATHID[[x]],aa$term))
})
names(EXTID2PATHID_add)<-names(EXTID2PATHID)
EXTID2PATHID_add<-c(EXTID2PATHID_add,addGene_set)

GO2ONT<-GO_DATA$GO2ONT
GO2ONT_add<-unique(add_geneset$term)
names(GO2ONT_add)<-GO2ONT_add
GO2ONT_add<-unlist(lapply(strsplit(GO2ONT_add,'_'),function(x)x[1]))
GO2ONT_add<-c(GO2ONT,GO2ONT_add)

PATHID2EXTID<-GO_DATA$PATHID2EXTID
addSet<-lapply(unique(add_geneset$term),function(x){
  add_geneset[add_geneset$term%in%x,'gene']
})
names(addSet)<-unique(add_geneset$term)
PATHID2EXTID_add<-c(PATHID2EXTID,addSet)

PATHID2NAME<-GO_DATA$PATHID2NAME
PATHID2NAME_add<-unique(add_geneset$term)
names(PATHID2NAME_add)<-PATHID2NAME_add
PATHID2NAME_add<-c(PATHID2NAME,PATHID2NAME_add)

GO_DATA$EXTID2PATHID<-EXTID2PATHID_add
GO_DATA$GO2ONT<-GO2ONT_add
GO_DATA$PATHID2EXTID<-PATHID2EXTID_add
GO_DATA$PATHID2NAME<-PATHID2NAME_add

saveRDS(GO_DATA,'E:/Mirror/ST_analysis/program/GO_DATA_add.rds')

##使用调整后的基因集进行GSEA
GO_DATA<-readRDS('E:/Mirror/ST_analysis/program/GO_DATA_add.rds')
trace('gseGO',edit = T)

#########################################################################################################################
#####按切片挑选对应的上升或下降的基因
dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_CAF_slingshot<-list.files(pattern = 'slingshot_CAF.rds',path = dir_sling,recursive = T)
file_CAF_startRes<-list.files(pattern = '_startRes_CAF.txt',path = dir_sling,recursive = T)
file_CAF_fitGAM<-list.files(pattern = 'fitGAM_CAF',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_CAF_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_CAF_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'
dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_GO/CAF/'

dataSlice_select<-c('gist01/slice1','hn-as02/slice3','lihc03/slice6','pdac03/slice2','skcm12/slice1')
dataSlice_select<-c('brca08/slice2','brca26/slice1','hn-as01/slice5','lihc03/slice6','lihc03/slice7','skcm13/slice1')
select_site<-lapply(dataSlice_select,function(x)grep(x,dataSlice)) %>% unlist()


for(ss in select_site){ ##1:length(dataSlice)
  #ss=56
  sim<-readRDS(paste0(dir_sling,file_CAF_slingshot[ss]))
  sce<-readRDS(paste0(dir_sling,file_CAF_fitGAM[ss]))
  startRes<-read.delim(paste0(dir_sling,file_CAF_startRes[ss]))
  #startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  CAF_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='CAF') rr<-'Y'
    return(rr)
  }) %>% unlist()
  CAF_site<-which(CAF_site=='Y')
  
  models <- sce; counts <- sim@assays@data$counts
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
  df$spot<-dm$spot
  df <- df[df$lineage %in% CAF_site,]
  counts<-counts[,df$spot]
  
  # cor_time<-apply(counts,1,function(x){cor.test(df$time,x)[["estimate"]][["cor"]]})
  # cor_time<-0-cor_time
  # ?cor.test
  
  # slice_assoGene<-data.frame(slice=dataSlice[ss],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  # geneList<-slice_assoGene$logFC
  # names(geneList)<-slice_assoGene$gene
  # gene_list_sorted <- sort(cor_time, decreasing = TRUE)
  # gsea_res <- clusterProfiler::gseGO(
  #   geneList = gene_list_sorted,
  #   ont = 'BP', # 可选 "BP" (生物过程), "MF" (分子功能), "CC" (细胞组分)
  #   OrgDb = org.Hs.eg.db,
  #   keyType = "SYMBOL",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "BH", # 多重检验校正方法
  #   minGSSize = 5,
  #   maxGSSize = 5000,
  #   verbose = FALSE,
  #   seed = 123
  # )
  # # ?clusterProfiler::gseGO
  # saveRDS(gsea_res,paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GSEAres.rds"))
  gsea_res<-readRDS(paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GSEAres.rds"))
  goresult <- data.frame(gsea_res@result)
  
  ###去除冗余通路，保留最具代表性的通路
  sim_res <- clusterProfiler::simplify(gsea_res, 
                                       cutoff = 0.8, # 相似性阈值（0.6-0.8）
                                       by = "p.adjust", 
                                       measure = "Wang", # 语义相似度算法
                                       semData = NULL) # 自动加载GO.db
  
  # 查看简化后通路数量
  nrow(sim_res@result)
  gosim_res <- data.frame(sim_res@result)
  
  # head(gosim_res,20)
  if(nrow(gosim_res)>5&length(which(gosim_res$NES>0))>5&length(which(gosim_res$NES<0))>5){
    ####构建通路网络
    ##方法1：基于基因重叠（Jaccard相似性）
    # 获取所有通路的基因列表
    
    ##up
    gosim_res_up<-gosim_res[which(gosim_res$NES>0),]
    pathway_genes <- lapply(gosim_res_up$ID, function(go_id){
      genes <- strsplit(gosim_res_up[gosim_res_up$ID == go_id, "core_enrichment"], "/")[[1]]
      unique(genes)
    })
    names(pathway_genes) <- gosim_res_up$ID
    
    # 计算Jaccard相似性系数
    jaccard_matrix <- matrix(nrow = length(pathway_genes), ncol = length(pathway_genes))
    rownames(jaccard_matrix) <- names(pathway_genes)
    colnames(jaccard_matrix) <- names(pathway_genes)
    
    for(i in 1:(length(pathway_genes)-1)){
      for(j in (i+1):length(pathway_genes)){
        a <- length(intersect(pathway_genes[[i]], pathway_genes[[j]]))
        b <- length(union(pathway_genes[[i]], pathway_genes[[j]]))
        jaccard_matrix[i,j] <- a/b
      }
    }
    
    # 转换为边列表
    edges <- as.data.frame(which(jaccard_matrix > 0.1 & upper.tri(jaccard_matrix), arr.ind = TRUE)) %>% # 设置Jaccard阈值
      mutate(
        from = rownames(jaccard_matrix)[row],
        to = colnames(jaccard_matrix)[col],
        weight = jaccard_matrix[cbind(row, col)]
      ) %>%
      # select(from, to, weight) %>%
      filter(weight > 0.05) # 过滤弱连接
    edges$from_Description<-gosim_res_up$Description[match(edges$from,gosim_res_up$ID)]
    edges$to_Description<-gosim_res_up$Description[match(edges$to,gosim_res_up$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GO_network_edges_up.txt"), sep = "\t", row.names = F,quote = F)
    
    
    ###down
    gosim_res_down<-gosim_res[which(gosim_res$NES<0),]
    pathway_genes <- lapply(gosim_res_down$ID, function(go_id){
      genes <- strsplit(gosim_res_down[gosim_res_down$ID == go_id, "core_enrichment"], "/")[[1]]
      unique(genes)
    })
    names(pathway_genes) <- gosim_res_down$ID
    
    # 计算Jaccard相似性系数
    jaccard_matrix <- matrix(nrow = length(pathway_genes), ncol = length(pathway_genes))
    rownames(jaccard_matrix) <- names(pathway_genes)
    colnames(jaccard_matrix) <- names(pathway_genes)
    
    for(i in 1:(length(pathway_genes)-1)){
      for(j in (i+1):length(pathway_genes)){
        a <- length(intersect(pathway_genes[[i]], pathway_genes[[j]]))
        b <- length(union(pathway_genes[[i]], pathway_genes[[j]]))
        jaccard_matrix[i,j] <- a/b
      }
    }
    
    # 转换为边列表
    edges <- as.data.frame(which(jaccard_matrix > 0.1 & upper.tri(jaccard_matrix), arr.ind = TRUE)) %>% # 设置Jaccard阈值
      mutate(
        from = rownames(jaccard_matrix)[row],
        to = colnames(jaccard_matrix)[col],
        weight = jaccard_matrix[cbind(row, col)]
      ) %>%
      # select(from, to, weight) %>%
      filter(weight > 0.05) # 过滤弱连接
    edges$from_Description<-gosim_res_down$Description[match(edges$from,gosim_res_down$ID)]
    edges$to_Description<-gosim_res_down$Description[match(edges$to,gosim_res_down$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GO_network_edges_down.txt"), sep = "\t", row.names = F,quote = F)
    # ####方法2：基于语义相似性（GOSemSim）
    # # 计算语义相似性矩阵（使用Wang方法）
    # sem_sim <- mgoSim(sim_res@result$ID, 
    #                   sim_res@result$ID,
    #                   semData=hsGO,
    #                   measure="Wang",
    #                   combine=NULL)
    # # 转换为长格式
    # edges_sem <- as.data.frame(as.table(sem_sim)) %>%
    #   filter(Var1 != Var2) %>% # 去除对角线
    #   rename(from=Var1, to=Var2, weight=Freq) %>%
    #   filter(weight > 0.3) # 设置相似性阈值
    
    # 绘图
    # data <- data.frame(
    #   Symbol = names(gene_list_sorted),
    #   cor_R = gene_list_sorted)
    # 
    # gosim_res_p <- gosim_res %>% 
    #   dplyr::arrange(qvalue) %>%
    #   head(20) %>%
    #   dplyr::mutate(log10P = -log10(pvalue)) %>%
    #   separate_rows(core_enrichment, sep = "/") %>% 
    #   left_join(data, by = c("core_enrichment" = "Symbol")) %>%
    #   dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
    # 
    # custom_colors <- colorRampPalette(c("#8075ad", "#f5edf0","#f5da73", "#ffb800"))(100)
    # # 使用 iris 数据集，根据不同物种绘制峰峦图
    # p <- ggplot(gosim_res_p, aes(x = cor_R, y = Description, fill = log10P)) +
    #   geom_density_ridges(alpha = 0.7, scale = 1.5) +
    #   labs(x = "cor_R", y = "") +
    #   scale_fill_gradientn(colors = custom_colors,name = "-Log10(pvalue)") +
    #   theme(
    #     panel.background = element_blank(),
    #     panel.grid = element_blank(),
    #     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    #     axis.text = element_text(color = "black", size = 12),
    #     axis.title = element_text(size = 12),
    #     legend.text = element_text(size = 12),
    #     legend.title = element_text(size = 12)
    #   )
    # # p
    # 
    # gosim_res_pp <- unique(gosim_res_p[,c(2,5)])
    # gosim_res_pp$order<-factor(paste0('A',1:nrow(gosim_res_pp)),levels = rev(paste0('A',1:nrow(gosim_res_pp))))
    # custom_colors1 <- colorRampPalette(c("#008bd0", "#eeeeee", "#fc4e00"))(100)
    # p1 <- ggplot(gosim_res_pp, aes(x = 1, y = order)) +
    #   geom_point(aes(size = abs(NES), color = NES)) + # 气泡大小和颜色映射 NES 数据
    #   # xlim(c(0.98,1.02))+
    #   scale_color_gradientn(colors = custom_colors1, name = "NES", 
    #                         limits = c(-max(abs(gosim_res_pp$NES)), max(abs(gosim_res_pp$NES)))) + # 颜色渐变
    #   scale_size_continuous(range = c(2, 6)) + # 气泡大小范围
    #   labs(x = NULL, y = NULL) + # 去除 x 和 y 轴标签
    #   theme_minimal() + # 简洁主题
    #   ggtitle(dataSlice[ss])+
    #   theme(
    #     panel.background = element_blank(), # 无背景
    #     panel.grid = element_blank(),
    #     panel.border = element_blank(), # 无边框
    #     axis.ticks = element_blank(), # 去除刻度
    #     axis.text.x = element_blank() # 去除 x 轴标签
    #   )
    # # p1
    # pdf(paste0('E:/Mirror/ST_analysis/pic/slingshot/GOenrich/CAF/',gsub('/','_',dataSlice[ss]),'_GOenrich.pdf'),width = 14,height = 6)
    # print(p+p1)
    # dev.off()
  }
  
  print(dataSlice[ss])
}


###CAF up和down合并
dir_GO<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_GO/TAM/'
file_updown<-list.files(pattern = 'network_edges_up.txt',path = dir_GO)

all_net_data<-c()
for(i in 1:length(file_updown)){
  #i=1
  net_data<-read.delim(paste0(dir_GO,file_updown[i]),stringsAsFactors = F,check.names = F)
  all_net_data<-rbind(all_net_data,net_data)
}
all_net_data<-all_net_data[order(all_net_data$weight,decreasing = T),]
table(duplicated(all_net_data[,1:2]))
all_net_data<-all_net_data[!duplicated(all_net_data[,1:2]),]
write.table(all_net_data,paste0(dir_GO,'all_edges_up.txt'),quote = F,sep = '\t',row.names = F)

name_path<-data.frame(path=c(all_net_data$from,all_net_data$to),
                      name=c(all_net_data$from_Description,all_net_data$to_Description))
name_path<-name_path[!duplicated(name_path$path),]
name_path$name<-gsub('_',' ',name_path$name)
write.table(name_path,paste0(dir_GO,'all_name_path_up.txt'),quote = F,sep = '\t',row.names = F)

dir_GO<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_GO/TAM/' ### CAF
CAF_updown<-read.delim(paste0(dir_GO,'all_edges_up.txt'),stringsAsFactors = F,check.names = F)
CAF_updown_name<-read.delim(paste0(dir_GO,'all_edges_up_name.txt'),stringsAsFactors = F,check.names = F)

net_site<-lapply(1:nrow(CAF_updown),function(x){#x=1
  from_lab<-CAF_updown$from[x]%in%CAF_updown_name$names
  to_lab<-CAF_updown$to[x]%in%CAF_updown_name$names
  x_site<-'N'
  if(from_lab==T&to_lab==T) x_site<-'Y' #### 看筛选出的结果是多少 多的话用& 少用|
  return(x_site)
}) %>% unlist()
table(net_site)
CAF_updown_select<-CAF_updown[which(net_site=='Y'),]
write.table(CAF_updown_select,paste0(dir_GO,'all_edges_up_select.txt'),quote = F,sep = '\t',row.names = F)




#?simplify
#########################################################################################################################
#####按切片挑选对应的上升或下降的基因
dir_sling<-'E:/Mirror/ST_analysis/data/10X Visium/slingshot/'
file_TAM_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_TAM_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_TAM_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_TAM_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
dir_pic<-'E:/Mirror/ST_analysis/pic/slingshot/'
dir_out<-'E:/Mirror/ST_analysis/data/pic_data/slingshot_GO/TAM/'

dataSlice_select<-c('gist01/slice1','hn-as02/slice3','lihc03/slice6','pdac03/slice2','skcm12/slice1')
dataSlice_select<-c('brca08/slice2','brca26/slice1','hn-as01/slice5','lihc03/slice6','lihc03/slice7','skcm13/slice1')
select_site<-lapply(dataSlice_select,function(x)grep(x,dataSlice)) %>% unlist()


for(ss in select_site){ ##1:length(dataSlice)
  #ss=83
  sim<-readRDS(paste0(dir_sling,file_TAM_slingshot[ss]))
  sce<-readRDS(paste0(dir_sling,file_TAM_fitGAM[ss]))
  startRes<-read.delim(paste0(dir_sling,file_TAM_startRes[ss]))
  #startRes<-startRes[which(startRes$pvalue<0.05),]
  
  lineages<-SlingshotDataSet(sim)@lineages
  TAM_site<-lapply(lineages,function(x){
    rr<-'N'
    if(x[1]=='Core'&&x[length(x)]=='TAM') rr<-'Y'
    return(rr)
  }) %>% unlist()
  TAM_site<-which(TAM_site=='Y')
  
  models <- sce; counts <- sim@assays@data$counts
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
  df$spot<-dm$spot
  df <- df[df$lineage %in% TAM_site,]
  counts<-counts[,df$spot]
  
  # cor_time<-apply(counts,1,function(x){cor.test(df$time,x)[["estimate"]][["cor"]]})
  # cor_time<-0-cor_time
  
  # slice_assoGene<-data.frame(slice=dataSlice[ss],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  # geneList<-slice_assoGene$logFC
  # names(geneList)<-slice_assoGene$gene
  # gene_list_sorted <- sort(cor_time, decreasing = TRUE)
  # gsea_res <- clusterProfiler::gseGO(
  #   geneList = gene_list_sorted,
  #   ont = 'BP', # 可选 "BP" (生物过程), "MF" (分子功能), "CC" (细胞组分)
  #   OrgDb = org.Hs.eg.db,
  #   keyType = "SYMBOL",
  #   pvalueCutoff = 0.05,
  #   pAdjustMethod = "BH", # 多重检验校正方法
  #   minGSSize = 5,
  #   maxGSSize = 5000,
  #   verbose = FALSE,
  #   seed = 123
  # )
  # # ?clusterProfiler::gseGO
  # saveRDS(gsea_res,paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GSEAres.rds"))
  gsea_res<-readRDS(paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GSEAres.rds"))
  goresult <- data.frame(gsea_res@result)
  
  ###去除冗余通路，保留最具代表性的通路
  sim_res <- clusterProfiler::simplify(gsea_res, 
                                       cutoff = 0.8, # 相似性阈值（0.6-0.8）
                                       by = "p.adjust", 
                                       measure = "Wang", # 语义相似度算法
                                       semData = NULL) # 自动加载GO.db
  
  # 查看简化后通路数量
  nrow(sim_res@result)
  gosim_res <- data.frame(sim_res@result)
  
  # head(gosim_res,20)
  if(nrow(gosim_res)>5&length(which(gosim_res$NES>0))>5&length(which(gosim_res$NES<0))>5){
    ####构建通路网络
    ##方法1：基于基因重叠（Jaccard相似性）
    # 获取所有通路的基因列表
    
    ##up
    gosim_res_up<-gosim_res[which(gosim_res$NES>0),]
    pathway_genes <- lapply(gosim_res_up$ID, function(go_id){
      genes <- strsplit(gosim_res_up[gosim_res_up$ID == go_id, "core_enrichment"], "/")[[1]]
      unique(genes)
    })
    names(pathway_genes) <- gosim_res_up$ID
    
    # 计算Jaccard相似性系数
    jaccard_matrix <- matrix(nrow = length(pathway_genes), ncol = length(pathway_genes))
    rownames(jaccard_matrix) <- names(pathway_genes)
    colnames(jaccard_matrix) <- names(pathway_genes)
    
    for(i in 1:(length(pathway_genes)-1)){
      for(j in (i+1):length(pathway_genes)){
        a <- length(intersect(pathway_genes[[i]], pathway_genes[[j]]))
        b <- length(union(pathway_genes[[i]], pathway_genes[[j]]))
        jaccard_matrix[i,j] <- a/b
      }
    }
    
    # 转换为边列表
    edges <- as.data.frame(which(jaccard_matrix > 0.1 & upper.tri(jaccard_matrix), arr.ind = TRUE)) %>% # 设置Jaccard阈值
      mutate(
        from = rownames(jaccard_matrix)[row],
        to = colnames(jaccard_matrix)[col],
        weight = jaccard_matrix[cbind(row, col)]
      ) %>%
      # select(from, to, weight) %>%
      filter(weight > 0.05) # 过滤弱连接
    edges$from_Description<-gosim_res_up$Description[match(edges$from,gosim_res_up$ID)]
    edges$to_Description<-gosim_res_up$Description[match(edges$to,gosim_res_up$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GO_network_edges_up.txt"), sep = "\t", row.names = F,quote = F)
    
    
    ###down
    gosim_res_down<-gosim_res[which(gosim_res$NES<0),]
    pathway_genes <- lapply(gosim_res_down$ID, function(go_id){
      genes <- strsplit(gosim_res_down[gosim_res_down$ID == go_id, "core_enrichment"], "/")[[1]]
      unique(genes)
    })
    names(pathway_genes) <- gosim_res_down$ID
    
    # 计算Jaccard相似性系数
    jaccard_matrix <- matrix(nrow = length(pathway_genes), ncol = length(pathway_genes))
    rownames(jaccard_matrix) <- names(pathway_genes)
    colnames(jaccard_matrix) <- names(pathway_genes)
    
    for(i in 1:(length(pathway_genes)-1)){
      for(j in (i+1):length(pathway_genes)){
        a <- length(intersect(pathway_genes[[i]], pathway_genes[[j]]))
        b <- length(union(pathway_genes[[i]], pathway_genes[[j]]))
        jaccard_matrix[i,j] <- a/b
      }
    }
    
    # 转换为边列表
    edges <- as.data.frame(which(jaccard_matrix > 0.1 & upper.tri(jaccard_matrix), arr.ind = TRUE)) %>% # 设置Jaccard阈值
      mutate(
        from = rownames(jaccard_matrix)[row],
        to = colnames(jaccard_matrix)[col],
        weight = jaccard_matrix[cbind(row, col)]
      ) %>%
      # select(from, to, weight) %>%
      filter(weight > 0.05) # 过滤弱连接
    edges$from_Description<-gosim_res_down$Description[match(edges$from,gosim_res_down$ID)]
    edges$to_Description<-gosim_res_down$Description[match(edges$to,gosim_res_down$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GO_network_edges_down.txt"), sep = "\t", row.names = F,quote = F)
    # ####方法2：基于语义相似性（GOSemSim）
    # # 计算语义相似性矩阵（使用Wang方法）
    # sem_sim <- mgoSim(sim_res@result$ID, 
    #                   sim_res@result$ID,
    #                   semData=hsGO,
    #                   measure="Wang",
    #                   combine=NULL)
    # # 转换为长格式
    # edges_sem <- as.data.frame(as.table(sem_sim)) %>%
    #   filter(Var1 != Var2) %>% # 去除对角线
    #   rename(from=Var1, to=Var2, weight=Freq) %>%
    #   filter(weight > 0.3) # 设置相似性阈值
    
    # 绘图
    # data <- data.frame(
    #   Symbol = names(gene_list_sorted),
    #   cor_R = gene_list_sorted)
    # 
    # gosim_res_p <- gosim_res %>% 
    #   dplyr::arrange(qvalue) %>%
    #   head(20) %>%
    #   dplyr::mutate(log10P = -log10(pvalue)) %>%
    #   separate_rows(core_enrichment, sep = "/") %>% 
    #   left_join(data, by = c("core_enrichment" = "Symbol")) %>%
    #   dplyr::mutate(Description = factor(Description, levels = rev(unique(Description))))
    # 
    # custom_colors <- colorRampPalette(c("#8075ad", "#f5edf0","#f5da73", "#ffb800"))(100)
    # # 使用 iris 数据集，根据不同物种绘制峰峦图
    # p <- ggplot(gosim_res_p, aes(x = cor_R, y = Description, fill = log10P)) +
    #   geom_density_ridges(alpha = 0.7, scale = 1.5) +
    #   labs(x = "cor_R", y = "") +
    #   scale_fill_gradientn(colors = custom_colors,name = "-Log10(pvalue)") +
    #   theme(
    #     panel.background = element_blank(),
    #     panel.grid = element_blank(),
    #     panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    #     axis.text = element_text(color = "black", size = 12),
    #     axis.title = element_text(size = 12),
    #     legend.text = element_text(size = 12),
    #     legend.title = element_text(size = 12)
    #   )
    # # p
    # 
    # gosim_res_pp <- unique(gosim_res_p[,c(2,5)])
    # gosim_res_pp$order<-factor(paste0('A',1:nrow(gosim_res_pp)),levels = rev(paste0('A',1:nrow(gosim_res_pp))))
    # custom_colors1 <- colorRampPalette(c("#008bd0", "#eeeeee", "#fc4e00"))(100)
    # p1 <- ggplot(gosim_res_pp, aes(x = 1, y = order)) +
    #   geom_point(aes(size = abs(NES), color = NES)) + # 气泡大小和颜色映射 NES 数据
    #   # xlim(c(0.98,1.02))+
    #   scale_color_gradientn(colors = custom_colors1, name = "NES", 
    #                         limits = c(-max(abs(gosim_res_pp$NES)), max(abs(gosim_res_pp$NES)))) + # 颜色渐变
    #   scale_size_continuous(range = c(2, 6)) + # 气泡大小范围
    #   labs(x = NULL, y = NULL) + # 去除 x 和 y 轴标签
    #   theme_minimal() + # 简洁主题
    #   ggtitle(dataSlice[ss])+
    #   theme(
    #     panel.background = element_blank(), # 无背景
    #     panel.grid = element_blank(),
    #     panel.border = element_blank(), # 无边框
    #     axis.ticks = element_blank(), # 去除刻度
    #     axis.text.x = element_blank() # 去除 x 轴标签
    #   )
    # # p1
    # pdf(paste0('E:/Mirror/ST_analysis/pic/slingshot/GOenrich/CAF/',gsub('/','_',dataSlice[ss]),'_GOenrich.pdf'),width = 14,height = 6)
    # print(p+p1)
    # dev.off()
  }
  
  print(dataSlice[ss])
}







# ?gseGO
# data(geneList, package = "DOSE")
# gene_list_sorted <- sort(geneList, decreasing = TRUE)
# 
# gsea_res <- clusterProfiler::gseGO(
#   geneList = gene_list_sorted,
#   ont = "BP", # 可选 "BP" (生物过程), "MF" (分子功能), "CC" (细胞组分)
#   OrgDb = org.Hs.eg.db,
#   eps = 1e-50,
#   keyType = "ENTREZID",
#   pvalueCutoff = 0.05,
#   pAdjustMethod = "BH", # 多重检验校正方法
#   minGSSize = 5,
#   maxGSSize = 5000,
#   verbose = FALSE
# )


