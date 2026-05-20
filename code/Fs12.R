dir_out1<-'/Fs12/new/CAF/'
dir_out<-'/Fs12/CAF/'
dataSlice_select <- c('OSCC/GSE220978/GSM6833484','GBM/GSE235672/GSM7507330','GC/GSE251950/GSM7990473',
                      'LIHC/lihc02/slice3','HGSC/GSE274657/GSM8454234','DSRCT/GSE263523/GSM8279108',
                      'MIBC/GSE246011/GSM7853988','PDAC/GSE254829/GSM8058244',"NPC/GSE206245/GSM6248650",
                      "HNSCC/GSE281978/GSM8633895")

select_site<-lapply(dataSlice_select,function(x)grep(x,dataSlice)) %>% unlist()


for(ss in select_site){ ##1:length(dataSlice)
  #ss=77
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
  
  cor_time<-apply(counts,1,function(x){cor.test(df$time,x)[["estimate"]][["cor"]]})
  # cor_time<-0-cor_time
  
  # slice_assoGene<-data.frame(slice=dataSlice[ss],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  # geneList<-slice_assoGene$logFC
  # names(geneList)<-slice_assoGene$gene
  
  # ?clusterProfiler::gseGO
  #saveRDS(gsea_res,paste0(dir_out,gsub('/','_',dataSlice[ss]),"_GSEAres.rds"))
  
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
      filter(weight > 0.1) 
    
    
    edges$from_Description<-gosim_res_up$Description[match(edges$from,gosim_res_up$ID)]
    edges$to_Description<-gosim_res_up$Description[match(edges$to,gosim_res_up$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out1,gsub('/','_',dataSlice[ss]),"_GO_network_edges_up.txt"), sep = "\t", row.names = F,quote = F)
    
    
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
      filter(weight > 0.1) 
    edges$from_Description<-gosim_res_down$Description[match(edges$from,gosim_res_down$ID)]
    edges$to_Description<-gosim_res_down$Description[match(edges$to,gosim_res_down$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out1,gsub('/','_',dataSlice[ss]),"_GO_network_edges_down.txt"), sep = "\t", row.names = F,quote = F)

  }
  
  print(dataSlice[ss])
}



dir_sling<-'/F6/slingshot/'
file_TAM_slingshot<-list.files(pattern = 'slingshot_TAM.rds',path = dir_sling,recursive = T)
file_TAM_startRes<-list.files(pattern = '_startRes_TAM.txt',path = dir_sling,recursive = T)
file_TAM_fitGAM<-list.files(pattern = 'fitGAM_TAM',path = dir_sling,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_TAM_slingshot,'_slingshot'),function(x)x[1]))
setdiff(dataSlice,unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
table(dataSlice==unlist(lapply(strsplit(file_TAM_startRes,'_star'),function(x)x[1])))
dir_out1<-'/Fs12/new/TAM/'
dir_out<-'/Fs12/TAM/'

#dataSlice_select<-c('gist01/slice1','hn-as02/slice3','lihc03/slice6','pdac03/slice2','skcm12/slice1')
#dataSlice_select<-c('brca08/slice2','brca26/slice1','hn-as01/slice5','lihc03/slice6','lihc03/slice7','skcm13/slice1')
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
  cor_time<-apply(counts,1,function(x){cor.test(df$time,x)[["estimate"]][["cor"]]})
  # cor_time<-0-cor_time
  
  # slice_assoGene<-data.frame(slice=dataSlice[ss],gene=rownames(startRes),logFC=startRes[,(CAF_site+3)],pvalue=startRes$pvalue)
  # geneList<-slice_assoGene$logFC
  # names(geneList)<-slice_assoGene$gene
  
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
      filter(weight > 0.1)
    edges$from_Description<-gosim_res_up$Description[match(edges$from,gosim_res_up$ID)]
    edges$to_Description<-gosim_res_up$Description[match(edges$to,gosim_res_up$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out1,gsub('/','_',dataSlice[ss]),"_GO_network_edges_up.txt"), sep = "\t", row.names = F,quote = F)
    
    
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
      filter(weight > 0.1) 
    edges$from_Description<-gosim_res_down$Description[match(edges$from,gosim_res_down$ID)]
    edges$to_Description<-gosim_res_down$Description[match(edges$to,gosim_res_down$ID)]
    write.table(edges[,-c(1,2)], paste0(dir_out1,gsub('/','_',dataSlice[ss]),"_GO_network_edges_down.txt"), sep = "\t", row.names = F,quote = F)

  }
  
  print(dataSlice[ss])
}


