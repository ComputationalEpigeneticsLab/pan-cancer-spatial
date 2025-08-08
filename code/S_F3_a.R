#######inferCNV
library(coda)
library(Seurat)
library(tidyverse)

#Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.1") 
library(rjags)
# BiocManager::install("infercnv",force=TRUE)
library(infercnv)
library(Matrix)





#######################相关文件构造
dir_in<-"E:/Mirror/ST_analysis/data/10X Visium/1ROGUE/"###/home/xujuan.hyd/ST_analysis/data/10X Visium/1ROGUE
file_anno<-list.files(pattern = "_cellAnnotation.txt",recursive = T,path = dir_in)
file_count<-list.files(pattern = "_ST_count.txt",recursive = T,path = dir_in)
table(unlist(lapply(strsplit(file_anno,"_"),function(x){x[1]}))==
        unlist(lapply(strsplit(file_count,"_"),function(x){x[1]})))
dataSet_name<-unlist(lapply(strsplit(file_anno,"/"),function(x){x[1]}))
slice_name<-unlist(lapply(strsplit(file_anno,"_"),function(x){x[1]}))
slice_name<-unlist(lapply(strsplit(slice_name,"/"),function(x){x[2]}))
file_count[165]

for(i in 1:length(file_anno)){
  #i=165
  #dir_out<-paste0("E:/Mirror/ST_analysis/data/10X Visium/inferCNV_R/",dataSet_name[i],"/")
  dir_out<-paste0("E:/Mirror/ST_analysis/data/10X Visium/inferCNV_R/",dataSet_name[i],"/")
  if(!dir.exists(dir_out)){dir.create(dir_out)}
  
  ###count表达矩阵
  ST_count<-read.delim(paste0(dir_in,file_count[i]),stringsAsFactors = F,check.names = F)
  
  ###注释文件
  cellAnnotations<-read.delim(paste0(dir_in,file_anno[i]),stringsAsFactors = F,check.names = F)
  cellAnnotations<-cellAnnotations[,-2]
  table(cellAnnotations$cellType)
  write.table(cellAnnotations,paste0(dir_in,dataSet_name[i],"/",slice_name[i],"_cellAnno_Rcnv.txt"),
              col.names = F,row.names = F,sep = "\t",quote = F)
  
  ###参考基因组
  #gene_anno <- read.delim("E:/Mirror/ST_analysis/data/gencode.v43.annotation.txt",
  #                        sep="\t",stringsAsFactors=F)
  gene_anno <- read.delim("E:/Mirror/ST_analysis/data//gencode.v43.annotation.txt",
                          sep="\t",stringsAsFactors=F)
  gene_anno1 <- gene_anno[gene_anno[,5] %in% rownames(ST_count),]
  gene_anno1 <- gene_anno1[!duplicated(gene_anno1[,5]),]
  write.table(gene_anno1[,c(5,1,2,3)],
              file=paste0(dir_in,dataSet_name[i],"/",slice_name[i],"_gene_ordering_file.txt"),
              sep="\t",row.names=F,col.names=F,quote=F)
  
  
  #第一步，根据上述的三个文件创建inferCNV对象
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(ST_count), # 可以直接提供矩阵对象
                                      annotations_file=paste0(dir_in,dataSet_name[i],"/",slice_name[i],"_cellAnno_Rcnv.txt"),###细胞类型文件
                                      delim="\t",
                                      gene_order_file=paste0(dir_in,dataSet_name[i],"/",slice_name[i],"_gene_ordering_file.txt"),##处理的参考基因组文件
                                      ref_group_names=NULL)
  
  #第二步，运行标准的inferCNV流程。
  # perform infercnv operations to reveal cnv signal
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir=paste0(dir_out,slice_name[i]),  # 输出文件夹
                               cluster_by_groups=F,   # 聚类
                               plot_steps=F,
                               no_plot = T,
                               output_format = "pdf",
                               scale_data=T,
                               noise_filter=0.12,
                               analysis_mode="subclusters",
                               #min_cells_per_gene = 0,
                               denoise=T, #去噪
                               HMM_type='i6',
                               tumor_subcluster_partition_method ='random_trees',
                               HMM=T,# 是否基于HMM预测CNV
                               #write_expr_matrix = T,
                               num_threads = 1
  )
}



####inferCNV的扩增和缺失在core和非core恶性spot中的差别
library(coda)
library(Seurat)
library(infercnv)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(magic)
library(rlang)
library(ggplot2)
library(ggraph)
library(ggpubr)
library(psych)


dir_cnvRe<-'E:/Mirror/ST_analysis/data/10X Visium/CNV_stem/'
###扩增或缺失统计结果
chr_allMean<-read.delim(paste0(dir_cnvRe,'chr_coreGreater_mean.txt'),stringsAsFactors = F,check.names = F)
chr_allMean$slices<-paste0(chr_allMean$cancer,'_',chr_allMean$slices)
chr_allMean_sig<-aggregate(chr_allMean[,3],by=list(chr_allMean$slices),function(x) length(which(x<0.05)))
chr_allMean$diff<-chr_allMean$meanCore-chr_allMean$meanBdyBud
chr_allMean_diff<-aggregate(chr_allMean[,7],by=list(chr_allMean$slices),function(x) sum(x[which(x>0)]))

chr_allMean_sig<-chr_allMean_sig[order(chr_allMean_sig$x,decreasing = T),]
chr_allMean_diff<-chr_allMean_diff[order(chr_allMean_diff$x,decreasing = T),]

select_slice1<-intersect(chr_allMean_sig$Group.1[1:50],chr_allMean_diff$Group.1[1:50])
select_slice1<-union(select_slice1,'cscc02_slice4')

select_slice1<-lapply(unique(chr_allMean$chr),function(x){#x='chr1'
  aa_chr<-chr_allMean[chr_allMean$chr%in%x,]
  aa_chr<-aa_chr[which(aa_chr$p<0.01),]
  aa_chr<-aa_chr[order(aa_chr$diff,decreasing = T),]
  return(aa_chr$slices[1:2])###aa_chr$slices[1:2]
})
names(select_slice1)<-unique(chr_allMean$chr)
select_slice1<-as.matrix(t(do.call(cbind,select_slice1)))
select_slice1<-reshape2::melt(select_slice1)

###参考基因组 
geneFile <- read.table("E:/Mirror/ST_analysis/data/gencode.v43.annotation.txt",header = T,sep = "\t",stringsAsFactors = F)
geneFile<-geneFile[!duplicated(geneFile$gene_name),]
rownames(geneFile)=geneFile$gene_name


dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
#file_bdy<-file_bdy[-grep('ovca02',file_bdy)]  ####还在跑
file_bdy<-file_bdy[-grep('luad01/slice2',file_bdy)] ###########################这个切片跑inferCNV时未知原因报错
data_slice<-unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x) x[1]))
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x) x[1]))

dir_CNV<-'E:/Mirror/ST_analysis/data/10X Visium/inferCNV_R2/'
file_clone<-list.files(pattern = 'hmm_mode-subclusters.cell_groupings',path = dir_CNV,recursive = T)
file_clone<-file_clone[-grep('luad01_slice2',file_clone)]
file_CNV<-list.files(pattern = 'run.final.infercnv_obj',path = dir_CNV,recursive = T)
table(unlist(lapply(strsplit(file_clone,'/'),function(x) x[1]))==unlist(lapply(strsplit(file_CNV,'/'),function(x) x[1])))
table(unlist(lapply(strsplit(file_clone,'/'),function(x) x[1]))==gsub('/','_',data_slice))

file_bdy<-paste0(gsub('_','/',select_slice1$value),'_BdyTumorCore.txt')
file_CNV<-paste0(select_slice1$value,'/run.final.infercnv_obj')
data_slice<-unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x) x[1]))
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x) x[1]))



##################################################################################################
###绝对变异
dir_cnvRe<-'E:/Mirror/ST_analysis/data/10X Visium/CNV_stem/'
###扩增或缺失统计结果
chr_allMean<-read.delim(paste0(dir_cnvRe,'chr_coreGreater_abs.txt'),stringsAsFactors = F,check.names = F)
chr_allMean$slices<-paste0(chr_allMean$cancer,'_',chr_allMean$slices)
chr_allMean_sig<-aggregate(chr_allMean[,3],by=list(chr_allMean$slices),function(x) length(which(x<0.05)))
chr_allMean$diff<-chr_allMean$meanCore-chr_allMean$meanBdyBud
chr_allMean_diff<-aggregate(chr_allMean[,7],by=list(chr_allMean$slices),function(x) sum(abs(x)))

chr_allMean_sig<-chr_allMean_sig[order(chr_allMean_sig$x,decreasing = T),]
chr_allMean_diff<-chr_allMean_diff[order(chr_allMean_diff$x,decreasing = T),]

select_slice2<-intersect(chr_allMean_sig$Group.1[1:50],chr_allMean_diff$Group.1[1:50])


####每个切片的扩增缺失条形图
library(coda)
library(Seurat)
library(infercnv)
library(Matrix)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(magic)


###参考基因组 
geneFile <- read.table("E:/Mirror/ST_analysis/data/gencode.v43.annotation.txt",header = T,sep = "\t",stringsAsFactors = F)
geneFile<-geneFile[!duplicated(geneFile$gene_name),]
rownames(geneFile)=geneFile$gene_name


dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
#file_bdy<-file_bdy[-grep('ovca02',file_bdy)]  ####还在跑
file_bdy<-file_bdy[-grep('luad01/slice2',file_bdy)] ###########################这个切片跑inferCNV时未知原因报错
data_slice<-unlist(lapply(strsplit(file_bdy,'_BdyTumorCore'),function(x) x[1]))
dataSet<-unlist(lapply(strsplit(data_slice,'/'),function(x) x[1]))

dir_CNV<-'E:/Mirror/ST_analysis/data/10X Visium/inferCNV_R2/'
file_clone<-list.files(pattern = 'hmm_mode-subclusters.cell_groupings',path = dir_CNV,recursive = T)
file_clone<-file_clone[-grep('luad01_slice2',file_clone)]
file_CNV<-list.files(pattern = 'run.final.infercnv_obj',path = dir_CNV,recursive = T)
table(unlist(lapply(strsplit(file_clone,'/'),function(x) x[1]))==unlist(lapply(strsplit(file_CNV,'/'),function(x) x[1])))
table(unlist(lapply(strsplit(file_clone,'/'),function(x) x[1]))==gsub('/','_',data_slice))


dir_pic<-'E:/Mirror/ST_analysis/pic/clone/CNV_bar/'
pdf(paste0(dir_pic,'CNV_core_bar.pdf'),width = 12,height = 5)
for(i in 1:length(file_bdy)){
  #i=108
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  table(st_bdy$FinalLocalType)
  
  
  if(length(which(st_bdy$FinalLocalType=="Core"))>10){
    infercnv_obj_example = readRDS(paste0(dir_CNV,file_CNV[i]))
    expr <- infercnv_obj_example@expr.data
    dim(expr)
    gn <- rownames(expr)
    
    sub_geneFile <-  geneFile[intersect(gn,geneFile$gene_name),c(5,1,2,3)]
    expr=expr[intersect(gn,geneFile$gene_name),rownames(st_bdy)[which(st_bdy$FinalLocalType=='Core')]]
    expr<-expr-1
    
    CNV_mean<-apply(expr,1,function(x){
      amp_x<-mean(x[which(x>0)])
      del_x<-mean(x[which(x<0)])
      return(c(amp_x,del_x))
    }) %>% t() 
    colnames(CNV_mean)<-c('amp','del')
    CNV_mean[is.nan(CNV_mean)]<-0
    CNV_mean<-as.data.frame(CNV_mean)
    
    plot_data<-data.frame(type=rep(c('amp','del'),each=nrow(CNV_mean)),
                          chr=c(sub_geneFile$chrom,sub_geneFile$chrom),
                          gene=c(rownames(CNV_mean),rownames(CNV_mean)),
                          value=c(CNV_mean$amp,CNV_mean$del))
    plot_data$gene<-factor(plot_data$gene,levels = unique(plot_data$gene))
    chr_bk<-as.data.frame(table(sub_geneFile$chrom))
    class(chr_bk$Var1)
    chr_bk$order<-as.numeric(unlist(lapply(strsplit(as.vector(chr_bk$Var1),'chr'),function(x)x[2])))
    chr_bk<-chr_bk[order(chr_bk$order),]
    chr_bk$site<-cumsum(chr_bk$Freq)+0.5
    chr_bk$lab_site<-chr_bk$site-chr_bk$Freq/2
    
    #plot_data<-plot_data[plot_data$chr%in%c('chr1','chr2'),]
    p <- ggplot(plot_data, aes(x=gene, y=value, fill=type)) +
      geom_bar(stat = "identity", width = 0.8) +
      labs(y="CNV", x="chr",title = paste0('Core_',data_slice[i])) +
      #geom_vline(xintercept=chr_bk$Freq,lty=4,col="grey70",lwd=0.6) +
      geom_vline(xintercept=chr_bk$site,lty=2,col="grey30",lwd=0.6) +
      # geom_text(aes(label = as.vector(chr_bk$Var1),x=chr_bk$lab_site,y=rep(c(0.3,0.6),11)),
      #           #position = position_fill(vjust = 0.5),
      #           #hjust = 3,vjust = 0.5,
      #           size=3) +
      annotate("text", x=chr_bk$lab_site,y=rep(c(0.05,0.1),11)[1:nrow(chr_bk)],label = as.vector(chr_bk$Var1),colour="black",size=3) +
      theme_classic(base_size = 15) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.y = element_line(),
            legend.position = "top"
      ) +
      scale_fill_manual(values = c('amp'='#ed0000','del'='#00448b'))
    print(p)
  }
  print(data_slice[i])
  
}
dev.off()




dir_pic<-'E:/Mirror/ST_analysis/pic/clone/CNV_bar/'
pdf(paste0(dir_pic,'CNV_BudBdy_bar.pdf'),width = 12,height = 5)
for(i in 1:length(file_bdy)){
  #i=1
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  table(st_bdy$FinalLocalType)
  
  
  if(length(which(st_bdy$FinalLocalType=="Boundary"|st_bdy$FinalLocalType=="Dispersion"))>10){
    infercnv_obj_example = readRDS(paste0(dir_CNV,file_CNV[i]))
    expr <- infercnv_obj_example@expr.data
    dim(expr)
    gn <- rownames(expr)
    
    sub_geneFile <-  geneFile[intersect(gn,geneFile$gene_name),c(5,1,2,3)]
    expr=expr[intersect(gn,geneFile$gene_name),
              rownames(st_bdy)[which(st_bdy$FinalLocalType=="Boundary"|st_bdy$FinalLocalType=="Dispersion")]]
    expr<-expr-1
    
    CNV_mean<-apply(expr,1,function(x){
      amp_x<-mean(x[which(x>0)])
      del_x<-mean(x[which(x<0)])
      return(c(amp_x,del_x))
    }) %>% t() 
    colnames(CNV_mean)<-c('amp','del')
    CNV_mean[is.nan(CNV_mean)]<-0
    CNV_mean<-as.data.frame(CNV_mean)
    
    plot_data<-data.frame(type=rep(c('amp','del'),each=nrow(CNV_mean)),
                          chr=c(sub_geneFile$chrom,sub_geneFile$chrom),
                          gene=c(rownames(CNV_mean),rownames(CNV_mean)),
                          value=c(CNV_mean$amp,CNV_mean$del))
    plot_data$gene<-factor(plot_data$gene,levels = unique(plot_data$gene))
    chr_bk<-as.data.frame(table(sub_geneFile$chrom))
    class(chr_bk$Var1)
    chr_bk$order<-as.numeric(unlist(lapply(strsplit(as.vector(chr_bk$Var1),'chr'),function(x)x[2])))
    chr_bk<-chr_bk[order(chr_bk$order),]
    chr_bk$site<-cumsum(chr_bk$Freq)+0.5
    chr_bk$lab_site<-chr_bk$site-chr_bk$Freq/2
    
    #plot_data<-plot_data[plot_data$chr%in%c('chr1','chr2'),]
    p <- ggplot(plot_data, aes(x=gene, y=value, fill=type)) +
      geom_bar(stat = "identity", width = 0.8) +
      labs(y="CNV", x="chr",title = paste0('BudBdy_',data_slice[i])) +
      #geom_vline(xintercept=chr_bk$Freq,lty=4,col="grey70",lwd=0.6) +
      geom_vline(xintercept=chr_bk$site,lty=2,col="grey30",lwd=0.6) +
      annotate("text", x=chr_bk$lab_site,y=rep(c(0.05,0.1),11)[1:nrow(chr_bk)],label = as.vector(chr_bk$Var1),colour="black",size=3) +
      theme_classic(base_size = 15) +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major.y = element_line(),
            legend.position = "top"
      ) +
      scale_fill_manual(values = c('amp'='#ed0000','del'='#00448b'))
    print(p)
  }
  print(data_slice[i])
}
dev.off()



