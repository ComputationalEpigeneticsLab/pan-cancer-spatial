####每种癌症类型取core的计算各个MP之间的相关性
library(corrplot)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(pheatmap)


integration_umap_MP<-read.delim('E:/Mirror/ST_analysis/data/pic_data/integration_umap_MP.txt',stringsAsFactors = F,check.names = F)
integration_umap_MP_core<-integration_umap_MP[which(integration_umap_MP$LocationType=='Core'),]
table(integration_umap_MP$LocationType)

cancer<-unlist(lapply(strsplit(integration_umap_MP_core$slice,'/'),function(x) x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()
cancer<-c(cancer,'slice')

color_pheatmap <- c(colorRampPalette(c("#0335FF","#809AFF"))(10),
                    colorRampPalette(c("#809AFF","white"))(10),
                    colorRampPalette(c("white","#DB878F"))(10),
                    colorRampPalette(c("#DB878F","#BE2331"))(10))
dir_pic<-'E:/Mirror/ST_analysis/pic/NMF_module/MP_corr/'

for(i in 1:length(cancer)){
  #i=20
  cancer_MP<-integration_umap_MP_core[grep(cancer[i],integration_umap_MP_core$slice),paste0('MP_',2:14)]
  
  cancer_MP<-cancer_MP[,]
  #tdc<-cor(cancer_MP, method="spearman")
  #testRes = corrplot::cor.mtest(cancer_MP, method="spearman",conf.level = 0.95)
  corr_re<-psych::corr.test(cancer_MP,cancer_MP, method="spearman", adjust = "fdr")
  dis_dot<-as.matrix(corr_re[["p.adj"]])
  dis_dot[as.matrix(corr_re[["p.adj"]]) >= 0.05 | abs(as.matrix(corr_re[["r"]]))<0.5]<-1
  diag(dis_dot)<-1
  
  me<-c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid")
  for(j in 1:length(me)){#j=7
    pdf(paste0(dir_pic,cancer[i],'_',me[j],'_corr_new2.pdf'),width = 8,height = 8)
    corrplot(corr_re$r, method = "circle", col = color_pheatmap, 
             order = "hclust",hclust.method = me[j],
             tl.col = "black", tl.cex = 0.8, tl.srt = 45,#tl.pos = "lt",
             #p.mat = dis_dot, ###一个显著性都没有时不能用
             diag = T, type = 'upper',###显示一半及对角线
             addCoef.col = 'grey10',number.cex = 0.8,###相关性数字颜色及大小
             title = cancer[i],
             sig.level = c(0.05), pch.cex = 1.5,
             insig = 'label_sig', pch.col = 'grey100')
    dev.off()
  }
  
}



###Misty
####在服务器中运行
library(dplyr)
library(mistyR)
library(future)
library(purrr)
library(distances)
#plotting
library(ggplot2)
#plan(multisession)


#dir_RCTD<-'/data/zhouweiwei/spatial/data/10X Visium/RCTD3/'
dir_MP<-'/data/zhouweiwei/spatial/data/10X Visium/moduleScore/'
dir_bdy<-'/data/zhouweiwei/spatial/data/10X Visium/copykat/'
file_MP<-list.files(pattern = 'MP_score.txt',path = dir_MP,recursive = T)
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_MP,'_'),function(x)x[1]))
table(dataSlice==unlist(lapply(strsplit(file_bdy,'_'),function(x)x[1])))

for(i in 1:229){
  #i=1
  dir_out<-paste0('/data/zhouweiwei/spatial/data/10X Visium/MistyR/MP2/',gsub('/','_',dataSlice[i]))
  if(!dir.exists(dir_out)) dir.create(dir_out)
  
  setwd(dir_out)
  st_MP<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType=='Core'),]
  tryCatch({
    st_MP<-st_MP[rownames(st_bdy),-1]
    #colnames(st_MP)<-gsub(' ','_',colnames(st_MP))
    expr<-st_MP############
    misty.intra <- create_initial_view(expr)
    pos<-data.frame(row=st_bdy$imagerow,col=st_bdy$imagecol)##############
    misty.views <- misty.intra %>% add_paraview(pos, l = 10)
    aa<-run_misty(misty.views,bypass.intra = T)
    misty.results <- collect_results(aa)
    saveRDS(misty.results,file = 'misty.results.rds')
  },error = function(e){
    stopMessage<-"something wrong"
    return(stopMessage)
  })
  print(dataSlice[i])
}



library(dplyr)
library(mistyR)
library(future)
library(purrr)
library(distances)
#plotting
library(ggplot2)
library(Seurat)
library(ComplexHeatmap)
library(circlize)


dir_Misty<-'E:/Mirror/ST_analysis/data/10X Visium/MistyR/MP/'
file_Misty<-list.files(pattern = 'misty.results.rds',path = dir_Misty,recursive = T)
data_slice<-unlist(lapply(strsplit(file_Misty,'/misty'),function(x)x[1]))
cancer<-unlist(lapply(strsplit(data_slice,'_'),function(x)x[1]))
cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()

dir_pic<-'E:/Mirror/ST_analysis/pic/Misty/'

####按所有切片求中值
misty.results<-readRDS(paste0(dir_Misty,file_Misty[1]))
plot.data2 <- misty.results$importances.aggregated
plot.data2<-plot.data2[which(plot.data2$view=='intra'),]
heat_data<-plot.data2[,c('Predictor','Target','Importance')]
all_heat_data<-heat_data

for(j in 2:length(file_Misty)){#j=2
  misty.results<-readRDS(paste0(dir_Misty,file_Misty[j]))
  plot.data2 <- misty.results$importances.aggregated
  plot.data2<-plot.data2[which(plot.data2$view=='intra'),]
  all_heat_data<-cbind(all_heat_data,plot.data2$Importance)
  
}

heat_data<-data.frame(Predictor=all_heat_data$Predictor,
                      Target=all_heat_data$Target,
                      Importance=apply(all_heat_data[,3:ncol(all_heat_data)],1,median))
heat_data<-reshape2::acast(heat_data,Predictor~Target)
#heat_data[which(heat_data<0)]<-0
heat_data<-heat_data[c('MP_11','MP_8','MP_10','MP_2','MP_14','MP_3','MP_12','MP_9','MP_6','MP_7','MP_4','MP_5','MP_13'),
                     c('MP_11','MP_8','MP_10','MP_2','MP_14','MP_3','MP_12','MP_9','MP_6','MP_7','MP_4','MP_5','MP_13')]
write.table(heat_data,'E:/Mirror/ST_analysis/data/pic_data/Misty/heat_data.txt',quote = F,sep = '\t')

heat_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/Misty/heat_data.txt',stringsAsFactors = F,check.names = F)

col_fun <- colorRamp2(
  c(min(heat_data,na.rm = T), max(heat_data,na.rm = T)),
  c("#FFFFFF",'#004D10')
)
p1<-Heatmap(as.matrix(heat_data),
            col = col_fun,
            #top_annotation =ha,#####顶部注释
            #left_annotation=la,#####左侧注释
            show_row_names = T,#####不显示行名
            show_column_names =T,#####不显示列名
            cluster_rows =F,
            na_col = 'white',
            cluster_columns = F,
            #name="TCGA_Immune_infiltration",
            column_title=paste0("Misty_",'allSlice'),
            name = 'Importance',
            rect_gp = gpar(col = 'grey85', lwd = 1.5)
)
print(p1)
pdf(paste0(dir_pic,'allslice_Misty_para.pdf'),width = 5,height = 4.2)
print(p1)
dev.off()



















