###8中T细胞状态
library(TCellSI)
library(dplyr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(monocle)
library(igraph)
library(RobustRankAggreg)
library(jsonlite)


file_RCTDstep<-list.files(pattern = 'public_step10',path = dir_step,recursive = T)
patient <- c("nr","r")

dir_Tcell <- "/TcellSI/"
file_Tcell<-list.files(pattern = 'TCellSI.txt',path = dir_Tcell,recursive = T)

dataSlice<-unlist(lapply(strsplit(file_Tcell,'_T'),function(x)x[1]))
#file_near<-file_near[match(paste0(dataSlice,'_BdyNearCore.txt'),file_near)]
file_Tcell_name <- gsub("/", "_", file_Tcell)
file_Tcell_name <- gsub("_TCellSI.txt", "", file_Tcell_name)


g_nr_positions <- which(file_Tcell_name %in% g_nr)
g_r_positions <- which(file_Tcell_name %in% g_r_noTIB)
g_nr_tib_positions <- which(file_Tcell_name %in% g_nr_TIB)

st_Tcell_nr_tib <- c()
for(i in g_nr_tib_positions){
  #i <- 1
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  rownames(st_Tcell) <- paste0(file_Tcell_name[i],"_",rownames(st_Tcell))
  st_Tcell_nr_tib <- rbind(st_Tcell_nr_tib,st_Tcell)
}
write.table(st_Tcell_nr_tib,file = "/public/NR_TIB_TCellSI.txt",quote = F,sep = '\t',row.names = T)

st_Tcell_nr <- c()
for(i in g_nr_positions){
  #i <- 1
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  rownames(st_Tcell) <- paste0(file_Tcell_name[i],"_",rownames(st_Tcell))
  st_Tcell_nr <- rbind(st_Tcell_nr,st_Tcell)
}
write.table(st_Tcell_nr,file = "/public/NR_TCellSI.txt",quote = F,sep = '\t',row.names = T)
st_Tcell_r <- c()
for(i in g_r_positions){
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  rownames(st_Tcell) <- paste0(file_Tcell_name[i],"_",rownames(st_Tcell))
  st_Tcell_r <- rbind(st_Tcell_r,st_Tcell)
  
}
write.table(st_Tcell_r,file = "/public/R_TCellSI.txt",quote = F,sep = '\t',row.names = T)

dir_Tcell <- "/public/"
file_Tcell<-list.files(pattern = 'TCellSI.txt',path = dir_Tcell,recursive = T)

all_step_Tcell<-c()
for(i in 1:2){
  #i=2
  st_near<-read.delim(paste0(dir_step,file_RCTDstep[i]),stringsAsFactors = F,check.names = F)
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  length(intersect(spot_step1,rownames(st_Tcell)))
  spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
  spot_step2<-union(spot_step1,spot_step2)
  spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
  spot_step3<-union(spot_step2,spot_step3)
  spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
  spot_step4<-union(spot_step3,spot_step4)
  spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
  spot_step5<-union(spot_step4,spot_step5)
  
  
  step_Tcell<-lapply(list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5),function(x){
    spot_Tcell<-intersect(x,rownames(st_Tcell))
    near_Tcell<-rep(0,8)
    names(near_Tcell)<-colnames(st_Tcell)
    if(length(spot_Tcell)>1){
      near_Tcell<-apply(st_Tcell[spot_Tcell,],2,mean)
    }
    return(near_Tcell)
  })
  step_Tcell<-do.call(rbind,step_Tcell)
  rownames(step_Tcell)<-paste0('step_',1:5)
  step_Tcell<-reshape2::melt(as.matrix(step_Tcell))
  step_Tcell$slice<-patient[i]
  
  all_step_Tcell<-rbind(all_step_Tcell,step_Tcell)
}

#all_step_Tcell$patient<-unlist(lapply(strsplit(all_step_Tcell$slice,'/'),function(x)x[1]))
#all_step_Tcell$cancer<-substr(all_step_Tcell$cancer,1,nchar(all_step_Tcell$cancer)-2) %>% toupper()
all_step_Tcell<-all_step_Tcell[which(all_step_Tcell$value!=0),]
#cancer<-unique(all_step_Tcell$cancer)



####绘制成热图 每个T细胞状态画一个图
###一行一个切片 一列是step
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

all_step_Tcell$Var2<-as.vector(all_step_Tcell$Var2)
all_step_Tcell$Var1<-as.vector(all_step_Tcell$Var1)
type <- unique(all_step_Tcell$Var2)


#### CAF 比较
pdf('\\TcellSI_compare_public1.pdf',width = 4,height = 4)
for( i in 1:length(unique(all_step_Tcell$Var2))){
  
  p_data<-all_step_Tcell[all_step_Tcell$Var2 %in% type[i],]
  TcellSI_name<-unique(all_step_Tcell$Var2)
  
  p_data<-p_data[p_data$Var1%in%paste0('step_',1:6),]
  p_box<-aggregate(p_data$value,by=list(p_data$slice,p_data$Var1),mean)
  colnames(p_box)<-c('slice','step','value')
  p_box<-p_box[order(p_box$slice),]
  
  library(ggpubr)
  p <- ggpaired(p_box, 
                x = "slice", 
                y = "value",
                color = "slice", 
                line.color = "gray", 
                line.size = 0.4,
                palette = "jco") +
    stat_compare_means(paired = TRUE)+
    labs(title = type[i])
  
  #print(p)
  p_box <- ggplot(p_box, aes(x = slice, y = value, fill = slice)) + 
    stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
    geom_boxplot(aes(fill = slice), color = 'black', width = 0.8,
                 position = position_dodge(0.9), alpha = 1, outlier.shape = NA) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = 'black')) +
    scale_fill_manual(values = c("nr_tib" = "#15629e", "r" = "#ddb424")) +
    scale_color_manual(values = c("nr_tib" = "#15629e", "r" = "#ddb424")) +
    ggtitle(paste0('ImmLigand_', i)) +
    # 添加p值（使用非配对检验，因为数据不平衡）
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                       method = "wilcox.test",  # 使用Wilcoxon秩和检验（非配对）
                       label.x = 1.5,
                       label.y = 1.8) +  # 设置p值位置在y=1.8
    # 设置y轴范围为0-2
    coord_cartesian(ylim = c(0, 2)) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10))
  print(p_box)
}
dev.off()

pdf('\\ESCC\\Deconvolution\\TcellSI_compare_public.pdf',width = 4,height = 4)
print(p)
dev.off()

##nr_tib vs r
all_step_Tcell<-c()
patient <- c("nr","nr_tib","r")
for(i in 2:3){
  #i=2
  st_near<-read.delim(paste0(dir_step,file_RCTDstep[i]),stringsAsFactors = F,check.names = F)
  st_Tcell<-read.delim(paste0(dir_Tcell,file_Tcell[i]),stringsAsFactors = F,check.names = F)
  
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  length(intersect(spot_step1,rownames(st_Tcell)))
  spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
  spot_step2<-union(spot_step1,spot_step2)
  spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
  spot_step3<-union(spot_step2,spot_step3)
  spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
  spot_step4<-union(spot_step3,spot_step4)
  spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
  spot_step5<-union(spot_step4,spot_step5)
  
  
  step_Tcell<-lapply(list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5),function(x){
    spot_Tcell<-intersect(x,rownames(st_Tcell))
    near_Tcell<-rep(0,8)
    names(near_Tcell)<-colnames(st_Tcell)
    if(length(spot_Tcell)>1){
      near_Tcell<-apply(st_Tcell[spot_Tcell,],2,mean)
    }
    return(near_Tcell)
  })
  step_Tcell<-do.call(rbind,step_Tcell)
  rownames(step_Tcell)<-paste0('step_',1:5)
  step_Tcell<-reshape2::melt(as.matrix(step_Tcell))
  step_Tcell$slice<-patient[i]
  
  all_step_Tcell<-rbind(all_step_Tcell,step_Tcell)
}

#all_step_Tcell$patient<-unlist(lapply(strsplit(all_step_Tcell$slice,'/'),function(x)x[1]))
#all_step_Tcell$cancer<-substr(all_step_Tcell$cancer,1,nchar(all_step_Tcell$cancer)-2) %>% toupper()
all_step_Tcell<-all_step_Tcell[which(all_step_Tcell$value!=0),]
#cancer<-unique(all_step_Tcell$cancer)



####绘制成热图 每个T细胞状态画一个图
###一行一个切片 一列是step
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

all_step_Tcell$Var2<-as.vector(all_step_Tcell$Var2)
all_step_Tcell$Var1<-as.vector(all_step_Tcell$Var1)
type <- unique(all_step_Tcell$Var2)


#### CAF 比较
pdf('\\TcellSI_compare_public_nrtib_vs_r2.pdf',width = 4,height = 4)
for( i in 1:length(unique(all_step_Tcell$Var2))){
  
  p_data<-all_step_Tcell[all_step_Tcell$Var2 %in% type[i],]
  TcellSI_name<-unique(all_step_Tcell$Var2)
  
  p_data<-p_data[p_data$Var1%in%paste0('step_',1:6),]
  p_box<-aggregate(p_data$value,by=list(p_data$slice,p_data$Var1),mean)
  colnames(p_box)<-c('slice','step','value')
  p_box<-p_box[order(p_box$slice),]
  
  library(ggpubr)
  p <- ggpaired(p_box, 
                x = "slice", 
                y = "value",
                color = "slice", 
                line.color = "gray", 
                line.size = 0.4,
                palette = "jco") +
    stat_compare_means(paired = TRUE)+
    labs(title = type[i])
  
  #print(p)
  p_box <- ggplot(p_box, aes(x = slice, y = value, fill = slice)) + 
    stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
    geom_boxplot(aes(fill = slice), color = 'black', width = 0.8,
                 position = position_dodge(0.9), alpha = 1, outlier.shape = NA) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = 'black')) +
    scale_fill_manual(values = c("nr_tib" = "#15629e", "r" = "#ddb424")) +
    scale_color_manual(values = c("nr_tib" = "#15629e", "r" = "#ddb424")) +
    ggtitle(type[i]) +
    # 添加p值（使用非配对检验，因为数据不平衡）
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                       method = "wilcox.test",  # 使用Wilcoxon秩和检验（非配对）
                       label.x = 1.5,
                       label.y = 0.0035) +  # 设置p值位置在y=1.8
    # 设置y轴范围为0-2
    #coord_cartesian(ylim = c(0, 2)) +
    theme(axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 10),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 10))
  print(p_box)
  
}
dev.off()

pdf('\\ESCC\\Deconvolution\\TcellSI_compare_public.pdf',width = 4,height = 4)
print(p)
dev.off()


