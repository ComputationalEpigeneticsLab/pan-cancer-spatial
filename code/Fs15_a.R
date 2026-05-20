#######不同step的细胞类型的变化
library(tidyverse)
library(Seurat)


hop_infor_imm<-read.table("\\new\\hop_infor_imm.txt",sep = "\t",header = T)

g_r_noTIB=c("GBM_GSE235672_GSM7507327","GBM_GSE235672_GSM7507330",
            "LIHC_GSE238264_GSM7661255","LIHC_GSE238264_GSM7661256","LIHC_GSE238264_GSM7661257","LIHC_GSE238264_GSM7661258",'LIHC_lihc03_slice7' )

g_nr_TIB=c("GBM_GSE235672_GSM7507312","GBM_GSE235672_GSM7507328","LIHC_GSE238264_GSM7661260",'LIHC_lihc03_slice5')

g_nr_noTIB=c("GBM_GSE235672_GSM7507323","GBM_GSE235672_GSM7507329","LIHC_GSE238264_GSM7661259","LIHC_GSE238264_GSM7661261")

g_nr <-c(g_nr_TIB,g_nr_noTIB)

hop_infor_imm_r <- hop_infor_imm[hop_infor_imm$slice %in% g_r_noTIB,]
hop_infor_imm_nr <- hop_infor_imm[hop_infor_imm$slice %in% g_nr,]

rctd_data<-read.table("\\new\\deconvolution_rctd_infor_imm.txt",sep = "\t",header = T)
rctd_data$cancer_type <- sapply(strsplit(rctd_data$X, "_"), function(x) x[1])
rctd_data$gse_id <- sapply(strsplit(rctd_data$X, "_"), function(x) x[2])
rctd_data$gsm_id <- sapply(strsplit(rctd_data$X, "_"), function(x) x[3])

rctd_data$slice <- paste0(rctd_data$cancer_type,"_",rctd_data$gse_id,"_",rctd_data$gsm_id)
rownames(rctd_data) <- rctd_data$X


rctd_data_r <- rctd_data[rctd_data$slice %in% g_r_noTIB,]
rctd_data_nr <- rctd_data[rctd_data$slice %in% g_nr,]
rctd_data_nr_tib <- rctd_data[rctd_data$slice %in% g_nr_TIB,]
colnames(hop_infor_imm)
#st_Deco<-st_Deco[,setdiff(colnames(st_Deco),c('Endothelial','Epithelial'))]


#######R组
st_step<-hop_infor_imm_r[,c("X","imagerow","imagecol","FinalLocalType",paste0('step_',1:10))]
#st_step<-st_step[rownames(rctd_data_r),]
rownames(st_step) <- st_step$X
st_step <- st_step[,-1]
#st_step<-st_step[which(st_step$FinalLocalType=='Normal'|st_step$FinalLocalType=='Immune'),]
st_Deco<-rctd_data_r[rownames(st_step)[st_step$FinalLocalType%in%c('Normal','Immune')],]
st_Deco <- st_Deco[,-1]
st_Deco <- st_Deco[,-c(13:16)]


# 提取每行的前缀（从行名中获取）
row_prefixes <- gsub("_[^_]+$", "", rownames(st_step))
#cat("行名前缀:\n")
#print(row_prefixes)

# 方法1: 为每个步骤补全细胞名
st_step_full <- st_step

# 处理step_1到step_10列
step_cols <- grep("^step_", colnames(st_step), value = TRUE)

for(ss in step_cols) {
  cat("\n处理", ss, "\n")
  
  # 为每一行补全细胞名
  full_names <- sapply(1:nrow(st_step), function(i) {
    current_cells <- as.character(st_step[i, ss])
    if(is.na(current_cells) || current_cells == "") return("")
    
    # 分割细胞名
    cells <- unlist(strsplit(current_cells, ','))
    cells <- cells[!is.na(cells) & cells != ""]
    
    if(length(cells) == 0) return("")
    
    # 补全为完整格式
    full_cells <- paste0(row_prefixes[i], "_", cells)
    return(paste(full_cells, collapse = ','))
  })
  
  st_step_full[, ss] <- full_names
}
write.table(st_step_full,file = paste0(dir_out,'public_step10_r.txt'),quote = F,sep = '\t',row.names = T)
st_step_CAF<-st_step_full[st_step_full$FinalLocalType%in%c('Boundary','Dispersion'),]

for(jj in ncol(st_step_CAF):6){#jj=6
  y_spot<-unlist(strsplit(st_step_CAF[,(jj-1)],',')) %>% unique()
  st_step_CAF[,jj]<-lapply(1:nrow(st_step_CAF),function(x){##x=1
    x_spot<-strsplit(st_step_CAF[x,jj],',') %>% unlist()
    x_spot<-setdiff(x_spot,y_spot)
    return(paste0(x_spot,collapse = ','))
  }) %>% unlist()
}
del_col<-which(apply(st_step_CAF,2,function(x){
  length(which(x==''))
})==nrow(st_step_CAF))
if(length(del_col)>0) st_step_CAF<-st_step_CAF[,-del_col]

step_cell<-c()
for(ss in 5:ncol(st_step_CAF)){##ss=9
  step_spot<-unlist(strsplit(st_step_CAF[,ss],',')) %>% unique()
  
  if(length(intersect(rownames(st_Deco),step_spot))>0){
    step_Deco<-st_Deco[intersect(rownames(st_Deco),step_spot),]
    #add_cell<-apply(step_Deco,2,mean)
    add_cell<-reshape2::melt(as.matrix(step_Deco))
    add_cell$step<-paste0('step_',(ss-4))
    step_cell<-rbind(step_cell,add_cell)
  }
  
}

write.table(step_cell,paste0(dir_out,'public_RCTD_step_r.txt'),quote = F,sep = '\t',row.names = F)

#######NR组
st_step<-hop_infor_imm_nr[,c("X","imagerow","imagecol","FinalLocalType",paste0('step_',1:10))]
#st_step<-st_step[rownames(rctd_data_r),]
rownames(st_step) <- st_step$X
st_step <- st_step[,-1]
#st_step<-st_step[which(st_step$FinalLocalType=='Normal'|st_step$FinalLocalType=='Immune'),]
st_Deco<-rctd_data_nr[rownames(st_step)[st_step$FinalLocalType%in%c('Normal','Immune')],]
st_Deco <- st_Deco[,-1]
st_Deco <- st_Deco[,-c(13:16)]


# 提取每行的前缀（从行名中获取）
row_prefixes <- gsub("_[^_]+$", "", rownames(st_step))
#cat("行名前缀:\n")
#print(row_prefixes)

# 方法1: 为每个步骤补全细胞名
st_step_full <- st_step

# 处理step_1到step_10列
step_cols <- grep("^step_", colnames(st_step), value = TRUE)

for(ss in step_cols) {
  cat("\n处理", ss, "\n")
  
  # 为每一行补全细胞名
  full_names <- sapply(1:nrow(st_step), function(i) {
    current_cells <- as.character(st_step[i, ss])
    if(is.na(current_cells) || current_cells == "") return("")
    
    # 分割细胞名
    cells <- unlist(strsplit(current_cells, ','))
    cells <- cells[!is.na(cells) & cells != ""]
    
    if(length(cells) == 0) return("")
    
    # 补全为完整格式
    full_cells <- paste0(row_prefixes[i], "_", cells)
    return(paste(full_cells, collapse = ','))
  })
  
  st_step_full[, ss] <- full_names
}

write.table(st_step_full,file = paste0(dir_out,'public_step10_nr.txt'),quote = F,sep = '\t',row.names = T)

st_step_CAF<-st_step_full[st_step_full$FinalLocalType%in%c('Boundary','Dispersion'),]

for(jj in ncol(st_step_CAF):6){#jj=6
  y_spot<-unlist(strsplit(st_step_CAF[,(jj-1)],',')) %>% unique()
  st_step_CAF[,jj]<-lapply(1:nrow(st_step_CAF),function(x){##x=1
    x_spot<-strsplit(st_step_CAF[x,jj],',') %>% unlist()
    x_spot<-setdiff(x_spot,y_spot)
    return(paste0(x_spot,collapse = ','))
  }) %>% unlist()
}
del_col<-which(apply(st_step_CAF,2,function(x){
  length(which(x==''))
})==nrow(st_step_CAF))
if(length(del_col)>0) st_step_CAF<-st_step_CAF[,-del_col]

step_cell<-c()
for(ss in 5:ncol(st_step_CAF)){##ss=9
  step_spot<-unlist(strsplit(st_step_CAF[,ss],',')) %>% unique()
  
  if(length(intersect(rownames(st_Deco),step_spot))>0){
    step_Deco<-st_Deco[intersect(rownames(st_Deco),step_spot),]
    #add_cell<-apply(step_Deco,2,mean)
    add_cell<-reshape2::melt(as.matrix(step_Deco))
    add_cell$step<-paste0('step_',(ss-4))
    step_cell<-rbind(step_cell,add_cell)
  }
  
}

write.table(step_cell,paste0(dir_out,'public_RCTD_step_nr.txt'),quote = F,sep = '\t',row.names = F)


#######NR_TIB组
st_step<-hop_infor_imm_nr_TIB[,c("X","imagerow","imagecol","FinalLocalType",paste0('step_',1:10))]
#st_step<-st_step[rownames(rctd_data_r),]
rownames(st_step) <- st_step$X
st_step <- st_step[,-1]
#st_step<-st_step[which(st_step$FinalLocalType=='Normal'|st_step$FinalLocalType=='Immune'),]
st_Deco<-rctd_data_nr_tib[rownames(st_step)[st_step$FinalLocalType%in%c('Normal','Immune')],]
st_Deco <- st_Deco[,-1]
st_Deco <- st_Deco[,-c(13:16)]


# 提取每行的前缀（从行名中获取）
row_prefixes <- gsub("_[^_]+$", "", rownames(st_step))
#cat("行名前缀:\n")
#print(row_prefixes)

# 方法1: 为每个步骤补全细胞名
st_step_full <- st_step

# 处理step_1到step_10列
step_cols <- grep("^step_", colnames(st_step), value = TRUE)

for(ss in step_cols) {
  cat("\n处理", ss, "\n")
  
  # 为每一行补全细胞名
  full_names <- sapply(1:nrow(st_step), function(i) {
    current_cells <- as.character(st_step[i, ss])
    if(is.na(current_cells) || current_cells == "") return("")
    
    # 分割细胞名
    cells <- unlist(strsplit(current_cells, ','))
    cells <- cells[!is.na(cells) & cells != ""]
    
    if(length(cells) == 0) return("")
    
    # 补全为完整格式
    full_cells <- paste0(row_prefixes[i], "_", cells)
    return(paste(full_cells, collapse = ','))
  })
  
  st_step_full[, ss] <- full_names
}

write.table(st_step_full,file = paste0(dir_out,'public_step10_nr_tib.txt'),quote = F,sep = '\t',row.names = T)

st_step_CAF<-st_step_full[st_step_full$FinalLocalType%in%c('Boundary','Dispersion'),]

for(jj in ncol(st_step_CAF):6){#jj=6
  y_spot<-unlist(strsplit(st_step_CAF[,(jj-1)],',')) %>% unique()
  st_step_CAF[,jj]<-lapply(1:nrow(st_step_CAF),function(x){##x=1
    x_spot<-strsplit(st_step_CAF[x,jj],',') %>% unlist()
    x_spot<-setdiff(x_spot,y_spot)
    return(paste0(x_spot,collapse = ','))
  }) %>% unlist()
}
del_col<-which(apply(st_step_CAF,2,function(x){
  length(which(x==''))
})==nrow(st_step_CAF))
if(length(del_col)>0) st_step_CAF<-st_step_CAF[,-del_col]

step_cell<-c()
for(ss in 5:ncol(st_step_CAF)){##ss=9
  step_spot<-unlist(strsplit(st_step_CAF[,ss],',')) %>% unique()
  
  if(length(intersect(rownames(st_Deco),step_spot))>0){
    step_Deco<-st_Deco[intersect(rownames(st_Deco),step_spot),]
    #add_cell<-apply(step_Deco,2,mean)
    add_cell<-reshape2::melt(as.matrix(step_Deco))
    add_cell$step<-paste0('step_',(ss-4))
    step_cell<-rbind(step_cell,add_cell)
  }
  
}

write.table(step_cell,paste0(dir_out,'public_RCTD_step_nr_tib.txt'),quote = F,sep = '\t',row.names = F)

file_RCTDstep<-list.files(pattern = 'public_RCTD_step',path = dir_step,recursive = T)
patient <- gsub("\\.txt$", "", file_RCTDstep)


all_RCTD_step<-c()
for(i in c(1,3)){
  RCTD_step<-read.delim(paste0(dir_step,file_RCTDstep[i]),stringsAsFactors = F,check.names = F)
  RCTD_step$patient<-patient[i]
  all_RCTD_step<-rbind(all_RCTD_step,RCTD_step)
}

###nr vs r
#### CAF 比较
p_data<-all_RCTD_step[all_RCTD_step$Var2%in%"CAF",]
table(p_data$step,p_data$patient)

p_data<-p_data[p_data$step%in%paste0('step_',1:6),]
p_box<-aggregate(p_data$value,by=list(p_data$patient,p_data$step),mean)
colnames(p_box)<-c('patient','step','value')
p_box<-p_box[order(p_box$patient),]

library(ggpubr)
p7 <- ggpaired(p_box, 
               x = "patient", 
               y = "value",
               color = "patient", 
               line.color = "gray", 
               line.size = 0.4,
               palette = "jco") +
  stat_compare_means(paired = TRUE)

p7 <- ggplot(p_box, aes(x = patient, y = value, fill = patient)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = patient), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 1, outlier.shape = NA) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("public_RCTD_step_nr" = "#15629e", "public_RCTD_step_r" = "#ddb424")) +
  scale_color_manual(values = c("public_RCTD_step_nr" = "#15629e", "public_RCTD_step_r" = "#ddb424")) +
  #ggtitle(paste0('ImmLigand_', i)) +
  # 添加p值（使用非配对检验，因为数据不平衡）
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",  # 使用Wilcoxon秩和检验（非配对）
                     label.x = 1.5,
                     label.y = 0.17) +  # 设置p值位置在y=1.8
  # 设置y轴范围为0-2
  coord_cartesian(ylim = c(0.125, 0.175)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p7)

pdf(paste0(dir_step,'CAF_compare_public1.pdf'),width = 5,height = 4)
print(p7)
dev.off()



##ESCC
dir_step<-'/ESCC/'
sample_group<-'public'
file_RCTDstep<-list.files(pattern = paste0('_RCTD_step_',sample_group,'.txt'),path = dir_step,recursive = T)
patient<-unlist(lapply(strsplit(file_RCTDstep,'/'),function(x) x[1]))


all_RCTD_step<-c()
for(i in 1:2){
  RCTD_step<-read.delim(paste0(dir_step,file_RCTDstep[i]),stringsAsFactors = F,check.names = F)
  RCTD_step$patient<-patient[i]
  all_RCTD_step<-rbind(all_RCTD_step,RCTD_step)
}


#### CAF 比较
p_data<-all_RCTD_step[all_RCTD_step$Var2%in%"CAF",]
table(p_data$step,p_data$patient)

p_data<-p_data[p_data$step%in%paste0('step_',1:6),]
p_box<-aggregate(p_data$value,by=list(p_data$patient,p_data$step),mean)
colnames(p_box)<-c('patient','step','value')
p_box<-p_box[order(p_box$patient),]


p8 <- ggplot(p_box, aes(x = patient, y = value, fill = patient)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = patient), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 1, outlier.shape = NA) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("BGM" = "#15629e", "XYZ" = "#ddb424")) +
  scale_color_manual(values = c("BGM" = "#15629e", "XYZ" = "#ddb424")) +
  #ggtitle(paste0('ImmLigand_', i)) +
  # 添加p值（使用非配对检验，因为数据不平衡）
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",  # 使用Wilcoxon秩和检验（非配对）
                     label.x = 1.5,
                     label.y = 0.3) +  # 设置p值位置在y=1.8
  # 设置y轴范围为0-2
  coord_cartesian(ylim = c(0, 0.35)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p8)

pdf(paste0(dir_step,'CAF_compare_ESCC1.pdf'),width = 5,height = 4)
print(p8)
dev.off()
