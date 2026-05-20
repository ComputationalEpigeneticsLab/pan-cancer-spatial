##所有癌症的各区域干性绘制箱式图####
#1.每个癌症切片的cytotrace合并
base_dir="/new_st_CytoTRACE"
output_dir='/stemness/'
all_files <- list.files(base_dir,recursive = TRUE,full.names = TRUE,pattern = "Stemness.txt")
length(all_files)

score_list <- lapply(all_files, function(x) {
  read.table(x, header = TRUE, sep = "\t", row.names = 1)
})
combined_data <- do.call(rbind, score_list)
write.table(combined_data,file = paste0(output_dir,'all_slice_stem.txt'),row.names = T,sep = '\t',quote = F)

# #2.每个切片的区域分类注释信息合并
base_dir='/new_st_copykat'
output_dir='/stemness/'
all_files <- list.files(base_dir,recursive = TRUE,full.names = TRUE,pattern = "BdyCoreBud.txt")
length(all_files)
sample_list <- lapply(all_files, function(x) {
  aa=read.table(x, header = TRUE, sep = "\t", row.names = 1)
  aa$cancer=strsplit(x,"\\/")[[1]][8]
  sampleid=paste0(strsplit(x,"\\/")[[1]][9],"_",strsplit(x,"\\/")[[1]][10])
  aa$dataset_slice=gsub("\\.txt","",sampleid)
  return(aa[,c("cell_name","FinalLocalType","cancer","dataset_slice")])
})
combined_data <- do.call(rbind, sample_list)
colnames(combined_data)[1]='cellname'
write.table(combined_data,file = paste0(output_dir,'all_slice_spottype.txt'),row.names = T,sep = '\t',quote = F)

##可视化
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)

all_slice_stem = read.table('/stemness/all_slice_stem.txt',header = T,row.names = 1,sep = '\t');gc()
all_slice_spottype = read.table('/stemness/all_slice_spottype.txt',header = T,row.names = 1,sep = '\t')
table(rownames(all_slice_spottype)%in%rownames(all_slice_stem))
table(all_slice_spottype$FinalLocalType)
all_slice_spottype[all_slice_spottype$FinalLocalType%in%c('Budding'),'FinalLocalType']='Dispersion'
all_slice_spottype[all_slice_spottype$FinalLocalType%in%c('Immune','Normal'),'FinalLocalType']='Stromal'
all_slice_stem = all_slice_stem[match(rownames(all_slice_spottype),rownames(all_slice_stem)),]
identical(rownames(all_slice_spottype),rownames(all_slice_stem))
all_slice_spottype$score = all_slice_stem$score
rm(all_slice_stem);gc()

all_box<-all_slice_spottype
all_box$cancer<-'all'
plotdata<-rbind(all_slice_spottype,all_box)
plotdata$FinalLocalType = factor(plotdata$FinalLocalType,levels = c('Core','Boundary','Dispersion','Stromal'))

p<-ggplot(plotdata, aes(x = cancer, y = score,fill=FinalLocalType))+
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = FinalLocalType), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                               "Dispersion"="#EE762D","Stromal"="#a0d5e8"))+
  scale_color_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                                "Dispersion"="#EE762D","Stromal"="#a0d5e8"))+
  #ggtitle("Stemness")+
  theme(plot.title = element_text(hjust = 0.4))+
  theme(plot.title = element_text(size = 6))+
  #stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =1)+
  stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =1.01)+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
p
ggsave("allCancer_stemness_box1.pdf",p,device = "pdf",path = "/stemness",height = 2.2,width = 35)



























