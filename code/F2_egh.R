########LUAD进展: AIS MIA IAC
####绘制stem得分在各阶段切片中的箱式图
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gghalves)

dir_locationType<-"/new_st_copykat/LUAD/GSE307534/"
dir_stemness <- "/new_st_CytoTRACE/LUAD/GSE307534/"
file_Type<-list.files(pattern = "BdyTumorCore.txt",path = dir_locationType,recursive = T)
file_stem<-list.files(pattern = "Stemness.txt",path = dir_stemness,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Type,"_"),function(x) x[1]))

pro_data<-c()
for(i in 1:length(dataset_slice)){
  #i <- 1
  st_Type<-read.delim(paste0(dir_locationType,file_Type[i]),stringsAsFactors = F,check.names = F)
  #st_Type<-st_Type[which(st_Type$FinalLocalType=="Core"),]
  st_stem<-read.delim(paste0(dir_stemness,file_stem[i]),stringsAsFactors = F,check.names = F)
  st_stem<-st_stem[st_Type$cell_name,]
  st_stem<-data.frame(stage=rep(dataset_slice[i],nrow(st_stem)),
                      stem_score=st_stem$score,
                      localType=st_Type$FinalLocalType)
  
  pro_data<-rbind(pro_data,st_stem)
  
}
table(pro_data$stage)

write.table(pro_data,file = "/stemness.txt",sep = "\t",col.names = T,row.names = F)

pro_data <- read.table(file = "/stemness.txt",sep = "\t",header = T)
stage_data <- read.table("/stage_data.txt",header = T,sep = "\t")
pro_data <- merge(pro_data,stage_data,by.x = "stage",by.y = "GSM")
pro_data <- pro_data[c(1,2,3,5,6)]
colnames(pro_data)[4] <- "stage2"

table(pro_data$stage,pro_data$stage2)


##固定顺序
spot_type<-c("AIS", "MIA", "IAC")
local_type<-c('Core','Boundary','Budding')

###选切片
plot_data<-pro_data[pro_data$stage2%in%spot_type,]
plot_data<-plot_data[plot_data$localType%in%local_type,]

plot_data<-mutate(plot_data,stage2 = factor(plot_data$stage2, levels = spot_type))
plot_data<-mutate(plot_data,localType = factor(plot_data$localType, levels = local_type))

e <- ggplot(plot_data, aes(x = stage, y = stem_score,fill=localType))+ 
  geom_violin(aes(color = localType), trim = T,position = position_dodge(0.8),alpha=0.6) +
  geom_boxplot(aes(color = localType), width = 0.2,position = position_dodge(0.8),alpha=0.5)+
  theme_classic(base_size = 10)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                               "Budding"="#EE762D","non_Mal"="#5477AF"))+
  scale_color_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                                "Budding"="#EE762D","non_Mal"="#5477AF"))+
  ggtitle('stem of LUAD progress')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12))+
  #ylim(c(0,9))+
  #xlab("")+ylab("CytoTRACE")+
  # geom_signif(comparisons = list(c("Boundary","Core"),
  #                                c("Dispersion","Core"),
  #                                c("Dispersion","Boundary")),
  # map_signif_level=F,
  # textsize=4,test=wilcox.test,step_increase=.2,test.args = "two.sided")+
  stat_compare_means(label = "p.format",label.x = 2,label.y = 1.1,size = 2)+
  theme_bw()+
  theme(axis.title.x = element_text(size=12),
        axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y =element_text(size=12),
        axis.text.y = element_text(size = 10) )
print(e)
pdf('/slice_progress2_CoreBdyDis.pdf',width = 14, height = 4)
print(e)
dev.off()



plot_data2<-plot_data
plot_data2$stage2 <- as.character(plot_data2$stage2)
# 赋值
plot_data2[plot_data2$stage2 %in% c("AIS","MIA"), 4] <- "AIS_MIA"
#再转回因子
plot_data2$stage2 <- as.factor(plot_data2$stage2)
table(plot_data2$stage2)

e <- ggplot(plot_data2, aes(x = stage2, y = stem_score,fill=localType))+ 
  geom_violin(aes(color = localType), trim = T,position = position_dodge(0.8),alpha=0.6) +
  geom_boxplot(aes(color = localType), width = 0.2,position = position_dodge(0.8),alpha=0.5)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#FF9900","Core"="#990033",
                               "Budding"="#CD5A5A","non_Mal"="#5477AF"))+
  scale_color_manual(values = c("Boundary"="#FF9900","Core"="#990033",
                                "Budding"="#CD5A5A","non_Mal"="#5477AF"))+
  ggtitle('stem of LUAD progress: AIS_MIA IAC')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12))+
  #ylim(c(0,9))+
  #xlab("")+ylab("CytoTRACE")+
  # geom_signif(comparisons = list(c("Boundary","Core"),
  #                                c("Dispersion","Core"),
  #                                c("Dispersion","Boundary")),
  # map_signif_level=F,
  # textsize=4,test=wilcox.test,step_increase=.2,test.args = "two.sided")+
  stat_compare_means(label.x = 2,label.y = 1.1)+
  theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(e)

pdf('/12slices_AISMIA_IAC_CoreBdyDis.pdf',width = 7, height = 5)
print(e)
dev.off()


###d数据
pro_data1 <- read.table(file = "/all_slice_stemScore.txt",sep = "\t",header = T)
old_luad_slices <- c("luad01/slice1","luad01/slice2","luad02/slice1","luad02/slice2","luad03/slice1","luad03/slice2")
old_plot_data<-pro_data1[pro_data1$dataset_slice%in%old_luad_slices,]
old_plot_data <- old_plot_data[,c(1,4,3)]
old_plot_data$stage2 <- NA
old_plot_data[old_plot_data$dataset_slice %in% c("luad01/slice1","luad01/slice2","luad02/slice1","luad02/slice2"), 4] <- "AIS_MIA"
old_plot_data[old_plot_data$dataset_slice %in% c("luad03/slice1","luad03/slice2"), 4] <- "IAC"


###d数据+f数据: 18slices
new_plot_data <- plot_data2[,c(1:4)]
colnames(new_plot_data) <- colnames(old_plot_data)
all_plot_data <- rbind(old_plot_data,new_plot_data)

all_plot_data[all_plot_data$stage2 %in% c("LUAD"), 4] <- "IAC"

##固定顺序
spot_type<-c("AIS_MIA", "IAC")
local_type<-c('Core','Boundary','Budding')

all_plot_data<-all_plot_data[all_plot_data$LocalType%in%local_type,]

all_plot_data<-mutate(all_plot_data,stage2 = factor(all_plot_data$stage2, levels = spot_type))
all_plot_data<-mutate(all_plot_data,LocalType = factor(all_plot_data$LocalType, levels = local_type))

e1 <- ggplot(all_plot_data, aes(x = stage2, y = CytoTRACE,fill=LocalType))+ 
  geom_violin(aes(color = LocalType), trim = T,position = position_dodge(0.8),alpha=0.6) +
  geom_boxplot(aes(color = LocalType), width = 0.2,position = position_dodge(0.8),alpha=0.5)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#FF9900","Core"="#990033",
                               "Budding"="#CD5A5A","non_Mal"="#5477AF"))+
  scale_color_manual(values = c("Boundary"="#FF9900","Core"="#990033",
                                "Budding"="#CD5A5A","non_Mal"="#5477AF"))+
  ggtitle('stem of LUAD progress: AIS_MIA IAC')+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(size = 12))+
  #ylim(c(0,9))+
  #xlab("")+ylab("CytoTRACE")+
  # geom_signif(comparisons = list(c("Boundary","Core"),
  #                                c("Dispersion","Core"),
  #                                c("Dispersion","Boundary")),
  # map_signif_level=F,
  # textsize=4,test=wilcox.test,step_increase=.2,test.args = "two.sided")+
  stat_compare_means(label.x = 2,label.y = 1.1)+
  theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(e1)

pdf('/18slices_AISMIA_IAC_CoreBdyDis.pdf',width = 7, height = 5)
print(e1)
dev.off()


