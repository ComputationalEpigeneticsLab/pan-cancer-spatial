#####LUAD阶段干性比较 绘制散点柱状图
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gghalves)



dir_locationType<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"

file_Type<-list.files(pattern = "BdyTumorCore.txt",path = dir_locationType,recursive = T)
file_stem<-list.files(pattern = "Stemness.txt",path = dir_locationType,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_Type,"_"),function(x) x[1]))
file_Type[164:169]
file_stem[164:169]
dataset_slice[164:169]

pro_data<-c()
for(i in 164:169){#i=164
  st_Type<-read.delim(paste0(dir_locationType,file_Type[i]),stringsAsFactors = F,check.names = F)
  #st_Type<-st_Type[which(st_Type$FinalLocalType=="Core"),]
  st_stem<-read.delim(paste0(dir_locationType,file_stem[i]),stringsAsFactors = F,check.names = F)
  st_stem<-st_stem[st_Type$cell_name,]
  st_stem<-data.frame(stage=rep(dataset_slice[i],nrow(st_stem)),
                      stem_score=st_stem$CytoTRACE,
                      localType=st_Type$FinalLocalType)
  
  pro_data<-rbind(pro_data,st_stem)
  
}
table(pro_data$stage)
#pro_data$stage<-'luad'
pro_data$stage[grep('luad01',pro_data$stage)]<-"AIS"
pro_data$stage[grep('luad02',pro_data$stage)]<-"MIA"
pro_data$stage[grep('luad03',pro_data$stage)]<-"IAC"
table(pro_data$slice,pro_data$stage)


##固定顺序
spot_type<-c("AIS", "MIA", "IAC")
spot_type<-c("MIA", "IAC")
spot_type<-c("luad01/slice1", "luad01/slice2", "luad02/slice1",'luad02/slice2','luad03/slice1','luad03/slice2')

local_type<-c('Core','Boundary','Dispersion')
plot_data<-pro_data[pro_data$stage%in%spot_type,]
plot_data<-plot_data[plot_data$localType%in%local_type,]
plot_data<-mutate(plot_data,stage = factor(plot_data$stage, levels = spot_type))
plot_data<-mutate(plot_data,localType = factor(plot_data$localType, levels = local_type))

#?geom_violin
e <- ggplot(plot_data, aes(x = stage, y = stem_score,fill=localType))+ 
  geom_jitter(aes(color = localType),size = 1,alpha=0.6,
              position=position_jitterdodge(jitter.width = 0.5, 
                                            jitter.height = 0, 
                                            dodge.width = .8)) + # 不重叠的散点图
  stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8)) + # 误差棒，中位数，25%和75%分位数
  stat_summary(aes(fill = localType), fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8)) + # 中位数水平线
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                               "Dispersion"="#EE762D","non_Mal"="#5477AF"))+
  scale_color_manual(values = c("Boundary"="#F6B86D","Core"="#D62D28",
                                "Dispersion"="#EE762D","non_Mal"="#5477AF"))+
  ggtitle('Stemness')+
  theme(plot.title = element_text(hjust = 0.4))+
  theme(plot.title = element_text(size = 12))+
  stat_compare_means(label.x = 1,label.y=1.1 #,aes(label = ..p.signif..)
  )+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(e)
pdf('E:/Mirror/ST_analysis/pic/re/3/stem_luad_progress2_CoreBdyDis2.pdf',width = 7, height = 5)
print(e)
dev.off()

?stat_compare_means


wilcox.test(plot_data$stem_score[which(plot_data$stage=='MIA'&plot_data$localType=='Core')],
            plot_data$stem_score[which(plot_data$stage=='MIA'&plot_data$localType!='Core')],alternative = 'greater')[["p.value"]]
wilcox.test(plot_data$stem_score[which(plot_data$stage=='IAC'&plot_data$localType=='Core')],
            plot_data$stem_score[which(plot_data$stage=='IAC'&plot_data$localType!='Core')],alternative = 'greater')[["p.value"]]



















