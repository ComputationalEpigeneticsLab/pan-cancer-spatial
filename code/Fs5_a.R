library(ggalluvial)
library(jsonlite)
library(plyr)
library(ggplot2)



new_MP<-readRDS("D:\\pan_cancer\\0_修稿\\F3\\F3a\\intra23_inter23Cluster_list(20-300).rds")
new_MP2<-unlist(new_MP)
new_MP2<-data.frame(cluster=rep(names(new_MP),as.numeric(unlist(lapply(new_MP,length)))),
                    slice=new_MP2)
table(new_MP2$cluster)


# program_MP<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/program_top50gene.txt',
#                        stringsAsFactors = F,check.names = F)
plot_data<-new_MP2
plot_data$cluster<-gsub('Cluster','MP',plot_data$cluster)
plot_data$cancer<-unlist(lapply(strsplit(plot_data$slice,'_'),function(x)x[1]))
table(plot_data$cancer)
length(unique(plot_data$cancer))
plot_data$program<-paste0(plot_data$cancer,'_program')
colnames(plot_data)[1]<-'MP'


####条形图
program_cancer<-as.data.frame.array(table(plot_data$cancer,plot_data$MP))
program_cancer<-reshape2::melt(as.matrix(program_cancer))
MP_order<-c('MP_1',
             'MP_2',
             'MP_3','MP_4','MP_5','MP_6',
             'MP_7', 
             'MP_8','MP_9', 'MP_10',  
             'MP_11', 'MP_12', 
             'MP_13',  
             'MP_14','MP_15','MP_16', 'MP_17')%>%toupper()
####需要排序的话改这些
program_cancer<-mutate(program_cancer,Legend = factor(program_cancer$Var2, levels = as.character(MP_order)))
table(program_cancer$Var1,program_cancer$Var2)
program_cancer<-program_cancer[which(program_cancer$value>0),]


cancer_col<-read.delim('D:\\pan_cancer\\0_修稿\\35cancer_color.txt')
cancer_col<-cancer_col[cancer_col$cancer%in%unique(plot_data$cancer),]
col2<-cancer_col$color
names(col2)<-cancer_col$cancer

p_compare<-ggplot(program_cancer,aes(x=Var2,y=value,fill=Var1)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  scale_fill_manual(values = col2)+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("MP")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('program_cancer')
print(p_compare)
pdf('D:\\pan_cancer\\0_修稿\\Fs5\\Fs5a\\Fs5a.pdf',height = 6,width = 9)
print(p_compare)
dev.off()



