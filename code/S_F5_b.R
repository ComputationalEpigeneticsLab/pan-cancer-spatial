####对每个spot分配MP
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)


integration_data<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/meta_order.txt',stringsAsFactors = F,check.names = F)
integration_data$LocationType[which(integration_data$LocationType=='Tumor')]<-'Dispersion'
integration_data$LocationType[which(integration_data$LocationType=='diploid')]<-'Normal'
integration_umap<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/umap.txt',stringsAsFactors = F,check.names = F)
integration_umap<-integration_umap[rownames(integration_data),]
integration_umap$slice<-integration_data$slice
integration_umap$cell_name<-integration_data$cell_name
integration_umap$LocationType<-integration_data$LocationType

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_bdy,'_Bdy'),function(x) x[1]))


dir_MPscore<-'E:/Mirror/ST_analysis/data/10X Visium/NMF_module/moduleScore/'
file_MPscore<-list.files(pattern = 'txt',path = dir_MPscore,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_MPscore,'_MP'),function(x) x[1]))

all_assign_MP<-c()
all_MP_score<-c()
for(i in 1:length(file_MPscore)){
  #i=1
  #st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  
  MP_score<-read.delim(paste0(dir_MPscore,file_MPscore[i]),stringsAsFactors = F,check.names = F)
  MP_score<-MP_score[,-1]
  rownames(MP_score)<-paste0(dataSlice[i],'_',rownames(MP_score))
  #all_MP_score<-rbind(all_MP_score,MP_score)
  MP_score<-scale(t(MP_score),center = T,scale = F) %>% t()
  assign_MP<-apply(MP_score,1,function(x){#x=MP_score[1,]
    names(x)<-colnames(MP_score)
    x<-x[order(x,decreasing = T)]
    y<-'unresolved'
    if(x[1]*0.95>=x[2]) y<-names(x)[1]
    return(y)
  })
  table(assign_MP)
  assign_MP<-data.frame(assign_MP=assign_MP,slice=dataSlice[i])
  #rownames(assign_MP)<-paste0(assign_MP$slice,'_',rownames(assign_MP))
  
  all_assign_MP<-rbind(all_assign_MP,assign_MP)
  print(dataSlice[i])
}

table(all_assign_MP$assign_MP)

all_assign_MP1<-all_assign_MP[paste0(integration_umap$slice,'_',integration_umap$cell_name),]
all_MP_score1<-all_MP_score[paste0(integration_umap$slice,'_',integration_umap$cell_name),]

integration_umap$assign_MP95<-all_assign_MP1$assign_MP
integration_umap<-cbind(integration_umap,all_MP_score1)

write.table(integration_umap,'E:/Mirror/ST_analysis/data/pic_data/integration_umap_MP.txt',quote = F,sep = '\t')
integration_umap_MP<-read.delim('E:/Mirror/ST_analysis/data/pic_data/integration_umap_MP.txt',stringsAsFactors = F,check.names = F)
table(integration_umap_MP$LocationType)

###找出每个spot的第一和第二大的MP类型
MP_top2<-apply(integration_umap_MP[,7:19],1,function(x){
  x<-x[order(x,decreasing = T)]
  return(names(x)[1:2])
}) %>% t()
integration_umap_MP$MP_top1<-MP_top2[,1]
integration_umap_MP$MP_top2<-MP_top2[,2]
integration_umap_MP<-integration_umap_MP[,c(1:6,20:ncol(integration_umap_MP),7:19)]
colnames(integration_umap_MP)[6]<-'assign_MP85'
write.table(integration_umap_MP,'E:/Mirror/ST_analysis/data/pic_data/integration_umap_MP.txt',quote = F,sep = '\t')

table(integration_umap_MP$MP_top1)
un_MP<-integration_umap_MP[which(integration_umap_MP$assign_MP85=='unresolved'),]
table(un_MP$MP_top1)
table(un_MP$MP_top2)
table(paste0(un_MP$MP_top1,':',un_MP$MP_top2))


colnames(integration_umap_MP)[1:2]<-c('x','y')
table(integration_umap_MP$assign_MP)

dir_pic<-'E:/Mirror/ST_analysis/pic/re/MP_score/'

spot_mal<-c('Core','Boundary','Dispersion')
plot_data<-integration_umap_MP[integration_umap_MP$LocationType%in%spot_mal,]
p1<-ggplot(data = plot_data,aes(x = x, y = y,color=MP_top1)) + 
  geom_point(size=0.2)+
  scale_color_manual(values = c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
                                "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
                                "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
                                'MP_14'='#FF9694','unresolved'='#7F7F7F'))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggtitle('assign_MP_top1')

pdf(paste0(dir_pic,'assign_MP.pdf'),width = 8,height = 6)
print(p1)
dev.off()

pdf(paste0(dir_pic,'assign_MP_top1.pdf'),width = 8,height = 6)
print(p1)
dev.off()

####MP丰度环形条状图
table(integration_umap_MP$assign_MP95)
table(integration_umap_MP$assign_MP)

plot_data<-as.data.frame.array(table(integration_umap_MP[integration_umap_MP$LocationType%in%spot_mal,'MP_top1']))
plot_data$assign_MP95<-rownames(plot_data)
colnames(plot_data)[1]<-'num'
plot_data$value<-plot_data$num/sum(plot_data$num)
#plot_data<-plot_data[order(plot_data$value,decreasing = T),]
plot_data<-mutate(plot_data,assign_MP95 = factor(plot_data$assign_MP95, levels = c(paste0("MP_",2:14))))

p5 <- ggplot(plot_data,aes(x = 1, y = value, fill = assign_MP95)) +
  geom_col(colour = "white")+ 
  coord_polar(theta = "y", start = 1.65) +
  geom_text(aes(label = paste0(round(value * 100, 2), "%")),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
                             "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
                             "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
                             'MP_14'='#FF9694','unresolved'='#7F7F7F'))+
  xlim(c(-0.2, 2)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
print(p5)
pdf(paste0(dir_pic,'assign_MP_pie_top1.pdf'),width = 7,height = 7)
print(p5)
dev.off()