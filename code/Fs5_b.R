###MP得分结果
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)



###数据来自357_newST_NMF_MPscore.R
dir_MP<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/'
file_MP<-list.files(pattern = '_MP_score.txt',path = dir_MP,recursive = T)
length(file_MP)
file_MP[1:10]
dataSlice<-unlist(lapply(strsplit(file_MP,'_MP_score'),function(x)x[1]))
dataSlice[1:10]
grep('OVCA_GSE189843_GSM5708493',gsub('/','_',dataSlice))

for(i in 1:length(dataSlice)){
  # i=496
  st_MP<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_MP[1:3,]
  MP_top2<-apply(st_MP[,1:17],1,function(x){
    x<-x[order(x,decreasing = T)]
    return(names(x)[1:2])
  }) %>% t()
  st_MP$MP_top1<-MP_top2[,1]
  st_MP$MP_top2<-MP_top2[,2]
  
  write.table(st_MP,paste0(dir_MP,file_MP[i]),quote = F,sep = '\t')
}


####绘制散点图##################
all_MP<-c()
for(i in 1:length(dataSlice)){
  # i=1
  st_MP<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_MP[1:3,]
  
  st_MP$cell_name<-rownames(st_MP)
  st_MP$dataSlice<-gsub('/','_',dataSlice[i])
  
  all_MP<-rbind(all_MP,st_MP)
  print(dataSlice[i])
}
all_MP[1:3,]
rownames(all_MP)<-paste0(all_MP$dataSlice,'_',all_MP$cell_name)
table(rownames(order_all_umap)==rownames(all_MP))
all_MP2<-all_MP[rownames(order_all_umap),]
table(rownames(order_all_umap)==rownames(all_MP2))
write.table(all_MP,'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/all_MP.txt',
            quote = F,sep = '\t')


all_MP<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/new_st_NMF/all_MP.txt',
                   stringsAsFactors = F,check.names = F)
all_MP<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/all_MP.txt',
                   stringsAsFactors = F,check.names = F)
rownames(all_MP)<-paste0(all_MP$dataSlice,'_',all_MP$cell_name)

all_umap<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/ST_merge_res0.8_metadatainfor_BdyCoreBud.txt',
                     stringsAsFactors = F,check.names = F)
all_umap<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/ST_merge_res0.8_metadatainfor_BdyCoreBud.txt',
                     stringsAsFactors = F,check.names = F)
length(interaction(rownames(all_umap),rownames(all_MP)))
all_MP<-all_MP[rownames(all_umap),]


length(intersect(unique(all_umap$batch),gsub('/','_',dataSlice)))
order_all_umap<-c()
for(i in 1:length(dataSlice)){
  #i=1
  st_MP<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_MP[1:3,]
  sub_all_umap<-all_umap[all_umap$batch%in%gsub('/','_',dataSlice[i]),]
  rownames(sub_all_umap)<-sub_all_umap$cell_name
  sub_all_umap<-sub_all_umap[rownames(st_MP),]
  sub_all_umap[1:3,]
  rownames(sub_all_umap)<-paste0(sub_all_umap$batch,'_',sub_all_umap$cell_name)
  order_all_umap<-rbind(order_all_umap,sub_all_umap)
  
  print(dataSlice[i])
}
order_all_umap[1:3,]
write.table(order_all_umap,'/data/zhouweiwei/ST_analysis/data/10X_Visium/all_umap_order.txt',
            quote = F,sep = '\t')


all_umap<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/all_umap_order.txt',
                     stringsAsFactors = F,check.names = F)
all_MP<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/all_MP.txt',
                   stringsAsFactors = F,check.names = F)
rownames(all_MP)<-paste0(all_MP$dataSlice,'_',all_MP$cell_name)
table(rownames(all_umap)==rownames(all_MP))
setdiff(rownames(all_umap),rownames(all_MP))
setdiff(rownames(all_MP),rownames(all_umap))

all_umap$MP_top1<-all_MP$MP_top1
write.table(all_umap,'/data/zhouweiwei/ST_analysis/data/10X_Visium/all_umap_order.txt',
            quote = F,sep = '\t')


###开始画图
all_umap<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/all_umap_order.txt',
                     stringsAsFactors = F,check.names = F)
all_MP<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/all_MP.txt',
                   stringsAsFactors = F,check.names = F)

# all_umap<-order_all_umap
all_umap$MP_top1<-all_MP$MP_top1

cancer_col<-read.delim('E:/Mirror/ST_analysis/data/35cancer_color.txt')
cancer_col<-read.delim('/data/zhouweiwei/ST_analysis/data/35cancer_color.txt')
col2<-cancer_col$color
names(col2)<-cancer_col$cancer

spot_mal<-c('Core','Boundary','Budding')
plot_data<-all_umap[all_umap$FinalLocalType%in%spot_mal,]
p1<-ggplot(data = plot_data,aes(x = UMAP_1, y = UMAP_2,color=MP_top1)) + 
  geom_point(size=0.2)+
  scale_color_manual(values = c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1",
                                "MP_5"="#6A8EC9","MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861",
                                "MP_9"="#652884","MP_10"="#444577","MP_11"="#8A7355","MP_12"="#B3BB61",
                                "MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
                                "MP_16"="#399335",'MP_17'='#96DD88'))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggtitle('assign_MP_top1')

# pdf(paste0(dir_pic,'assign_MP.pdf'),width = 8,height = 6)
dir_pic<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/00pic/'
png(paste0(dir_pic,'assign_MP.png'),width = 540, height = 480)
print(p1)
dev.off()



####MP丰度环形条状图
plot_data<-all_umap[all_umap$FinalLocalType%in%spot_mal,]
plot_data2<-as.data.frame.array(table(plot_data[,'MP_top1']))
plot_data2$MP_top1<-rownames(plot_data2)
colnames(plot_data2)[1]<-'num'
plot_data2$value<-plot_data2$num/sum(plot_data2$num)
#plot_data2<-plot_data2[order(plot_data2$value,decreasing = T),]
plot_data2<-mutate(plot_data2,MP_top1 = factor(plot_data2$MP_top1, levels = c(paste0("MP_",1:17))))

p5 <- ggplot(plot_data2,aes(x = 1, y = value, fill = MP_top1)) +
  geom_col(colour = "white")+ 
  coord_polar(theta = "y", start = 2) +
  geom_text(aes(label = paste0(round(value * 100, 2), "%")),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1",
                             "MP_5"="#6A8EC9","MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861",
                             "MP_9"="#652884","MP_10"="#444577","MP_11"="#8A7355","MP_12"="#B3BB61",
                             "MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
                             "MP_16"="#399335",'MP_17'='#96DD88'))+
  xlim(c(-1, 2)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
print(p5)
dir_pic<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/00pic/'
pdf(paste0(dir_pic,'assign_MP_pie.pdf'),width = 7,height = 7)
# png(paste0(dir_pic,'assign_MP_pie.png'),width = 520, height = 480)
print(p5)
dev.off()



