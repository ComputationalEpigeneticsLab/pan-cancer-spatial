#####cellchat结果统计
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(reshape2)


####bdy与step1_CAF
dir_bdyCAF<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_bdyCAF,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for (i in 1:length(file_LR)) {
  #i=3
  LR_data<-read.csv(paste0(dir_bdyCAF,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    return('CAF'%in%a&&setdiff(unique(c(LR_data$source,LR_data$target)),'CAF')%in%a)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  
  if(nrow(LR_data)>0){
    LR_data$LR<-paste0(LR_data$ligand,'_of_',LR_data$receptor)
    path_LR<-as.data.frame(table(LR_data$pathway_name,LR_data$LR))
    path_LR<-path_LR[which(path_LR$Freq!=0),]
    path<-as.data.frame(table(LR_data$pathway_name))
    path_LR<-merge(path_LR,path,by='Var1',all=T)
    path_LR<-path_LR[,c(1,4,2,3)]
    colnames(path_LR)<-c('pathway_name','path_num','LR','LR_num')
    path_LR$slice<-dataSlice[i]
    all_path_LR<-rbind(all_path_LR,path_LR)
  }
  
  print(dataSlice[i])
}
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF.txt',quote = F,sep = '\t',row.names = F)



####bdy与step1_TAM
dir_bdyTAM<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_bdyTAM,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for (i in 1:length(file_LR)) {
  #i=3
  LR_data<-read.csv(paste0(dir_bdyTAM,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    return('TAM'%in%a&&setdiff(unique(c(LR_data$source,LR_data$target)),'TAM')%in%a)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  
  if(nrow(LR_data)>0){
    LR_data$LR<-paste0(LR_data$ligand,'_of_',LR_data$receptor)
    path_LR<-as.data.frame(table(LR_data$pathway_name,LR_data$LR))
    path_LR<-path_LR[which(path_LR$Freq!=0),]
    path<-as.data.frame(table(LR_data$pathway_name))
    path_LR<-merge(path_LR,path,by='Var1',all=T)
    path_LR<-path_LR[,c(1,4,2,3)]
    colnames(path_LR)<-c('pathway_name','path_num','LR','LR_num')
    path_LR$slice<-dataSlice[i]
    all_path_LR<-rbind(all_path_LR,path_LR)
  }
  
  print(dataSlice[i])
}
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM.txt',quote = F,sep = '\t',row.names = F)





####基质中CAF与Immune
dir_matrixCAF<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_CAF/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_matrixCAF,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for (i in 1:length(file_LR)) {
  #i=1
  LR_data<-read.csv(paste0(dir_matrixCAF,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    other<-lapply(setdiff(unique(c(LR_data$source,LR_data$target)),'CAF'),function(y){
      return(y%in%a)
    }) %>% unlist()
    return('CAF'%in%a&&TRUE%in%other)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  
  if(nrow(LR_data)>0){
    LR_data$LR<-paste0(LR_data$ligand,'_of_',LR_data$receptor)
    path_LR<-as.data.frame(table(LR_data$pathway_name,LR_data$LR))
    path_LR<-path_LR[which(path_LR$Freq!=0),]
    path<-as.data.frame(table(LR_data$pathway_name))
    path_LR<-merge(path_LR,path,by='Var1',all=T)
    path_LR<-path_LR[,c(1,4,2,3)]
    colnames(path_LR)<-c('pathway_name','path_num','LR','LR_num')
    path_LR$slice<-dataSlice[i]
    all_path_LR<-rbind(all_path_LR,path_LR)
  }
  
  print(dataSlice[i])
}
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF.txt',quote = F,sep = '\t',row.names = F)




####基质中TAM与Immune
dir_matrixTAM<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_matrixTAM,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for (i in 1:length(file_LR)) {
  #i=1
  LR_data<-read.csv(paste0(dir_matrixTAM,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    other<-lapply(setdiff(unique(c(LR_data$source,LR_data$target)),'TAM'),function(y){
      return(y%in%a)
    }) %>% unlist()
    return('TAM'%in%a&&TRUE%in%other)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  
  if(nrow(LR_data)>0){
    LR_data$LR<-paste0(LR_data$ligand,'_of_',LR_data$receptor)
    path_LR<-as.data.frame(table(LR_data$pathway_name,LR_data$LR))
    path_LR<-path_LR[which(path_LR$Freq!=0),]
    path<-as.data.frame(table(LR_data$pathway_name))
    path_LR<-merge(path_LR,path,by='Var1',all=T)
    path_LR<-path_LR[,c(1,4,2,3)]
    colnames(path_LR)<-c('pathway_name','path_num','LR','LR_num')
    path_LR$slice<-dataSlice[i]
    all_path_LR<-rbind(all_path_LR,path_LR)
  }
  
  print(dataSlice[i])
}
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM.txt',quote = F,sep = '\t',row.names = F)




####统计出现最多的通路及其对应的LR
bdy_step1CAF<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF.txt',stringsAsFactors = F,check.names = F)
bdy_step1TAM<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM.txt',stringsAsFactors = F,check.names = F)

matrix_CAF<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF.txt',stringsAsFactors = F,check.names = F)
matrix_TAM<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM.txt',stringsAsFactors = F,check.names = F)


bdy_step1CAFpath<-as.data.frame(table(bdy_step1CAF$pathway_name))####不能用table计算频率 重复出现的LR会导致pathway少
bdy_step1CAFpath<-bdy_step1CAFpath[order(bdy_step1CAFpath$Freq,decreasing = T),]
bdy_step1CAFLR<-bdy_step1CAF[bdy_step1CAF$pathway_name%in%bdy_step1CAFpath$Var1[1],]

bdy_step1TAMpath<-as.data.frame(table(bdy_step1TAM$pathway_name))
bdy_step1TAMpath<-bdy_step1TAMpath[order(bdy_step1TAMpath$Freq,decreasing = T),]
bdy_step1TAMLR<-bdy_step1TAM[bdy_step1TAM$pathway_name%in%bdy_step1TAMpath$Var1[1],]


matrixCAFpath<-as.data.frame(table(matrix_CAF$pathway_name))
matrixCAFpath<-matrixCAFpath[order(matrixCAFpath$Freq,decreasing = T),]
matrixCAFLR<-matrix_CAF[matrix_CAF$pathway_name%in%matrixCAFpath$Var1[1],]

matrixTAMpath<-as.data.frame(table(matrix_TAM$pathway_name))
matrixTAMpath<-matrixTAMpath[order(matrixTAMpath$Freq,decreasing = T),]
matrixTAMLR<-matrix_TAM[matrix_TAM$pathway_name%in%matrixTAMpath$Var1[1],]



dir_pic<-'E:/Mirror/ST_analysis/pic/cellchat/'


#######################################################################################################################
###网络图#########################################
library(ggraph)
library(tidygraph)
library(igraph)
library(dplyr)

library(svglite)
svglite("I:/2023课题/3.细菌呈递的/12.browse/3.细菌病毒图片/Bacteria_tree.svg", width = 12, height = 12)
g_bacteria
dev.off()

dir_pic<-'E:/Mirror/ST_analysis/pic/cellchat/'
# library(showtext)
# library(extrafont)
# font_import("C:/Windows/Fonts")
# font_add('TNR','C:/Windows/Fonts/arial.ttf')
# showtext_auto()


bdy_step1CAFLR<-bdy_step1CAF[bdy_step1CAF$pathway_name%in%'LAMININ',]##'LAMININ' 'COLLAGEN'
net_data<-aggregate(bdy_step1CAFLR$LR_num,by=list(bdy_step1CAFLR$LR),sum)
colnames(net_data)[2]<-'count'
net_data$ligand<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[1]))
net_data$receptor<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[2]))

write.table(net_data,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_CAF_net_data_LAMININ.txt',quote = F,sep = '\t',row.names = F)
net_data<-net_data[,c(3,4,2,1)]
write.csv(net_data,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_CAF_net_data_LAMININ.csv',quote = F,row.names = F)


net_data2<-data.frame(id=c(unique(net_data$ligand),unique(net_data$receptor)),
                      label=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))
write.csv(net_data2,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_CAF_net_data_LAMININ2.csv',quote = F,row.names = F)

net_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_CAF_net_data_COLLAGEN.txt',stringsAsFactors = F,check.names = F)
##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor","count")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)

#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1, 3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=1.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('bdy_step1CAF_COLLAGEN')
p_net
svglite(paste0(dir_pic,'2/net_bdy_step1CAF_COLLAGEN.svg'),width = 3.5,height = 2.6)
print(p_net)
dev.off()
?theme_graph


bdy_step1TAMLR<-bdy_step1TAM[bdy_step1TAM$pathway_name%in%'COLLAGEN',]##'LAMININ' 'COLLAGEN'
net_data<-aggregate(bdy_step1TAMLR$LR_num,by=list(bdy_step1TAMLR$LR),sum)
colnames(net_data)[2]<-'count'
net_data$ligand<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[1]))
net_data$receptor<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[2]))

write.table(net_data,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_TAM_net_data_COLLAGEN.txt',quote = F,sep = '\t',row.names = F)

net_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_TAM_net_data_LAMININ.txt',stringsAsFactors = F,check.names = F)
##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor","count")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)

#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1,3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=1.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('bdy_step1TAM_LAMININ')
p_net
svglite(paste0(dir_pic,'2/net_bdy_step1TAM_LAMININ.svg'),width = 3.5,height = 2.6)
print(p_net)
dev.off()



matrix_CAFLR<-matrix_CAF[matrix_CAF$pathway_name%in%'LAMININ',]##'LAMININ' 'COLLAGEN'
net_data<-aggregate(matrix_CAFLR$LR_num,by=list(matrix_CAFLR$LR),sum)
colnames(net_data)[2]<-'count'
net_data$ligand<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[1]))
net_data$receptor<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[2]))

write.table(net_data,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_net_data_LAMININ.txt',quote = F,sep = '\t',row.names = F)

net_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_net_data_COLLAGEN.txt',stringsAsFactors = F,check.names = F)

##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor","count")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)

#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'   'fabric'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1,3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=1.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('matrix_CAF_COLLAGEN')
p_net
svglite(paste0(dir_pic,'2/net_matrix_CAF_COLLAGEN.svg'),width = 3.5,height = 2.6)
print(p_net)
dev.off()




matrix_TAMLR<-matrix_TAM[matrix_TAM$pathway_name%in%'COLLAGEN',]##'LAMININ' 'COLLAGEN'
net_data<-aggregate(matrix_TAMLR$LR_num,by=list(matrix_TAMLR$LR),sum)
colnames(net_data)[2]<-'count'
net_data$ligand<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[1]))
net_data$receptor<-unlist(lapply(strsplit(net_data$Group.1,'_of_'),function(x)x[2]))

write.table(net_data,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_net_data_COLLAGEN.txt',quote = F,sep = '\t',row.names = F)

net_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_net_data_LAMININ.txt',stringsAsFactors = F,check.names = F)

##点数据
nodes <- data.frame(name = c(unique(net_data$ligand),unique(net_data$receptor)),
                    cluster=c(rep('ligand',length(unique(net_data$ligand))),rep('receptor',length(unique(net_data$receptor)))))

de<-rbind(as.data.frame(table(net_data$ligand)),
          as.data.frame(table(net_data$receptor)))
nodes$degree<-de$Freq[match(nodes$name,de$Var1)]

###边数据
edges <- net_data[c("ligand","receptor","count")]

g <- tbl_graph(nodes = nodes, edges = edges)

colors <- colorRampPalette(c("red", "orange", "blue"),space = "rgb")(3)

#?ggraph
p_net<-ggraph(g,layout='linear',circular = TRUE) +###'stress'  'linear'   'fabric'
  # geom_edge_bend(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  # geom_edge_arc(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  geom_edge_link(mapping = aes(edge_width = count),strength = 0.02,alpha = 0.5,color="#99B8D7") +
  #scale_edge_colour_manual(values = c("lightblue")) +
  scale_edge_width_continuous(range = c(0.5,1)) +
  geom_node_point(aes(size = degree,colour = cluster),alpha = 0.9) +
  scale_size_continuous(range = c(1,3)) +  #设置点大小范围，可以设置值越小，点越大
  scale_color_manual(values = c("#B992D3","#2C85BB")) +
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),angle=0,hjust=0,size=1.5) + # 设置点的注释
  theme_graph(title_size = 8,strip_text_size = 4,caption_size=3)+
  ggtitle('matrix_TAM_LAMININ')
p_net
svglite(paste0(dir_pic,'2/net_matrix_TAM_LAMININ.svg'),width = 3.5,height = 2.6)
print(p_net)
dev.off()

























