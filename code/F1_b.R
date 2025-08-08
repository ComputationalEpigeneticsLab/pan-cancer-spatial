####地球图
library(plot1cell)
library(Seurat)
library(tidyverse)
library(stringr)
library(RColorBrewer)



integration_data<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/meta_order.txt',stringsAsFactors = F,check.names = F)
integration_data$LocationType[which(integration_data$LocationType=='Tumor')]<-'Dispersion'
integration_data$LocationType[which(integration_data$LocationType=='diploid')]<-'Normal'
#integration_data<-integration_data[integration_data$LocationType%in%c(c('Core','Boundary','Dispersion')),]

integration_umap<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/umap.txt',stringsAsFactors = F,check.names = F)
integration_umap<-integration_umap[rownames(integration_data),]
integration_umap$leiden<-integration_data$leiden
integration_umap$slice<-integration_data$slice
#integration_umap$cell_name<-integration_data$cell_name
integration_umap$LocationType<-integration_data$LocationType
#integration_umap<-integration_umap[c(1:1000,200000:210000,300000:350000),]

integration_umap$cells<-rownames(integration_umap)
colnames(integration_umap)[1:3]<-c('dim1','dim2','Cluster')
integration_umap$Cluster<-paste0('Cluster_',integration_umap$Cluster)
integration_umap$Cluster<-factor(integration_umap$Cluster,levels = paste0('Cluster_',0:12))
table(integration_umap$Cluster)
integration_umap$cancer<-unlist(lapply(strsplit(integration_umap$slice,'/'),function(x)x[1]))
integration_umap$cancer<-substr(integration_umap$cancer,1,nchar(integration_umap$cancer)-2) %>% toupper()
#integration_umap<-integration_umap[integration_umap$LocationType%in%c('Core','Boundary','Dispersion'),]
integration_umap<-integration_umap[which(integration_umap$LocationType!='not.defined'),]
integration_umap<-integration_umap[order(integration_umap$Cluster),]
table(integration_umap$Cluster)
integration_umap<-cell_order(integration_umap)####每种Cluster按其个数从1排序


integration_umap$x_polar2 <- log10(integration_umap$x_polar)
integration_umap<-integration_umap[,c(3,1,2,4:ncol(integration_umap))]
integration_umap$x<-transform_coordinates(integration_umap$dim1, zoom = 0.8)
integration_umap$y<-transform_coordinates(integration_umap$dim2, zoom = 0.8)

integration_umap<-mutate(integration_umap,Cluster = factor(integration_umap$Cluster, 
                                                           levels = paste0('Cluster_',0:12)))
class(integration_umap)
#circ_data<-circ_data[1:1000,]
color_celltype_big<-rand_color(length(levels(integration_umap$Cluster)))
color_celltype_big<-c("#B005C5", "#EFD347", "#0DC8B9", "#9EF06F", "#46A90C", "#6DEA4B", "#A0291F",
                      "#E9FEB0", "#FB2B92", "#F1B0A0", "#E704C4", "#DF83E6", "#91060D")

color_loca2<-c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
               'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
               'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
               'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b')
cluster_color<-c("Boundary"="#F6B86D","Core"="#D62D28","Dispersion"="#EE762D",'Immune'='#a9d38a','Normal'='#2375ae')

#png(filename =  'E:/Mirror/ST_analysis/pic/re/circlize_plot.png', width = 6, height = 6,units = 'in', res = 300)
pdf('E:/Mirror/ST_analysis/pic/re/circlize_plot3.pdf', width = 8, height = 8)
plot_circlize(integration_umap,do.label = F, pt.size = 0.02,
              col.use = color_celltype_big ,
              bg.color = 'white',#'#f8f2e4',
              kde2d.n = 200, repel = F, label.cex = 0.7)

add_track(integration_umap, group = "cancer", colors = color_loca2, track_num = 2) 
add_track(integration_umap, group = "LocationType",colors = cluster_color, track_num = 3)
dev.off()


?plot_circlize


transform_coordinates <- function(
    coord_data, 
    zoom
){
  center_data<-coord_data-mean(c(min(coord_data),max(coord_data)))
  max_data<-max(center_data)
  new_data<-center_data*zoom/max_data
  new_data
}

data_plot<-integration_umap
kde2d.n = 1000

data_plot %>%
  dplyr::group_by(Cluster) %>%
  summarise(x = median(x = x), y = median(x = y)) -> centers
z <- MASS::kde2d(data_plot$x, data_plot$y, n=kde2d.n)
celltypes<-names(table(data_plot$Cluster))
cell_colors <- scales::hue_pal()(length(celltypes))
if(!is.null(col.use)){
  cell_colors=col.use
  col_df<-data.frame(Cluster=celltypes, color2=col.use)
  cells_order<-rownames(data_plot)
  data_plot<-merge(data_plot, col_df, by="Cluster")
  rownames(data_plot)<-data_plot$cells
  data_plot<-data_plot[cells_order,]
  data_plot$Colors<-data_plot$color2
}
circos.clear()
bg.color='#F9F2E4'
par(bg = bg.color)
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA,panel.fun = function(x, y) {
  # circos.text(CELL_META$xcenter,
  #             CELL_META$cell.ylim[2]+ mm_y(4),
  #             CELL_META$sector.index,
  #             cex=0.5, col = 'black', facing = "bending.inside", niceFacing = T)
  #circos.axis(labels.cex = 0.3, col = 'black', labels.col =  'black')
})
for(i in 1:length(celltypes)){
  dd<-data_plot[data_plot$Cluster==celltypes[i],]
  circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = cell_colors[i],  lwd=3, sector.index = celltypes[i])
}
text(x = 1, y=0.1, labels = "Cluster", cex = 0.4, col = 'black',srt=-90)
pt.size = 0.5
points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.2), cex = pt.size);
contour.nlevels = 100
contour.levels = c(0.2,0.3)
contour(z, drawlabels=F, nlevels= 100, levels = contour.levels,col = '#ae9c76', add=TRUE)
do.label = T
repel=FALSE
label.cex = 0.5
if(do.label){
  if(repel){
    textplot(x=centers$x, y=centers$y, words =  centers$Cluster,cex = label.cex, new = F,show.lines=F)
  } else {
    text(centers$x,centers$y, labels=centers$Cluster, cex = label.cex, col = 'black')
  }
}

add_track(data_plot, group = "cancer", colors = color_loca2, track_num = 2) 
add_track(data_plot, group = "LocationType",colors = cluster_color, track_num = 3)





plot_circlize <- function(
    data_plot,
    do.label = T,
    contour.levels = c(0.2,0.3),
    pt.size = 0.5,
    kde2d.n = 1000,
    contour.nlevels = 100,
    bg.color='#F9F2E4',
    col.use=NULL,
    label.cex = 0.5,
    repel=FALSE
) {
  data_plot %>%
    dplyr::group_by(Cluster) %>%
    summarise(x = median(x = x), y = median(x = y)) -> centers
  z <- MASS::kde2d(data_plot$x, data_plot$y, n=kde2d.n)
  celltypes<-names(table(data_plot$Cluster))
  cell_colors <- scales::hue_pal()(length(celltypes))
  if(!is.null(col.use)){
    cell_colors=col.use
    col_df<-data.frame(Cluster=celltypes, color2=col.use)
    cells_order<-rownames(data_plot)
    data_plot<-merge(data_plot, col_df, by="Cluster")
    rownames(data_plot)<-data_plot$cells
    data_plot<-data_plot[cells_order,]
    data_plot$Colors<-data_plot$color2
  }
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12),points.overflow.warning=FALSE)
  circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA,panel.fun = function(x, y) {
    # circos.text(CELL_META$xcenter,
    #             CELL_META$cell.ylim[2]+ mm_y(4),
    #             CELL_META$sector.index,
    #             cex=0.5, col = 'black', facing = "bending.inside", niceFacing = T)
    #circos.axis(labels.cex = 0.3, col = 'black', labels.col =  'black')
  })
  for(i in 1:length(celltypes)){
    dd<-data_plot[data_plot$Cluster==celltypes[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = cell_colors[i],  lwd=3, sector.index = celltypes[i])
  }
  text(x = 1, y=0.1, labels = "Cluster", cex = 0.4, col = 'black',srt=-90)
  points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.2), cex = pt.size);
  contour(z, drawlabels=F, nlevels= 100, levels = contour.levels,col = '#ae9c76', add=TRUE)
  if(do.label){
    if(repel){
      textplot(x=centers$x, y=centers$y, words =  centers$Cluster,cex = label.cex, new = F,show.lines=F)
    } else {
      text(centers$x,centers$y, labels=centers$Cluster, cex = label.cex, col = 'black')
    }
  } 
}




meta_order<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/meta_order.txt',stringsAsFactors = F,check.names = F)
meta_umap<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/intergration/umap.txt',stringsAsFactors = F,check.names = F,row.names = 1)
meta_umap$data_slice<-meta_order[rownames(meta_umap),'slice']
meta_umap$cancer<-unlist(lapply(strsplit(meta_umap$data_slice,'/'),function(x) x[1]))
meta_umap$cancer<-toupper(substr(meta_umap$cancer,1,nchar(meta_umap$cancer)-2))
table(meta_umap$cancer)

meta_umap$cluster<-meta_order[rownames(meta_umap),'leiden']
colnames(meta_umap)[1:2]<-c('UMAP_1','UMAP_2')

p1<-ggplot(data = meta_umap,aes(x = UMAP_1, y = UMAP_2,color=cancer)) + 
  geom_point(size=0.1)+
  scale_color_manual(values = c('BRCA'='#b7996d','CESC'='#e32427','CRC'='#8bc96d','CSCC'='#b05a28','GBM'='#a4cde1',
                                'GIST'='#96cb8f','HGSC'='#277fb8','HN-AS'='#f38989','IPMN'='#5c9e43','LIHC'='#c6b598',
                                'LUAD'='#7A9AC5','MIBC'='#60592E','OSCC'='#C5BE97','OVCA'='#C89192','PCNSL'='#44637F',
                                'PDAC'='#549da3','PRAD'='#f9b769','RCC'='#af93c4','SKCM'='#d4a55b'))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle('')
png('E:/Mirror/ST_analysis/pic/re/1/integration_cancer.png',width = 800,height = 600)
print(p1)
dev.off()
pdf('E:/Mirror/ST_analysis/pic/re/1/integration_cancer.pdf',width = 8,height = 6)
print(p1)
dev.off()
#?png

meta_umap$cluster<-paste0('Cluster',meta_umap$cluster)
unique(meta_umap$cluster)
p1<-ggplot(data = meta_umap,aes(x = UMAP_1, y = UMAP_2,color=cluster)) + 
  geom_point(size=0.1)+
  scale_color_manual(values = c('Cluster0'='#D8D9C1','Cluster1'='#D7B7DF','Cluster2'='#ECBEAA','Cluster3'='#5F91C7',
                                'Cluster4'='#DC98BF','Cluster5'='#E4CFB2','Cluster6'='#EFA7A9','Cluster7'='#7CBFB6',
                                'Cluster8'='#C2E3EC','Cluster9'='#B196C1','Cluster10'='#F9D977','Cluster11'='#A7CBEF',
                                'Cluster12'='#F6D6DA'))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle('')
png('E:/Mirror/ST_analysis/pic/re/1/integration_cluster.png',width = 800,height = 600)
print(p1)
dev.off()



library(entropy)
###熵值计算
###1、对每个簇的样本来源算熵值
sample_entropy<-lapply(unique(meta_umap$cluster),function(x){#x=unique(meta_umap$cluster)[1]
  aa<-meta_umap[which(meta_umap$cluster==x),]
  bb<-as.data.frame(table(aa$data_slice))
  bb$Freq<-bb$Freq/sum(bb$Freq)
  entropy(bb$Freq, unit = 'log2')
}) %>% unlist()
sample_entropy<-data.frame(cluster=unique(meta_umap$cluster),sample_entropy=sample_entropy)

meta_umap_en<-merge(meta_umap,sample_entropy,by='cluster')
p1<-ggplot(data = meta_umap_en,aes(x = UMAP_1, y = UMAP_2,color=sample_entropy)) + 
  geom_point(size=0.1)+
  scale_color_gradientn(colours = c(colorRampPalette(c("#E1CFC4","#F39C67"))(50),###数字可以改
                                    colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle('')
png('E:/Mirror/ST_analysis/pic/re/1/integration_sample_entropy.png',width = 800,height = 600)
print(p1)
dev.off()
pdf('E:/Mirror/ST_analysis/pic/re/1/integration_sample_entropy.pdf',width = 8,height = 6)
print(p1)
dev.off()


###2、对每个簇的癌症来源算熵值
cancer_entropy<-lapply(unique(meta_umap$cluster),function(x){#x=unique(meta_umap$cluster)[1]
  aa<-meta_umap[which(meta_umap$cluster==x),]
  bb<-as.data.frame(table(aa$cancer))
  bb$Freq<-bb$Freq/sum(bb$Freq)
  entropy(bb$Freq, unit = 'log2')
}) %>% unlist()
cancer_entropy<-data.frame(cluster=unique(meta_umap$cluster),cancer_entropy=cancer_entropy)

meta_umap_en<-merge(meta_umap,cancer_entropy,by='cluster')
p1<-ggplot(data = meta_umap_en,aes(x = UMAP_1, y = UMAP_2,color=cancer_entropy)) + 
  geom_point(size=0.1)+
  scale_color_gradientn(colours = c(colorRampPalette(c("#E1CFC4","#F39C67"))(50),###数字可以改
                                    colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=2)))+
  ggtitle('')
png('E:/Mirror/ST_analysis/pic/re/1/integration_cancer_entropy.png',width = 800,height = 600)
print(p1)
dev.off()
pdf('E:/Mirror/ST_analysis/pic/re/1/integration_cancer_entropy.pdf',width = 8,height = 6)
print(p1)
dev.off()


meta_umap$copykat<-meta_order$copykat.pred
meta_umap_p<-meta_umap[which(meta_umap$copykat!='not.defined'),]
table(meta_umap_p$copykat)
p1<-ggplot(data = meta_umap_p,aes(x = UMAP_1, y = UMAP_2,color=copykat)) + 
  geom_point(size=0.1)+
  scale_color_manual(values = c("aneuploid"="#CD3333","diploid"="#37A0A0",
                                "not.defined"="#CCCCCC"####在这里not.defined的点需要提前先删了
  ))+
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  ggtitle('')
png('E:/Mirror/ST_analysis/pic/re/1/integration_copykat.png',width = 800,height = 600)
print(p1)
dev.off()
pdf('E:/Mirror/ST_analysis/pic/re/1/integration_copykat.pdf',width = 8,height = 6)
print(p1)
dev.off()





