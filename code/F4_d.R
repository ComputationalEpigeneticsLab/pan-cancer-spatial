####不同MP类别的core的bdy邻居的Imm或Normal邻居的免疫类型组成
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)
####绘制成环形热图
library(OmicCircos)
library(S4Vectors)
library(tidyverse)

library(tidyverse)
library(magrittr)
library(circlize)
library(ComplexHeatmap)
library(dendextend)
library(reshape2)


#####################################################################################################
###bdy_step1
dir_re<-'E:/Mirror/ST_analysis/data/pic_data/bdy_step1_MP/bdy_step1_core_bdy_step1_mp_ImmuneInfiltration_immune/'
file_re<-list.files(pattern = 'txt',path = dir_re,recursive = T)

re_slice<-unlist(lapply(strsplit(file_re,'_M'),function(x)x[1]))
re_MP<-paste0('M',unlist(lapply(strsplit(file_re,'_M'),function(x)x[2])))
re_MP<-gsub('.txt','',re_MP)

allSliceMP<-c()
for(i in 1:length(file_re)){
  #i=1
  re_data<-read.delim(paste0(dir_re,file_re[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  apply(re_data,2,mean)
  allSliceMP<-rbind(allSliceMP,
                    data.frame(celltype=colnames(re_data),
                               value=apply(re_data,2,mean),
                               slice=re_slice[i],
                               MP=re_MP[i]))
}
allSliceMP<-allSliceMP[which(allSliceMP$value%in%NaN==F),]
allSliceMP$cancer<-unlist(lapply(strsplit(allSliceMP$slice,'_'),function(x)x[1]))
allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()


dir_re<-'E:/Mirror/ST_analysis/data/pic_data/bdy_MP/bdy_mp_step1_ImmuneInfiltration_immune_normal/'
file_re<-list.files(pattern = 'txt',path = dir_re,recursive = T)

re_slice<-unlist(lapply(strsplit(file_re,'_M'),function(x)x[1]))
re_MP<-paste0('M',unlist(lapply(strsplit(file_re,'_M'),function(x)x[2])))
re_MP<-gsub('.txt','',re_MP)

allSliceMP<-c()
for(i in 1:length(file_re)){
  #i=1
  re_data<-read.delim(paste0(dir_re,file_re[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  #apply(re_data,2,mean)
  allSliceMP<-rbind(allSliceMP,
                    data.frame(celltype=colnames(re_data),
                               value=apply(re_data,2,mean),
                               slice=re_slice[i],
                               MP=re_MP[i],
                               s_MP=file_re[i]))
}
allSliceMP<-allSliceMP[which(allSliceMP$value%in%NaN==F),]
allSliceMP$cancer<-unlist(lapply(strsplit(allSliceMP$slice,'_'),function(x)x[1]))
allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()


allSliceMP$slice_MP<-paste0(allSliceMP$slice,'of',allSliceMP$MP)

unique(file_re) %>% length()
length(unique(allSliceMP$slice_MP))
cir_data<-as.data.frame(acast(allSliceMP[,c('slice_MP','celltype','value')],slice_MP~celltype))

cir_data<-cir_data[,c('CAF','TAM',setdiff(colnames(cir_data),c('CAF','TAM')))]
cir_data$MP<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[2]))
cir_data$slice<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[1]))
cir_data$cancer<-unlist(lapply(strsplit(cir_data$slice,'_'),function(x)x[1]))
cir_data$cancer<-substr(cir_data$cancer,1,nchar(cir_data$cancer)-2) %>% toupper()


cir_data$order<-unlist(lapply(strsplit(cir_data$MP,'_'),function(x) x[2]))%>%as.numeric()
cir_data<-cir_data[order(cir_data$cancer),]
cir_data<-cir_data[order(cir_data$order),]



range(cir_data[,1:14],na.rm = T)
plot_col=colorRamp2(seq(0,0.4,length.out=5),c("#463380",'#365d8f',"#6dc060",'#85c551',"#bfd631"))##"#133868",'#238DC1',
plot_col=colorRamp2(seq(0,1,length.out=3),c("white",'#D45848',"#6C0E24"))##"#133868",'#238DC1',

MP_col=c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
         "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
         "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
         'MP_14'='#FF9694')
names(MP_col)<-unique(cir_data$MP)

lev_split2<-factor(cir_data$order)
#pdf("E:/Mirror/ST_analysis/pic/re/1/circos/circos_cluster.pdf",width = 10,height = 10)
pdf("E:/Mirror/ST_analysis/pic/NMF_module/coreBdyImm/cir_celltype_MP.pdf",width = 10,height = 10)
circos.clear()
#circos.par(gap.degree=1)
#circos.par(gap.after = c(rep(1,12), 10))
circos.par(gap.after = c(rep(1,13)))
#?circos.par
circos.heatmap(cir_data['MP'], col = MP_col, na.col = "#D9D9D9",
               track.height = 0.02,###热图条带的宽度
               track.margin=c(0.02,0),###前面的值是与向内的下一个热图轨迹的间隔，后面的值是与向外的上一个热图轨迹的间隔
               #rownames.side = "outside",
               cluster = F,
               cell.border=NA,
               cell.lwd=0.02,
               split = lev_split2 , 
               rownames.cex = 0.6)
clusterleg=Legend(title = "MP", at = names(MP_col), 
                  legend_gp = gpar(fill = MP_col ))
grid.draw(clusterleg)
dev.off()





