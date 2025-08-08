####干性得分在spot中的得分
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)


dir_pic<-"E:/Mirror/ST_analysis/pic/MP_score/"
dir_pic<-"E:/Mirror/ST_analysis/pic/MP_score_2/"
dir_copykat<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_stem<-list.files(pattern = 'Stemness.txt',path = dir_copykat,recursive = T)
data_slice<-unlist(lapply(strsplit(file_stem,'_Stemness'),function(x) x[1]))
file_bdy<-list.files(pattern = "BdyTumorCore.txt",path = dir_copykat,recursive = T)
table(data_slice==unlist(lapply(strsplit(file_bdy,'_Bdy'),function(x) x[1])))



pdf(paste0(dir_pic,'st_stem.pdf'),width = 4.8, height = 4)
for(i in 1:length(file_MPscore)){
  #i=1
  st_stem<-read.delim(paste0(dir_copykat,file_stem[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[rownames(st_stem),]
  
  plot_data<-data.frame(imagerow=st_bdy$imagerow,imagecol=st_bdy$imagecol,
                        stem=st_stem$CytoTRACE)
  # plot_data$stem[which(st_bdy$FinalLocalType!='Core')]<--0.001
  p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
    geom_point(aes(colour=stem),size=1) +
    # scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
    #                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
    scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(1),
                                      colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                      colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                      colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
    )+
    theme_classic()+
    labs(title = data_slice[i],x = "",y = "")
  print(p)
}
dev.off()
