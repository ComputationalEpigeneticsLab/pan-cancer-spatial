####spot标签改名字
library(tidyverse)
library(ggplot2)
library(ggpubr)

dir_bdy<-"E:/Mirror/ST_analysis/data/10X Visium/copykat/"
file_bdy<-list.files(pattern = "_BdyTumorCore.txt",path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_bdy,"_"),function(x) x[1]))

pdf("E:/Mirror/ST_analysis/pic/copykat/TumorCoreBdy_4.pdf",width = 5, height = 4)
for(i in 1:length(file_bdy)){
  #i=168
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  #st_bdy<-st_nearSpot
  # table(st_bdy$LocationType)
  # st_bdy$FinalLocalType<-st_bdy$LocationType
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="Tumor")]<-"Dispersion"
  # st_bdy$FinalLocalType[which(st_bdy$LocationType=="diploid")]<-"Normal"
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
  table(st_bdy$FinalLocalType)
  #write.table(st_bdy,paste0(dir_bdy,file_bdy[i]),quote = F,sep = '\t')
  
  p_spot2<-ggplot(data = st_bdy,aes(x = imagerow, y = imagecol,color=FinalLocalType)) + 
    geom_point(size=1)+
    # scale_color_manual(values =  c("Dispersion"="#CD5A5A","Boundary"="#FF9900","Normal"="#669933","Core"="#990033",
    #                                "Immune"="#5477AF","not.defined"="#CCCCCC"))+
    scale_color_manual(values =  c('Core'='#d62d28','Boundary'='#f6b86d','Dispersion'='#ee762d',
                                   'Immune'='#a9d38a','Normal'='#2375ae'))+
    theme_classic()+
    ggtitle(dataset_slice[i])+
    guides(colour = guide_legend(override.aes = list(size=2)))
  print(p_spot2)
  # pdf('E:/Mirror/ST_analysis/pic/re/4/bdy_core_dis.pdf',width = 4.6, height = 4)
  # print(p_spot2)
  # dev.off()
  
}
dev.off()



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



