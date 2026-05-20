####MP14在spot中的得分
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)



dir_pic<-"/data/10X_Visium/plot/MPscore/F3d/"

dir_MPscore<-'/data/10X_Visium/new_st_NMF/MPscore_new/'
file_MPscore<-list.files(pattern = 'MP_score.txt',path = dir_MPscore,recursive = T)
data_slice<-unlist(lapply(strsplit(file_MPscore,'_MP'),function(x) x[1]))

dir_copykat<-'/data/10X_Visium/new_st_copykat/'
file_bdy<-list.files(pattern = "BdyTumorCore.txt",path = dir_copykat,recursive = T)
table(data_slice==unlist(lapply(strsplit(file_bdy,'_Bdy'),function(x) x[1])))

MP<-paste0('MP_',1:17)


for(j in 1:length(MP)){
  #j=1
  pdf(paste0(dir_pic,MP[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(file_MPscore)){
    #i=1
    st_MPscore<-read.delim(paste0(dir_MPscore,file_MPscore[i]),stringsAsFactors = F,check.names = F)
    st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
    st_bdy<-st_bdy[rownames(st_MPscore),]
    
    plot_data<-data.frame(imagerow=st_bdy$imagerow,imagecol=st_bdy$imagecol,
                          MPscore=st_MPscore[,MP[j]])
    plot_data$MPscore[which(st_bdy$FinalLocalType!='Core')]<--0.001
    p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
      geom_point(aes(colour=MPscore),size=1) +
      # scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
      #                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
      scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(1),
                                        colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                        colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                        colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
      )+
      theme_classic()+
      labs(title = data_slice[i],x = "",y = "",col=MP[j])
    print(p)
  }
  dev.off()
}


