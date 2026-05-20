####MP14在spot中的得分
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)



dir_pic<-"/MP/"

dir_MPscore<-'/MPscore/MPscore/ESCC/'
file_MPscore<-list.files(pattern = 'MP_score.txt',path = dir_MPscore,recursive = T)
data_slice<-unlist(lapply(strsplit(file_MPscore,'_MP'),function(x) x[1]))

dir_copykat<-'/new_st_copykat/ESCC/'
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

##eco
eco<-paste0('ecosystem',1:6)

for(j in 1:6){
  #j=1
  pdf(paste0(dir_pic,eco[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(file_MPscore)){
    #i=1
    st_MPscore<-read.delim(paste0(dir_MPscore,file_MPscore[i]),stringsAsFactors = F,check.names = F)
    st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
    st_bdy<-st_bdy[rownames(st_MPscore),]
      
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_8')]<-'ecosystem1'
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_7','MP_17')]<-'ecosystem3'
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_1','MP_4')]<-'ecosystem4'
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
    st_MPscore$ecosystem[st_MPscore$MP_top1%in%c('MP_15','MP_11')]<-'ecosystem6'
    
    st_MPscore$ecosystem1_score <- st_MPscore$MP_8
    st_MPscore$ecosystem2_score <- rowSums(st_MPscore[, c('MP_6', 'MP_3', 'MP_12', 'MP_13', 'MP_14', 'MP_16')], na.rm = TRUE)
    st_MPscore$ecosystem3_score <- rowSums(st_MPscore[, c('MP_7', 'MP_17')], na.rm = TRUE)
    st_MPscore$ecosystem4_score <- rowSums(st_MPscore[, c('MP_1', 'MP_4')], na.rm = TRUE)
    st_MPscore$ecosystem5_score <- rowSums(st_MPscore[, c('MP_9', 'MP_10', 'MP_2', 'MP_5')], na.rm = TRUE)
    st_MPscore$ecosystem6_score <- rowSums(st_MPscore[, c('MP_15', 'MP_11')], na.rm = TRUE)
      
    ecosystem_scores <- st_MPscore[, c('ecosystem1_score', 'ecosystem2_score', 
                                       'ecosystem3_score', 'ecosystem4_score', 
                                       'ecosystem5_score', 'ecosystem6_score')]
    # 行归一化
    ecosystem_norm <- t(apply(ecosystem_scores, 1, function(x) x / sum(x, na.rm = TRUE)))
    
    # 将归一化后的结果添加到数据框中
    colnames(ecosystem_norm) <- paste0(colnames(ecosystem_scores), "_norm")
    st_MPscore <- cbind(st_MPscore, ecosystem_norm)
    
    plot_data<-data.frame(imagerow=st_bdy$imagerow,imagecol=st_bdy$imagecol,
                          ecoscore=st_MPscore[,26+j])
    plot_data$ecoscore[which(st_bdy$FinalLocalType!='Core')]<--0.001
    p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
      geom_point(aes(colour=ecoscore),size=1) +
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

for(j in 1:6){
  #i=1
  st_MPscore<-read.delim(paste0(dir_MPscore,file_MPscore[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[rownames(st_MPscore),]
  st_MPscore$ecosystem<-''
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_8')]<-'ecosystem1'
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_7','MP_17')]<-'ecosystem3'
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_1','MP_4')]<-'ecosystem4'
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
  st_MPscore$ecosystem[st_MPscore$top_MP%in%c('MP_15','MP_11')]<-'ecosystem6'
  
  plot_data<-data.frame(imagerow=st_bdy$imagerow,imagecol=st_bdy$imagecol,
                        MPscore=st_MPscore[,eco[j]])
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




