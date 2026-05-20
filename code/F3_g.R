####MP 分5类的各个免疫通路得分绘制
library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)


dir_pic<-'/10X_Visium/plot/'
dir_18Imm<-'/10X_Visium/new_st_18Imm/'
file_metas<-list.files(pattern = '_18Imm.txt',path = dir_18Imm,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_metas,'_'),function(x) x[1]))

#pathway<-colnames(st_MPscore)[1:18]

pathway<-c("Cytolytic_markers","HLA_molecules","IFN_y_pathway_genes","chemokines","adhesion_molecules","Tcell_inflamed",
           "APC_co_stimulation","APC_co_inhibition","CCR","Check_point","HLA","Inflammation_promoting",
           "MHC_class_I","Parainflammation","T_cell_co_inhibition","T_cell_co_stimulation","Type_I_IFN_Reponse","Type_II_IFN_Reponse")

dir_bdy<-'/10X_Visium/new_st_copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)

for(j in 1:length(pathway)){
  #j=1
  pdf(paste0(dir_pic,pathway[j],'.pdf'),width = 4.8, height = 4)
  for(i in 1:length(dataset_slice)){
    #i=1
    st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
    
    st_MPscore<-read.delim(paste0(dir_18Imm,file_metas[i]),stringsAsFactors = F,check.names = F)

    st_MPscore$LocalType<-st_bdy$FinalLocalType
    st_MPscore$imagerow<-st_bdy$imagerow
    st_MPscore$imagecol<-st_bdy$imagecol
    
    plot_data<-data.frame(imagerow=st_MPscore$imagerow,imagecol=st_MPscore$imagecol,
                          MPscore=st_MPscore[,pathway[j]])
    plot_data$MPscore<-plot_data$MPscore-min(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    plot_data$MPscore<-plot_data$MPscore/max(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    range(plot_data$MPscore[which(st_MPscore$LocalType=='Core')])
    plot_data$MPscore[which(st_MPscore$LocalType!='Core')]<--0.001
    p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
      geom_point(aes(colour=MPscore),size=.6) +
      # scale_color_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
      #                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
      scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(1),
                                        colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                        colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                        colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
      )+
      theme_classic()+
      labs(title = paste0(dataset_slice[i],'_of_',pathway[j]),x = "",y = "",col='score')
    print(p)
  }
  dev.off()
}
















