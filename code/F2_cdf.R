dir_out<-'/new_st_CytoTRACE/'
dir_copykat<-'/new_st_copykat/'
file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_copykat,recursive = T)
# file_bdy<-file_bdy[-grep('slice',file_bdy)]
length(file_bdy)

dir_rds<-'/new_st/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
length(file_rds)
library(stringr)
all_slice_stemScore <- c()
for(i in 1:length(file_bdy)){
  #i=1
  parts <- str_split(file_bdy[i],"/")[[1]]
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  stem_score<-read.delim(paste0(dir_out,parts[1],"/",parts[2],"/",str_split(parts[3],"_")[[1]][1],"_Stemness.txt"),stringsAsFactors = F,check.names = F)
  
  score<-data.frame(dataset_slice=paste0(parts[2],"/",str_split(parts[3],"_")[[1]][1]),
                        cell_name=st_bdy$cell_name,
                        LocalType=st_bdy$FinalLocalType,
                        CytoTRACE=stem_score$score,
                        cancer=parts[1])
  all_slice_stemScore<-rbind(all_slice_stemScore,score)
  
  print(parts)
}
write.table(all_slice_stemScore,paste0(dir_out,'all_slice_stemScore.txt'),quote = F,sep = '\t')


library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)

dir_dotplot <- "/luad_stemscore"
dir_copykat<-'/pan_cancer/copykat/'


###########slice_stem_score###########
all_slice_stem <- read.table("/all_slice_stemScore.txt",header = T,sep = "\t")
####luad
dir_select_slice <- "/dotplot"
select_slice <- list.files(path = dir_select_slice)
dataset_name<-gsub('_.*','',select_slice)
slice_name <- gsub(".*_","",select_slice)  
slice_name <- gsub("\\.pdf$", "", slice_name)  
dataset_slice_name <- paste(dataset_name,slice_name,sep = "/")
luad_slice <- paste(dataset_name,slice_name,sep = "_")
for (i in 1:length(dataset_slice_name)) {
  #i <- 1
  one_slice_stem <- all_slice_stem[all_slice_stem$dataset_slice==dataset_slice_name[i],]
  
  one_slice_st_bdy<-read.delim(paste0(dir_copykat,dataset_name[i],"\\",slice_name[i],"_BdyTumorCore.txt"),stringsAsFactors = F,check.names = F)
  one_slice_st_bdy<-one_slice_st_bdy[one_slice_stem$cell_name,]
  
  plot_data<-data.frame(imagerow=one_slice_st_bdy$imagerow,imagecol=one_slice_st_bdy$imagecol,
                        stemscore=one_slice_stem$CytoTRACE)
  p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
    geom_point(aes(colour=stemscore),size=.6) +
    scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(25),
                                      colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                      colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                      colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
    )+ #设置填充颜色
    theme_classic()+
    labs(title = dataset_slice_name[i],x = "",y = "")
  print(p)
  pdf(paste0(dir_dotplot,"\\",luad_slice[i],".pdf"),width = 4.8, height = 4)
  print(p)
  dev.off()
  
  
}

####10_slice
dir_select_slice <- "/plot_stemness/new_10slice"
select_slice <- list.files(path = dir_select_slice)

dir_dotplot <- "/plot_stemness/10slice_stemness_new"
dir_copykat<-'/new_st_copykat/' 
all_slice_stem <- read.table("/all_slice_stemScore.txt",header = T,sep = "\t")
for (i in 1:length(select_slice)) {
  #i <- 1
  dataset_name<-select_slice[i]
  cancer_name <- str_split(dataset_name,"_")[[1]][1]
  slice_name <- paste0(str_split(dataset_name,"_")[[1]][2],"/",str_split(dataset_name,"_")[[1]][3])
    
  one_slice_stem <- all_slice_stem[all_slice_stem$dataset_slice==slice_name,]
  
  one_slice_st_bdy<-read.delim(paste0(dir_copykat,cancer_name,"/",slice_name,"_BdyTumorCore.txt"),stringsAsFactors = F,check.names = F)
  one_slice_st_bdy<-one_slice_st_bdy[one_slice_stem$cell_name,]
  
  plot_data<-data.frame(imagerow=one_slice_st_bdy$imagerow,imagecol=one_slice_st_bdy$imagecol,
                        stemscore=one_slice_stem$CytoTRACE)
  p<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
    geom_point(aes(colour=stemscore),size=.6) +
    scale_color_gradientn(colours = c(colorRampPalette(c("#04040b","#5b2177"))(25),
                                      colorRampPalette(c("#5b2177","#b93c6d"))(25),
                                      colorRampPalette(c("#b93c6d","#eb7d60"))(25),
                                      colorRampPalette(c("#eb7d60","#f6f0b7"))(25))
    )+ #设置填充颜色
    theme_classic()+
    labs(title = dataset_name,x = "",y = "")
  #print(p)
  pdf(paste0(dir_dotplot,"/",dataset_name,".pdf"),width = 4.8, height = 4)
  print(p)
  dev.off()
  
  
}