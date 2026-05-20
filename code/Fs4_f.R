library(Seurat)
library(dplyr)
library(jsonlite)
library(stringr)
library(ggplot2)
library(ggpubr)

###########slice_gene_exp###########
RBP_gene <- c("HNRNPF","HNRNPH1","ELAVL1")
TF_gene <- c("TRIM28","RXRB","MLX","CREB1","CTCF")

dir_rds<-'/10X_Visium/new_st/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)

dir_copykat<-'/10X_Visium/new_st_copykat/'
file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_copykat,recursive = T)

dir_exp <- "/10X_Visium/plot/TFRBP_exp/"


dir_select_slice <- "/10X_Visium/plot/10slice_stemness"
select_slice <- list.files(path = dir_select_slice)
slice_name <- gsub("\\.pdf$", "", select_slice)  

######RBP#######
for(i in 1:length(slice_name)){
  #i=1
  dataset_name<-slice_name[i]
  dataset_name<-str_split(dataset_name,"_")[[1]]
  selset_dataset_name <- paste0(dataset_name[1],"/",dataset_name[2],"/",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,dataset_name[1],"/",dataset_name[2],"/processed_",dataset_name[3],".rds"))
  st_bdy<-read.delim(paste0(dir_copykat,selset_dataset_name,"_BdyCoreBud.txt"),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[RBP_gene,]
  
  RBP_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                             LocalType=st_bdy$FinalLocalType,
                             dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                             imagerow=st_bdy$imagerow,
                             imagecol=st_bdy$imagecol)
  RBP_gene_exp <- cbind(RBP_gene_exp,t(st_conut))
  write.table(RBP_gene_exp,paste0(dir_exp,paste0(dataset_name[2],"_",dataset_name[3]),'_RBP_gene.txt'),quote = F,sep = '\t')
  print(selset_dataset_name)
}


#####TF_gene######
for(i in 1:length(slice_name)){
  #i=1
  dataset_name<-slice_name[i]
  dataset_name<-str_split(dataset_name,"_")[[1]]
  selset_dataset_name <- paste0(dataset_name[1],"/",dataset_name[2],"/",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,dataset_name[1],"/",dataset_name[2],"/processed_",dataset_name[3],".rds"))
  st_bdy<-read.delim(paste0(dir_copykat,selset_dataset_name,"_BdyCoreBud.txt"),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[TF_gene,]
  
  TF_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                            LocalType=st_bdy$FinalLocalType,
                            dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                            imagerow=st_bdy$imagerow,
                            imagecol=st_bdy$imagecol)
  TF_gene_exp <- cbind(TF_gene_exp,t(st_conut))
  write.table(TF_gene_exp,paste0(dir_exp,paste0(dataset_name[2],"_",dataset_name[3]),'_TF_gene.txt'),quote = F,sep = '\t')
  print(selset_dataset_name)
}

######plot####
dir_TF_plot <- "/10X_Visium/plot/TFRBP_plot/TF"
dir_exp <- "/10X_Visium/plot/TFRBP_exp/TF/"
for(i in 1:length(slice_name)){
  dataset_name<-slice_name[i]
  dataset_name<-str_split(dataset_name,"_")[[1]]
  selset_dataset_name <- paste0(dataset_name[1],"/",dataset_name[2],"/",dataset_name[3])
  
  TF_gene_exp <- read.delim(paste0(dir_exp,dataset_name[2],"_",dataset_name[3],"_TF_gene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  TF_gene_exp <- removeColsAllNa(TF_gene_exp)
  if(dim(TF_gene_exp)[2]>6){
    gene_data <- TF_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    
    for(j in 1:length(gene_name)){
      gene_name0 <- gene_name[j]
      plot_data<-data.frame(imagerow=TF_gene_exp$imagerow,imagecol=TF_gene_exp$imagecol,TF_gene_exp=TF_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","TF_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=TF_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = paste0(dataset_name[2],"_",dataset_name[3]),x = "",y = "")
      
      pdf(paste0(dir_TF_plot,"/",dataset_name[2],"_",dataset_name[3],"_TF_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(paste0(dataset_name[2],"_",dataset_name[3]))
}


dir_RBP_plot <- "/10X_Visium/plot/TFRBP_plot/RBP"
dir_exp <- "/10X_Visium/plot/TFRBP_exp/RBP/"
for(i in 1:length(slice_name)){
  dataset_name<-slice_name[i]
  dataset_name<-str_split(dataset_name,"_")[[1]]
  selset_dataset_name <- paste0(dataset_name[1],"/",dataset_name[2],"/",dataset_name[3])
  
  RBP_gene_exp <- read.delim(paste0(dir_exp,dataset_name[2],"_",dataset_name[3],"_RBP_gene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  RBP_gene_exp <- removeColsAllNa(RBP_gene_exp)
  if(dim(RBP_gene_exp)[2]>6){
    gene_data <- RBP_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    for (g in 1:length(gene_name)) {
      #g <- 1
      gene_name0 <- gene_name[g]
      plot_data<-data.frame(imagerow=RBP_gene_exp$imagerow,imagecol=RBP_gene_exp$imagecol,RBP_gene_exp=RBP_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","RBP_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=RBP_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = paste0(dataset_name[2],"_",dataset_name[3]),x = "",y = "")
      
      pdf(paste0(dir_RBP_plot,"/",paste0(dataset_name[2],"_",dataset_name[3]),"_RBP_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(paste0(dataset_name[2],"_",dataset_name[3]))
}  

