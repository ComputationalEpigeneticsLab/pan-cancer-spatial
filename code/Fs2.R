###########slice_gene_exp###########
dir_rds<-'/new_st/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))

dir_copykat<-'/new_st_copykat/'
file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_copykat,recursive = T)

dir_exp <- "/geneexp/"
######core#######
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[core_gene,]
  
  core_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                             LocalType=st_bdy$FinalLocalType,
                             dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                             imagerow=st_bdy$imagerow,
                             imagecol=st_bdy$imagecol)
  core_gene_exp <- cbind(core_gene_exp,t(st_conut))
  write.table(core_gene_exp,paste0(dir_exp,slice_name,'_coregene.txt'),quote = F,sep = '\t')
  print(slice_name)
}
#####immune######
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[immune_gene,]
  
  immune_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                               LocalType=st_bdy$FinalLocalType,
                               dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                               imagerow=st_bdy$imagerow,
                               imagecol=st_bdy$imagecol)
  immune_gene_exp <- cbind(immune_gene_exp,t(st_conut))
  write.table(immune_gene_exp,paste0(dir_exp,slice_name,'_immunegene.txt'),quote = F,sep = '\t')
  print(slice_name)
}

###########boundary_gene/dispersion_gene############
#############boundary#############
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[boundary_gene,]
  
  boundary_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                                 LocalType=st_bdy$FinalLocalType,
                                 dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                                 imagerow=st_bdy$imagerow,
                                 imagecol=st_bdy$imagecol)
  boundary_gene_exp <- cbind(boundary_gene_exp,t(st_conut))
  write.table(boundary_gene_exp,paste0(dir_exp,slice_name,'_boundarygene.txt'),quote = F,sep = '\t')
  print(slice_name)
}

#############budding#############
for(i in 1:length(dataset_slice)){
  #i=1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])
  
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_copykat,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[colnames(st_rds),]
  
  
  st_conut<-st_rds@assays[["Spatial"]]@counts%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[dispersion_gene,]
  
  dispersion_gene_exp <-data.frame(cell_name=st_bdy$cell_name,
                                   LocalType=st_bdy$FinalLocalType,
                                   dataSlice=paste0(dataset_name[2],"_",dataset_name[3]),
                                   imagerow=st_bdy$imagerow,
                                   imagecol=st_bdy$imagecol)
  dispersion_gene_exp <- cbind(dispersion_gene_exp,t(st_conut))
  write.table(dispersion_gene_exp,paste0(dir_exp,slice_name,'_dispersiongene.txt'),quote = F,sep = '\t')
  print(slice_name)
}

######plot####
dir_core_plot <- "/plot_core"
for(i in 1:length(dataset_slice)){
  #i <- 760
  #i <- 1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])

  core_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_coregene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  core_gene_exp <- removeColsAllNa(core_gene_exp)
  if(dim(core_gene_exp)[2]>6){
    gene_data <- core_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    
    for(j in 1:length(gene_name)){
      gene_name0 <- gene_name[j]
      plot_data<-data.frame(imagerow=core_gene_exp$imagerow,imagecol=core_gene_exp$imagecol,core_gene_exp=core_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","core_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=core_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = dataset_he_slice[i],x = "",y = "")
      
      pdf(paste0(dir_core_plot,"/",slice_name,"_core_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(slice_name)
}

#######immune######## 
dir_immune_plot <- "/plot_immune"
for(i in 1:length(dataset_slice)){
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])

  immune_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_immunegene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  immune_gene_exp <- removeColsAllNa(immune_gene_exp)
  if(dim(immune_gene_exp)[2]>6){
    gene_data <- immune_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    for (g in 1:length(gene_name)) {
      #g <- 1
      gene_name0 <- gene_name[g]
      plot_data<-data.frame(imagerow=immune_gene_exp$imagerow,imagecol=immune_gene_exp$imagecol,immune_gene_exp=immune_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","immune_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=immune_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = slice_name,x = "",y = "")
      
      pdf(paste0(dir_immune_plot,"/",slice_name,"_immune_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  print(slice_name)
}  

########boundary##########
dir_boundary_plot <- "/plot_boundary"
for(i in 1:length(dataset_slice)){
  #i <- 1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])

  boundary_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_boundarygene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  boundary_gene_exp <- removeColsAllNa(boundary_gene_exp)
  if(dim(boundary_gene_exp)[2]>6){
    gene_data <- boundary_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    for (g in 1:length(gene_name)) {
      #g <- 1
      gene_name0 <- gene_name[g]
      plot_data<-data.frame(imagerow=boundary_gene_exp$imagerow,imagecol=boundary_gene_exp$imagecol,boundary_gene_exp=boundary_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","boundary_gene_exp")
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=boundary_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = slice_name,x = "",y = "")
      
      pdf(paste0(dir_boundary_plot,"/",slice_name,"_boundary_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
    }
  }
  
  print(slice_name)
}  
#########budding#####
dir_budding_plot <- "/plot_budding"
for(i in 1:length(dataset_slice)){
  #i <- 1
  dataset_he_slice<-file_bdy[i]
  dataset_name<-str_split(dataset_he_slice,"_")[[1]][1]
  dataset_name<-str_split(dataset_name,"/")[[1]]
  slice_name <- paste0(dataset_name[1],"_",dataset_name[2],"_",dataset_name[3])

  dispersion_gene_exp <- read.delim(paste0(dir_exp,slice_name,"_dispersiongene.txt"),stringsAsFactors = F,check.names = F)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
  dispersion_gene_exp <- removeColsAllNa(dispersion_gene_exp)
  if(dim(dispersion_gene_exp)[2]>6){
    gene_data <- dispersion_gene_exp[,-c(1:5)]
    gene_name <- colnames(gene_data)
    for (g in 1:length(gene_name)) {
      #g <- 1
      gene_name0 <- gene_name[g]
      plot_data<-data.frame(imagerow=dispersion_gene_exp$imagerow,imagecol=dispersion_gene_exp$imagecol,dispersion_gene_exp=dispersion_gene_exp[gene_name0])
      colnames(plot_data) <- c("imagerow","imagecol","dispersion_gene_exp")
      
      
      p1<-ggplot(plot_data, aes(x=imagerow, y=imagecol)) +
        geom_point(aes(colour=dispersion_gene_exp),size=.6) +
        scale_color_gradientn(colours = c(colorRampPalette(c("#393b84","#3a6f8a"))(25),
                                          colorRampPalette(c("#3a6f8a","#3b9283"))(25),
                                          colorRampPalette(c("#3b9283","#6ab062"))(25),
                                          colorRampPalette(c("#6ab062","#b8d32f"))(25))
        )+ #设置填充颜色
        theme_classic()+
        labs(title = slice_name,x = "",y = "")
      
      pdf(paste0(dir_budding_plot,"/",slice_name,"_budding_",gene_name0,".pdf"),width = 4.8, height = 4)
      print(p1)
      dev.off()
      
    }
    
  }
  
  print(slice_name)
  
}  


