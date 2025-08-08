####CAF TAM亚型展示
####绘制成环形热图
####不同MP类型的core的bdy邻居的多步ImmNor邻居的细胞类型比例
library(OmicCircos)
library(S4Vectors)
library(tidyverse)

library(tidyverse)
library(magrittr)
library(circlize)
library(ComplexHeatmap)
library(dendextend)

subCAFTAM<-c('apCAF','iCAF',"ifnCAF","rCAF","mCAF","vCAF","dCAF","tCAF",
             "blood_TAM","Microglial_TAM",'C1Q_TAM','FCN1_TAM',"SPP1_TAM")
####step1-5

####core的bdy的邻居
dir_re<-'E:/Mirror/ST_analysis/data/pic_data/bdy_step1_MP/CytoSpace/immNormal/'
file_step1<-list.files(pattern = 'step_1_',path = dir_re,recursive = T)
file_step2<-list.files(pattern = 'step_2_',path = dir_re,recursive = T)
file_step3<-list.files(pattern = 'step_3_',path = dir_re,recursive = T)
file_step4<-list.files(pattern = 'step_4_',path = dir_re,recursive = T)
file_step5<-list.files(pattern = 'step_5_',path = dir_re,recursive = T)

file_step<-list(file_step1=file_step1,
                file_step2=file_step2,
                file_step3=file_step3,
                file_step4=file_step4,
                file_step5=file_step5)
step_n<-c('step1','step2','step3','step4','step5')


plot_col=colorRamp2(seq(0,0.3,length.out=5),c("#463380",'#365d8f',"#6dc060",'#85c551',"#bfd631"))##"#133868",'#238DC1',
#plot_col=colorRamp2(seq(0,1,length.out=3),c("white",'#D45848',"#6C0E24"))##"#133868",'#238DC1',

MP_col=c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
         "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
         "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
         'MP_14'='#FF9694')
#names(MP_col)<-unique(cir_data$MP)

dend_col <- structure(MP_col, names = 2:14)

dir_pic<-'E:/Mirror/ST_analysis/pic/subCAFTAM/cir_heatmap/coreMP_bdy_step/'
for(j in 1:5){#j=1
  file_re<-file_step[[j]]
  
  re_slice<-unlist(lapply(strsplit(file_re,'_M'),function(x)x[1]))
  re_MP<-paste0('M',unlist(lapply(strsplit(file_re,'_M'),function(x)x[2])))
  re_MP<-unlist(lapply(strsplit(re_MP,'_step'),function(x)x[1]))
  
  allSliceMP<-c()
  for(i in 1:length(file_re)){
    #i=1
    re_data<-read.delim(paste0(dir_re,file_re[i]),stringsAsFactors = F,check.names = F,row.names = 1)
    #apply(re_data,2,mean)
    allSliceMP<-rbind(allSliceMP,
                      data.frame(celltype=colnames(re_data),
                                 value=apply(re_data,2,mean),
                                 slice=re_slice[i],
                                 MP=re_MP[i]))
  }
  allSliceMP<-allSliceMP[which(allSliceMP$value%in%NaN==F),]
  allSliceMP$cancer<-unlist(lapply(strsplit(allSliceMP$slice,'_'),function(x)x[1]))
  allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()
  
  
  allSliceMP$slice_MP<-paste0(allSliceMP$slice,'of',allSliceMP$MP)
  allSliceMP$step<-step_n[j]
  
  allSliceMP<-allSliceMP[allSliceMP$celltype%in%subCAFTAM,]
  #all_step_MP<-rbind(all_step_MP,allSliceMP)
  
  # unique(file_re) %>% length()
  # length(unique(allSliceMP$slice_MP))
  cir_data<-as.data.frame(reshape2::acast(allSliceMP[,c('slice_MP','celltype','value')],slice_MP~celltype))
  
  cir_data<-cir_data[,c("mCAF",'iCAF',"tCAF","ifnCAF",'apCAF',"rCAF","vCAF","dCAF",
                        'C1Q_TAM',"SPP1_TAM",'FCN1_TAM',"blood_TAM","Microglial_TAM")]
  #range(cir_data[,1:13],na.rm = T)
  cir_data$MP<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[2]))
  cir_data$slice<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[1]))
  cir_data$cancer<-unlist(lapply(strsplit(cir_data$slice,'_'),function(x)x[1]))
  cir_data$cancer<-substr(cir_data$cancer,1,nchar(cir_data$cancer)-2) %>% toupper()
  
  
  cir_data$order<-unlist(lapply(strsplit(cir_data$MP,'_'),function(x) x[2]))%>%as.numeric()
  cir_data<-cir_data[order(cir_data$cancer),]
  cir_data<-cir_data[order(cir_data$order),]
  #range(cir_data[,1:14],na.rm = T)
  lev_split2<-factor(cir_data$order)
  
  pdf(paste0(dir_pic,step_n[j],'_cir_subCAFTAM_coreMP_heatmap.pdf'),width = 6,height = 6)
  circos.clear()
  #circos.par(gap.degree=1)
  #circos.par(gap.after = c(rep(1,12), 10))
  circos.par(gap.after = c(rep(1,12),18))
  #?circos.par
  circos.heatmap(cir_data[,c(1:13)], col = plot_col, na.col = "white",
                 track.height = 0.4,###热图条带的宽度
                 track.margin=c(0.02,0),###前面的值是与向内的下一个热图轨迹的间隔，后面的值是与向外的上一个热图轨迹的间隔
                 #rownames.side = "outside",
                 cluster = T,dend.side = 'inside',dend.track.height = 0.2,
                 dend.callback = function(dend, m, si) {
                   color_branches(dend, k = 15, col = dend_col[si])
                 },
                 cell.border=NA,
                 cell.lwd=0.02,
                 split = lev_split2 , 
                 rownames.cex = 0.5)
  lgd <- Legend(title = paste0(step_n[j],"_RCTD"), col_fun = plot_col)
  grid.draw(lgd)
  dev.off()
  circos.clear()
}


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
####bdy的邻居
dir_re<-'E:/Mirror/ST_analysis/data/pic_data/bdy_MP/CytoSpace/immNormal/'
file_step1<-list.files(pattern = 'step_1_',path = dir_re,recursive = T)
file_step2<-list.files(pattern = 'step_2_',path = dir_re,recursive = T)
file_step3<-list.files(pattern = 'step_3_',path = dir_re,recursive = T)
file_step4<-list.files(pattern = 'step_4_',path = dir_re,recursive = T)
file_step5<-list.files(pattern = 'step_5_',path = dir_re,recursive = T)

file_step<-list(file_step1=file_step1,
                file_step2=file_step2,
                file_step3=file_step3,
                file_step4=file_step4,
                file_step5=file_step5)
step_n<-c('step1','step2','step3','step4','step5')


plot_col=colorRamp2(seq(0,0.3,length.out=5),c("#463380",'#365d8f',"#6dc060",'#85c551',"#bfd631"))##"#133868",'#238DC1',
#plot_col=colorRamp2(seq(0,1,length.out=3),c("white",'#D45848',"#6C0E24"))##"#133868",'#238DC1',

MP_col=c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
         "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
         "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
         'MP_14'='#FF9694')
#names(MP_col)<-unique(cir_data$MP)

dend_col <- structure(MP_col, names = 2:14)

dir_pic<-'E:/Mirror/ST_analysis/pic/subCAFTAM/cir_heatmap/bdyMP_step/'
for(j in 1:5){#j=1
  file_re<-file_step[[j]]
  
  re_slice<-unlist(lapply(strsplit(file_re,'_M'),function(x)x[1]))
  re_MP<-paste0('M',unlist(lapply(strsplit(file_re,'_M'),function(x)x[2])))
  re_MP<-unlist(lapply(strsplit(re_MP,'_step'),function(x)x[1]))
  
  allSliceMP<-c()
  for(i in 1:length(file_re)){
    #i=1
    re_data<-read.delim(paste0(dir_re,file_re[i]),stringsAsFactors = F,check.names = F,row.names = 1)
    #apply(re_data,2,mean)
    allSliceMP<-rbind(allSliceMP,
                      data.frame(celltype=colnames(re_data),
                                 value=apply(re_data,2,mean),
                                 slice=re_slice[i],
                                 MP=re_MP[i]))
  }
  allSliceMP<-allSliceMP[which(allSliceMP$value%in%NaN==F),]
  allSliceMP$cancer<-unlist(lapply(strsplit(allSliceMP$slice,'_'),function(x)x[1]))
  allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()
  
  
  allSliceMP$slice_MP<-paste0(allSliceMP$slice,'of',allSliceMP$MP)
  allSliceMP$step<-step_n[j]
  
  allSliceMP<-allSliceMP[allSliceMP$celltype%in%subCAFTAM,]
  #all_step_MP<-rbind(all_step_MP,allSliceMP)
  
  # unique(file_re) %>% length()
  # length(unique(allSliceMP$slice_MP))
  cir_data<-as.data.frame(reshape2::acast(allSliceMP[,c('slice_MP','celltype','value')],slice_MP~celltype))
  
  cir_data<-cir_data[,c("mCAF",'iCAF',"tCAF","ifnCAF",'apCAF',"rCAF","vCAF","dCAF",
                        'C1Q_TAM',"SPP1_TAM",'FCN1_TAM',"blood_TAM","Microglial_TAM")]
  #range(cir_data[,1:13],na.rm = T)
  cir_data$MP<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[2]))
  cir_data$slice<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[1]))
  cir_data$cancer<-unlist(lapply(strsplit(cir_data$slice,'_'),function(x)x[1]))
  cir_data$cancer<-substr(cir_data$cancer,1,nchar(cir_data$cancer)-2) %>% toupper()
  
  
  cir_data$order<-unlist(lapply(strsplit(cir_data$MP,'_'),function(x) x[2]))%>%as.numeric()
  cir_data<-cir_data[order(cir_data$cancer),]
  cir_data<-cir_data[order(cir_data$order),]
  #range(cir_data[,1:14],na.rm = T)
  lev_split2<-factor(cir_data$order)
  
  pdf(paste0(dir_pic,step_n[j],'_cir_subCAFTAM_bdyMP_heatmap.pdf'),width = 6,height = 6)
  circos.clear()
  #circos.par(gap.degree=1)
  #circos.par(gap.after = c(rep(1,12), 10))
  circos.par(gap.after = c(rep(1,12),18))
  #?circos.par
  circos.heatmap(cir_data[,c(1:13)], col = plot_col, na.col = "white",
                 track.height = 0.4,###热图条带的宽度
                 track.margin=c(0.02,0),###前面的值是与向内的下一个热图轨迹的间隔，后面的值是与向外的上一个热图轨迹的间隔
                 #rownames.side = "outside",
                 cluster = T,dend.side = 'inside',dend.track.height = 0.2,
                 dend.callback = function(dend, m, si) {
                   color_branches(dend, k = 15, col = dend_col[si])
                 },
                 cell.border=NA,
                 cell.lwd=0.02,
                 split = lev_split2 , 
                 rownames.cex = 0.5)
  lgd <- Legend(title = paste0(step_n[j],"_RCTD"), col_fun = plot_col)
  grid.draw(lgd)
  dev.off()
  circos.clear()
}