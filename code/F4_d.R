####不同MP的bdy的step1-5邻居的细胞类型组成
library(OmicCircos)
library(S4Vectors)
library(tidyverse)
library(magrittr)
library(circlize)
library(ComplexHeatmap)
library(dendextend)
library(reshape2)


dir_MP<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/'
file_MP<-list.files(pattern = '_MP_score.txt',path = dir_MP,recursive = T)
length(file_MP)
file_MP[1:10]
dataSlice<-unlist(lapply(strsplit(file_MP,'_MP_score'),function(x)x[1]))
dataSlice[1:10]
cancer<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[1]))
dataset<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[2]))
slice<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[3]))

dir_near<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_hop10/'
file_near<-paste0(dataSlice,'_nearSpotStep1to10.txt')
file_near[1:10]

dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/RCTD3/'
file_RCTD<-paste0(dataSlice,'_Deconvolution.txt')
file_RCTD[1:10]

dir_RCTD<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_Cytospace/'
file_RCTD<-list.files(pattern = '_cellFreq.txt',path = dir_RCTD,recursive = T)
length(file_RCTD)
file_RCTD[1:10]
dataSlice<-unlist(lapply(strsplit(file_RCTD,'_cellFreq'),function(x)x[1]))
dataSlice[1:10]
cancer<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[1]))
dataset<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[2]))
slice<-unlist(lapply(strsplit(dataSlice,'/'),function(x)x[3]))

dir_MP<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_NMF/MPscore_new/'
dir_near<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_hop10/'
file_MP<-paste0(dataSlice,'_MP_score.txt')
file_near<-paste0(dataSlice,'_nearSpotStep1to10.txt')


dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstepRCTD/'
dir_out<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstepCytoSpace/'

for(i in 1:length(dataSlice)){
  # i=1
  st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
  spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
  spot_step1<-intersect(spot_step1,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
  spot_step2<-intersect(spot_step2,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
  spot_step3<-intersect(spot_step3,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
  spot_step4<-intersect(spot_step4,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
  spot_step5<-intersect(spot_step5,st_near$cell_name[st_near$FinalLocalType%in%c('Immune','Normal')])
  
  spot_step5<-setdiff(spot_step5,spot_step4)
  spot_step4<-setdiff(spot_step4,spot_step3)
  spot_step3<-setdiff(spot_step3,spot_step2)
  spot_step2<-setdiff(spot_step2,spot_step1)
  spot_list<-list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5)
  names(spot_list)<-c('step_1','step_2','step_3','step_4','step_5')
  
  st_MP<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  st_MP<-st_MP[rownames(st_near),]
  st_MP$FinalLocalType<-st_near$FinalLocalType
  st_MP<-st_MP[st_MP$FinalLocalType%in%'Boundary',]
  st_MP[1:3,]
  table(st_MP$MP_top1)
  
  st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD[i]),stringsAsFactors = F,check.names = F)
  
  if(nrow(st_MP)>10){
    bdy_MP<-unique(st_MP$MP_top1)
    
    for(j in 1:length(bdy_MP)){#j=1
      MP_spot<-rownames(st_MP)[which(st_MP$MP_top1==bdy_MP[j])]
      
      for(nn in c('step_1','step_5')){#nn='step_1'
        MP_spot_step<-unique(intersect(unlist(spot_list[nn]),unlist(strsplit(st_near[MP_spot,nn],','))))
        MP_spot_step<-intersect(rownames(st_RCTD),MP_spot_step)
        if(length(MP_spot_step)>0){
          write.table(st_RCTD[MP_spot_step,],
                      paste0(dir_out,gsub('/','_',dataSlice[i]),'_',bdy_MP[j],'_',nn,'.txt'),
                      quote = F,sep = '\t')
        }
        
      }
      
    }
    
  }
  
}


####绘制环形热图##################
###RCTD
dir_re<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstepRCTD/'
file_step1<-list.files(pattern = 'step_1',path = dir_re,recursive = T)
# file_step2<-list.files(pattern = 'step_2_',path = dir_re,recursive = T)
# file_step3<-list.files(pattern = 'step_3_',path = dir_re,recursive = T)
# file_step4<-list.files(pattern = 'step_4_',path = dir_re,recursive = T)
file_step5<-list.files(pattern = 'step_5',path = dir_re,recursive = T)

file_step<-list(file_step1=file_step1,
                # file_step2=file_step2,
                # file_step3=file_step3,
                # file_step4=file_step4,
                file_step5=file_step5)
step_n<-c('step1','step5')##'step2','step3','step4',


plot_col=colorRamp2(seq(0,0.4,length.out=5),c("#463380",'#365d8f',"#6dc060",'#85c551',"#bfd631"))##"#133868",'#238DC1',
#plot_col=colorRamp2(seq(0,1,length.out=3),c("white",'#D45848',"#6C0E24"))##"#133868",'#238DC1',

MP_col=c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1",
         "MP_5"="#6A8EC9","MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861",
         "MP_9"="#652884","MP_10"="#444577","MP_11"="#8A7355","MP_12"="#B3BB61",
         "MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
         "MP_16"="#399335",'MP_17'='#96DD88')
#names(MP_col)<-unique(cir_data$MP)

dend_col <- structure(MP_col, names = 1:17)

dir_pic<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstep_pic/RCTD/'

#all_step_MP<-c()
for(j in 1:2){#j=2
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
  # allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()
  
  
  allSliceMP$slice_MP<-paste0(allSliceMP$slice,'of',allSliceMP$MP)
  allSliceMP$step<-step_n[j]
  allSliceMP[1:3,]
  allSliceMP<-allSliceMP[which(allSliceMP$celltype!='unknown'),]
  #all_step_MP<-rbind(all_step_MP,allSliceMP)
  
  # unique(file_re) %>% length()
  # length(unique(allSliceMP$slice_MP))
  cir_data<-as.data.frame(reshape2::acast(allSliceMP[,c('slice_MP','celltype','value')],slice_MP~celltype,mean))
  cir_data[1:3,]
  # ?acast
  
  cir_data<-cir_data[,c('CAF','TAM',setdiff(colnames(cir_data),c('CAF','TAM')))]
  cir_data<-cir_data[,setdiff(colnames(cir_data),c('Endothelial','Epithelial'))]
  
  cir_data$MP<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[2]))
  cir_data$slice<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[1]))
  cir_data$cancer<-unlist(lapply(strsplit(cir_data$slice,'_'),function(x)x[1]))
  # cir_data$cancer<-substr(cir_data$cancer,1,nchar(cir_data$cancer)-2) %>% toupper()
  
  
  cir_data$order<-unlist(lapply(strsplit(cir_data$MP,'_'),function(x) x[2]))%>%as.numeric()
  cir_data<-cir_data[order(cir_data$cancer),]
  cir_data<-cir_data[order(cir_data$order),]
  #range(cir_data[,1:14],na.rm = T)
  lev_split2<-factor(cir_data$order)
  
  pdf(paste0(dir_pic,step_n[j],'_cir_celltype_MP_heatmap.pdf'),width = 6,height = 6)
  circos.clear()
  #circos.par(gap.degree=1)
  #circos.par(gap.after = c(rep(1,12), 10))
  circos.par(gap.after = c(rep(1,16),18))
  #?circos.par
  circos.heatmap(cir_data[,c(1:12)], col = plot_col, na.col = "white",####去除Endothelial  Epithelial
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




###CytoSpace
dir_re<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstepCytoSpace//'
file_step1<-list.files(pattern = 'step_1',path = dir_re,recursive = T)
file_step5<-list.files(pattern = 'step_5',path = dir_re,recursive = T)

file_step<-list(file_step1=file_step1,
                file_step5=file_step5)
step_n<-c('step1','step5')


plot_col=colorRamp2(seq(0,0.3,length.out=5),c("#463380",'#365d8f',"#6dc060",'#85c551',"#bfd631"))##"#133868",'#238DC1',
#plot_col=colorRamp2(seq(0,1,length.out=3),c("white",'#D45848',"#6C0E24"))##"#133868",'#238DC1',

MP_col=c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1",
         "MP_5"="#6A8EC9","MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861",
         "MP_9"="#652884","MP_10"="#444577","MP_11"="#8A7355","MP_12"="#B3BB61",
         "MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
         "MP_16"="#399335",'MP_17'='#96DD88')
#names(MP_col)<-unique(cir_data$MP)

dend_col <- structure(MP_col, names = 1:17)

dir_pic<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_bdyMPstep_pic/CytoSpace/'
subCAFTAM<-c('apCAF','iCAF',"ifnCAF","rCAF","mCAF","vCAF","dCAF","tCAF",
             "blood_TAM","Microglial_TAM",'C1Q_TAM','FCN1_TAM',"SPP1_TAM")

for(j in 1:2){#j=2
  file_re<-file_step[[j]]
  
  re_slice<-unlist(lapply(strsplit(file_re,'_M'),function(x)x[1]))
  re_MP<-paste0('M',unlist(lapply(strsplit(file_re,'_M'),function(x)x[2])))
  re_MP<-unlist(lapply(strsplit(re_MP,'_step'),function(x)x[1]))
  
  allSliceMP<-c()
  for(i in 1:length(file_re)){
    #i=1
    tryCatch({
      re_data<-read.delim(paste0(dir_re,file_re[i]),stringsAsFactors = F,check.names = F,row.names = 1)
      #apply(re_data,2,mean)
      allSliceMP<-rbind(allSliceMP,
                        data.frame(celltype=colnames(re_data),
                                   value=apply(re_data,2,mean),
                                   slice=re_slice[i],
                                   MP=re_MP[i]))
    },error = function(e){
      print(file_re[i])
    })
    
  }
  allSliceMP<-allSliceMP[which(allSliceMP$value%in%NaN==F),]
  allSliceMP$cancer<-unlist(lapply(strsplit(allSliceMP$slice,'_'),function(x)x[1]))
  # allSliceMP$cancer<-substr(allSliceMP$cancer,1,nchar(allSliceMP$cancer)-2) %>% toupper()

  allSliceMP$slice_MP<-paste0(allSliceMP$slice,'of',allSliceMP$MP)
  allSliceMP$step<-step_n[j]
  allSliceMP<-allSliceMP[which(allSliceMP$celltype!='unknown'),]

  allSliceMP<-allSliceMP[allSliceMP$celltype%in%subCAFTAM,]
  #all_step_MP<-rbind(all_step_MP,allSliceMP)
  allSliceMP[1:3,]
  
  # unique(file_re) %>% length()
  # length(unique(allSliceMP$slice_MP))
  cir_data<-as.data.frame(reshape2::acast(allSliceMP[,c('slice_MP','celltype','value')],slice_MP~celltype,mean))
  cir_data[1:3,]
  
  cir_data<-cir_data[,c("mCAF",'iCAF',"tCAF","ifnCAF",'apCAF',"rCAF","vCAF","dCAF",
                        'C1Q_TAM',"SPP1_TAM",'FCN1_TAM',"blood_TAM","Microglial_TAM")]
  #range(cir_data[,1:13],na.rm = T)
  cir_data$MP<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[2]))
  cir_data$slice<-unlist(lapply(strsplit(rownames(cir_data),'of'),function(x)x[1]))
  cir_data$cancer<-unlist(lapply(strsplit(cir_data$slice,'_'),function(x)x[1]))

  cir_data$order<-unlist(lapply(strsplit(cir_data$MP,'_'),function(x) x[2]))%>%as.numeric()
  cir_data<-cir_data[order(cir_data$cancer),]
  cir_data<-cir_data[order(cir_data$order),]
  #range(cir_data[,1:14],na.rm = T)
  lev_split2<-factor(cir_data$order)
  
  pdf(paste0(dir_pic,step_n[j],'_cir_subCAFTAM_bdyMP_heatmap.pdf'),width = 6,height = 6)
  circos.clear()
  #circos.par(gap.degree=1)
  #circos.par(gap.after = c(rep(1,12), 10))
  circos.par(gap.after = c(rep(1,16),18))
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
  lgd <- Legend(title = paste0(step_n[j],"_CytoSpace"), col_fun = plot_col)
  grid.draw(lgd)
  dev.off()
  circos.clear()
}




















