###8中T细胞状态
library(TCellSI)
library(dplyr)
library(tidyverse)
library(magrittr)
library(Seurat)
library(monocle)
library(igraph)
library(RobustRankAggreg)
library(jsonlite)


dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
sample_group<-'1_2_9_10_12_public'
file_RCTD<-list.files(pattern = paste0('Deconvolution_',sample_group),path = dir_bdy,recursive = T) 
# Deconvolution_CAFTAM_ESCC_integrate.txt
# Deconvolution_CAFTAM_public.txt
# Deconvolution_CAFTAM.txt

dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
patient<-unlist(lapply(strsplit(file_rds,'/'),function(x)x[1]))


i=2
st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
table(st_bdy$FinalLocalType)
st_RCTD<-read.delim(paste0(dir_bdy,file_RCTD[i]),stringsAsFactors = F,check.names = F)
st_RCTD<-st_RCTD[,setdiff(colnames(st_RCTD),c("Endothelial","Epithelial"))]
st_bdy<-st_bdy[rownames(st_RCTD),]
st_RCTD$celltype<-apply(st_RCTD,1,function(x){
  names(x)<-colnames(st_RCTD)
  x<-x[order(x,decreasing = T)]
  return(names(x)[1])
})
st_RCTD$FinalLocalType<-st_bdy$FinalLocalType
table(st_RCTD$celltype)
st_RCTD<-st_RCTD[which(st_RCTD$celltype=='T lymphocytes'),]
T_spot<-rownames(st_RCTD)[which(st_RCTD$FinalLocalType=='Immune'|st_RCTD$FinalLocalType=='Normal')]

tryCatch({
  if(length(T_spot)>=5){
    st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
    st_count<-st_rds@assays[["Spatial"]]@counts[,T_spot] %>% as.matrix()
    st_scores <- TCSS_Calculate(st_count) #%>% t()
    st_scores<-t(st_scores)
    write.table(st_scores,paste0(dir_bdy,patient[i],'/slice1_TCellSI_',sample_group,'.txt'),quote = F,sep = '\t')
  }
  print(patient[i])
},error = function(e){
  stopMessage_aaa<-"Unable to convert data"
  print(paste0(patient[i],'_stop'))
})

####比较 Quiescence 得分在两个切片里的高低 画箱式图
p_data<-all_step_Tcell[all_step_Tcell$Var2%in%'Quiescence',]
TcellSI_name<-unique(all_step_Tcell$Var2)

pdf('E:/Mirror/ST_analysis/pic/ESCC/TcellSI/TcellSI_compare_1_2_9_10_12_public.pdf',width = 8,height = 4)
for(j in TcellSI_name){
  all_step_Tcell<-c()
  for(i in 1:2){
    #i=1
    st_near<-read.delim(paste0(dir_near,file_near[i]),stringsAsFactors = F,check.names = F)
    st_Tcell<-read.delim(paste0(dir_near,file_Tcell[i]),stringsAsFactors = F,check.names = F)
    
    spot_step1<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_1'],',')) %>% unique()
    length(intersect(spot_step1,rownames(st_Tcell)))
    spot_step2<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_2'],',')) %>% unique()
    spot_step2<-union(spot_step1,spot_step2)
    spot_step3<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_3'],',')) %>% unique()
    spot_step3<-union(spot_step2,spot_step3)
    spot_step4<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_4'],',')) %>% unique()
    spot_step4<-union(spot_step3,spot_step4)
    spot_step5<-unlist(strsplit(st_near[st_near$FinalLocalType%in%'Boundary','step_5'],',')) %>% unique()
    spot_step5<-union(spot_step4,spot_step5)
    
    step_spot<-list(spot_step1,spot_step2,spot_step3,spot_step4,spot_step5)
    names(step_spot)<-paste0('step_',1:5)
    
    step_Tcell<-lapply(1:5,function(x){
      spot_Tcell<-intersect(step_spot[[x]],rownames(st_Tcell))
      if(length(spot_Tcell)>1){
        near_Tcell<-data.frame(step=names(step_spot)[x],
                               score=st_Tcell[spot_Tcell,j])
      }
      return(near_Tcell)
    })
    step_Tcell<-do.call(rbind,step_Tcell)
    step_Tcell$patient<-patient[i]
    
    all_step_Tcell<-rbind(all_step_Tcell,step_Tcell)
  }
  
  #all_step_Tcell$patient<-unlist(lapply(strsplit(all_step_Tcell$slice,'/'),function(x)x[1]))
  #all_step_Tcell$cancer<-substr(all_step_Tcell$cancer,1,nchar(all_step_Tcell$cancer)-2) %>% toupper()
  all_step_Tcell<-all_step_Tcell[which(all_step_Tcell$score!=0),]
  
  
  p_box<-ggplot(all_step_Tcell, aes(x = step, y = score,fill=patient))+ 
    # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
    stat_boxplot(geom = 'errorbar',width=0.4,position = position_dodge(0.9))+
    geom_boxplot(aes(fill = patient), color='black',width = 0.7,lwd=0.1,#fatten=0.9,
                 position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    theme_classic(base_size = 12)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
    scale_color_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
    ggtitle(j)+
    # theme(plot.title = element_text(hjust = 0.4))+
    # theme(plot.title = element_text(size = 12))+
    stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =max(all_step_Tcell$score))+
    #stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
    #theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  print(p_box)
  
}
dev.off()


