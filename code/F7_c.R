#######不同step的细胞类型的变化
library(tidyverse)
library(Seurat)




dir_step<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_step<-list.files(pattern = 'nearSpotStep1to20.txt',path = dir_step,recursive = T)
sample_group<-'1_2_9_10_12_public'
file_Deco<-list.files(pattern = paste0('Deconvolution_',sample_group,'.txt'),path = dir_step,recursive = T)
#file_Deco<-file_Deco[c(1:3,67,164)]
dataSlice<-unlist(lapply(strsplit(file_step,'_'),function(x) x[1]))
table(dataSlice==unlist(lapply(strsplit(file_Deco,'_'),function(x) x[1])))
dataset<-unlist(lapply(strsplit(dataSlice,'/'),function(x) x[1]))

dir_out<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'

for(i in 1:length(file_Deco)){
  #i=1
  #if(!dir.exists(paste0(dir_out,dataset[i]))) dir.create(paste0(dir_out,dataset[i]))
  
  st_step<-read.delim(paste0(dir_step,file_step[i]),stringsAsFactors = F,check.names = F)
  st_Deco<-read.delim(paste0(dir_step,file_Deco[i]),stringsAsFactors = F,check.names = F)
  colnames(st_step)
  #st_Deco<-st_Deco[,setdiff(colnames(st_Deco),c('Endothelial','Epithelial'))]
  st_step<-st_step[,c("cell.names","imagerow","imagecol","FinalLocalType",paste0('step_',1:20))]
  st_step<-st_step[rownames(st_Deco),]
  
  #st_step<-st_step[which(st_step$FinalLocalType=='Normal'|st_step$FinalLocalType=='Immune'),]
  st_Deco<-st_Deco[rownames(st_step)[st_step$FinalLocalType%in%c('Normal','Immune')],]
  
  st_step_CAF<-st_step[st_step$FinalLocalType%in%c('Boundary','Dispersion'),]
  
  for(jj in ncol(st_step_CAF):6){#jj=6
    y_spot<-unlist(strsplit(st_step_CAF[,(jj-1)],',')) %>% unique()
    st_step_CAF[,jj]<-lapply(1:nrow(st_step_CAF),function(x){##x=1
      x_spot<-strsplit(st_step_CAF[x,jj],',') %>% unlist()
      x_spot<-setdiff(x_spot,y_spot)
      return(paste0(x_spot,collapse = ','))
    }) %>% unlist()
  }
  del_col<-which(apply(st_step_CAF,2,function(x){
    length(which(x==''))
  })==nrow(st_step_CAF))
  if(length(del_col)>0) st_step_CAF<-st_step_CAF[,-del_col]
  
  step_cell<-c()
  for(ss in 5:ncol(st_step_CAF)){##ss=9
    step_spot<-unlist(strsplit(st_step_CAF[,ss],',')) %>% unique()
    
    if(length(intersect(rownames(st_Deco),step_spot))>0){
      step_Deco<-st_Deco[intersect(rownames(st_Deco),step_spot),]
      #add_cell<-apply(step_Deco,2,mean)
      add_cell<-reshape2::melt(as.matrix(step_Deco))
      add_cell$step<-paste0('step_',(ss-4))
      step_cell<-rbind(step_cell,add_cell)
    }
    
  }
  
  write.table(step_cell,paste0(dir_out,dataSlice[i],'_RCTD_step_',sample_group,'.txt'),quote = F,sep = '\t',row.names = F)
}



dir_step<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_RCTDstep<-list.files(pattern = paste0('_RCTD_step_',sample_group,'.txt'),path = dir_step,recursive = T)
patient<-unlist(lapply(strsplit(file_RCTDstep,'/'),function(x) x[1]))


all_RCTD_step<-c()
for(i in 1:2){
  RCTD_step<-read.delim(paste0(dir_step,file_RCTDstep[i]),stringsAsFactors = F,check.names = F)
  RCTD_step$patient<-patient[i]
  all_RCTD_step<-rbind(all_RCTD_step,RCTD_step)
}


#### CAF 比较
p_data<-all_RCTD_step[all_RCTD_step$Var2%in%"CAF",]
table(p_data$step,p_data$patient)

p_data<-p_data[p_data$step%in%paste0('step_',1:6),]
p_box<-aggregate(p_data$value,by=list(p_data$patient,p_data$step),mean)
colnames(p_box)<-c('patient','step','value')
p_box<-p_box[order(p_box$patient),]

library(ggpubr)
ggpaired(p_box, 
         x = "patient", 
         y = "value",
         color = "patient", 
         line.color = "gray", 
         line.size = 0.4,
         palette = "jco") +
  stat_compare_means(paired = TRUE)


p_box<-ggplot(p_data, aes(x = step, y = value,fill=patient))+ 
  # geom_violin(aes(color = near), trim = T,position = position_dodge(0.8),alpha=0.6) +
  stat_boxplot(geom = 'errorbar',width=0.4,position = position_dodge(0.9))+
  geom_boxplot(aes(fill = patient), color='black',width = 0.7,lwd=0.1,#fatten=0.9,
               position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
  theme_classic(base_size = 12)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
  scale_color_manual(values = c("BGM"="#d2e7cf","XYZ"="#179192"))+
  ggtitle('CAF')+
  # theme(plot.title = element_text(hjust = 0.4))+
  # theme(plot.title = element_text(size = 12))+
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),label.x = 2,label.y =0.68)+
  #stat_compare_means(aes(label = ..p.signif..),label.x = 2,label.y =0.6)+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )+
  coord_cartesian(ylim = c(0,0.7))
print(p_box)


pdf('E:/Mirror/ST_analysis/pic/ESCC/Deconvolution/CAF_compare_public.pdf',width = 9,height = 4)
print(p_box)
dev.off()

