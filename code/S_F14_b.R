####免疫受体的表达
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

imm_check<-read.delim('E:/Mirror/ST_analysis/other_data/imm-check-gene.txt',stringsAsFactors = F,check.names = F)

dir_rds<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/ST_expression/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)

dir_bdy<-'E:/Mirror/ST_analysis/data/ESCC/1re_data/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)
patient<-unlist(lapply(strsplit(file_bdy,'/'),function(x)x[1]))

all_immLigand<-c()
for(i in 1:length(file_bdy)){
  #i=1
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[st_bdy$FinalLocalType%in%c("Core",'Boundary','Dispersion'),]
  
  st_conut<-st_rds@assays[["Spatial"]]@counts[,st_bdy$cell.names]%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[imm_check$ligand,]
  rownames(st_conut)<-imm_check$ligand
  st_conut[is.na(st_conut)]<-0
  
  imm_score<-data.frame(cell_name=st_bdy$cell.names,
                        LocalType=st_bdy$FinalLocalType,
                        patient=patient[i]
  )
  imm_score<-cbind(imm_score,t(st_conut)) %>% as.data.frame()
  all_immLigand<-rbind(all_immLigand,imm_score)
}
#write.table()

all_exp_mean<-aggregate(all_immLigand[,4:16],by=list(all_immLigand$patient,all_immLigand$LocalType),mean)
all_exp_mean<-all_exp_mean[all_exp_mean$Group.2%in%'Core',]
rownames(all_exp_mean)<-all_exp_mean$Group.1
all_exp_mean<-all_exp_mean[,-c(1,2)]

all_p_value<-c()
for(i in patient){#i=patient[1]
  p_exp<-all_immLigand[all_immLigand$patient%in%i,]
  core_site<-which(p_exp$LocalType=='Core')
  BdyDis_site<-which(p_exp$LocalType=='Boundary'|p_exp$LocalType=='Dispersion')
  
  p_value<-apply(p_exp[,4:16],2,function(x){
    wilcox.test(x[core_site],x[BdyDis_site],alternative = 'greater')[["p.value"]]
  }) %>% t() %>% as.data.frame()
  rownames(p_value)<-i
  all_p_value<-rbind(all_p_value,p_value)
}


all_exp_mean<-t(all_exp_mean)
all_p_value<-t(all_p_value)

dis_dot<-ifelse(as.matrix(all_p_value) < 0.05, "*", "")
# dis_dot<-ifelse(as.matrix(all_p_value) < 0.01, "**", "")
# dis_dot<-ifelse(as.matrix(all_p_value) < 0.001, "***", "")
#diag(dis_dot)<-""
dis_dot[is.na(dis_dot)]<-""
range(all_exp_mean,na.rm = T)
#heat_data<-scale(all_exp_mean)
heat_data<-apply(t(all_exp_mean), 1,function(x){
  x<-x-min(x)
  x<-x/max(x)
})
range(heat_data)
bk<-seq(0,1,length.out=100)
# color_pheatmap<-colorRampPalette(c(#"#313B93","#4578B6",'#82B9D8','#C8E7F1',
#                                    "#F7F7C7",
#                                    '#FADB8C',"#F58B52",'#DA3627',"#A31D2A"))(255) ###"#CC281B"
color_pheatmap<-c(colorRampPalette(c("white",'#E9C1C6'))(5),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(95)) ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(heat_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = T,
                      cluster_cols = T,
                      treeheight_row = T,treeheight_col = T,
                      display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "white",###black
                      fontsize = 10,
                      cellwidth=15,
                      cellheight=15,
                      main = "Imm_ligand_CoreVsBdyDis",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/ESCC/ImmLigand_CoreVsBdyDis.pdf',width = 3.6,height = 3.7)
print(p)
dev.off()
