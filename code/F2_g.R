############免疫受体在一个癌症的所有切片的core区的表达均值
##以及core区与非core恶性区的秩和检验显著性
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

imm_check<-read.delim('E:/Mirror/ST_analysis/other_data/imm-check-gene.txt',stringsAsFactors = F,check.names = F)



dir_rds<-'E:/Mirror/ST_analysis/data/10X Visium/ST_expression/'
file_rds<-list.files(pattern = '.rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('/ST/expression_position/','/',file_rds)
dataset_slice<-gsub('.rds','',dataset_slice)
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x) x[1]))
cancer<-substr(dataSet,1,nchar(dataSet)-2) %>% toupper()

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/copykat/'
file_bdy<-list.files(pattern = 'BdyTumorCore.txt',path = dir_bdy,recursive = T)

all_immLigand<-c()
for(i in 1:229){
  #i=1
  st_rds<-readRDS(paste0(dir_rds,file_rds[i]))
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType%in%c('Core','Boundary','Dispersion')),]
  
  st_conut<-st_rds@assays[["Spatial"]]@counts[,st_bdy$cell_name]%>%as.matrix()%>%as.data.frame()
  st_conut<-st_conut[imm_check$ligand,]
  rownames(st_conut)<-imm_check$ligand
  st_conut[is.na(st_conut)]<-0
  
  imm_score<-data.frame(cell_name=st_bdy$cell_name,
                        LocalType=st_bdy$FinalLocalType,
                        dataSlice=dataset_slice[i],
                        cancer=cancer[i])
  imm_score<-cbind(imm_score,t(st_conut)) %>% as.data.frame()
  all_immLigand<-rbind(all_immLigand,imm_score)
  print(dataset_slice[i])
}
write.table(all_immLigand,'E:/Mirror/ST_analysis/data/pic_data/all_immLigand_CoreBdyDis.txt',quote = F,sep = '\t')

all_immLigand<-read.delim('E:/Mirror/ST_analysis/data/pic_data/all_immLigand_CoreBdyDis.txt',stringsAsFactors = F,check.names = F)
apply(all_immLigand[,5:18],2,function(x){
  length(which(x==0))
})
all_immLigand<-all_immLigand[,-18]

cancer<-unique(all_immLigand$cancer)

all_exp_mean<-c()
all_p_value<-c()
for(i in 1:length(cancer)){#i=1
  cancer_Imm<-all_immLigand[which(all_immLigand$cancer==cancer[i]),]
  exp_mean<-aggregate(cancer_Imm[,5:17],by=list(cancer_Imm$LocalType),mean)
  exp_mean<-exp_mean[which(exp_mean$Group.1=='Core'),]
  rownames(exp_mean)<-cancer[i]
  exp_mean<-exp_mean[,-1]
  all_exp_mean<-rbind(all_exp_mean,exp_mean)
  
  core_site<-which(cancer_Imm$LocalType=='Core')
  BdyDis_site<-which(cancer_Imm$LocalType=='Boundary'|cancer_Imm$LocalType=='Dispersion')
  
  p_value<-apply(cancer_Imm[,5:17],2,function(x){
    wilcox.test(x[core_site],x[BdyDis_site],alternative = 'greater')[["p.value"]]
  }) %>% t() %>% as.data.frame()
  rownames(p_value)<-cancer[i]
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
heat_data<-scale(all_exp_mean)
heat_data<-apply(all_exp_mean, 1,function(x){
  x<-x-min(x)
  x<-x/max(x)
}) %>% t()
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
pdf('E:/Mirror/ST_analysis/pic/re/3/ImmLigand_CoreVsBdyDis.pdf',width = 5.5,height = 4)
print(p)
dev.off()









