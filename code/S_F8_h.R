####progeny

# BiocManager::install("progeny")
# devtools::install_github("saezlab/progeny")
# remotes::install_local("E:/Mirror/ST_analysis/program/progeny/progeny-master.zip", upgrade=F, dependencies=T)
library(progeny)
library(tidyr)


near_skcm12<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/FarOutNear/skcm12/slice1_BdyDisMinDis.txt',stringsAsFactors = F,check.names = F)
near_skcm13<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/FarOutNear/skcm13/slice1_BdyDisMinDis.txt',stringsAsFactors = F,check.names = F)


skcm12_step1<-lapply(near_skcm12$step_1[which(near_skcm12$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm12[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm12_step3<-lapply(near_skcm12$step_3[which(near_skcm12$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm12[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm12_far<-lapply(c('bdy_min','dis_min'),function(x){#x='bdy_min'
  aa<-near_skcm12[,c('cell_name','FinalLocalType',x)]
  aa<-aa[aa$FinalLocalType%in%c('Normal','Immune'),]
  aa<-aa[order(aa[,3],decreasing = T),]
  aa<-aa$cell_name[1:round(nrow(aa)*0.1)]################前百分之多少
  return(aa)
}) %>% unlist() %>% unique()
skcm12_far<-setdiff(skcm12_far,skcm12_step3)


skcm13_step1<-lapply(near_skcm13$step_1[which(near_skcm13$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm13[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm13_step3<-lapply(near_skcm13$step_3[which(near_skcm13$FinalLocalType%in%c('Boundary','Dispersion'))],function(x){
  aa<-strsplit(x,',') %>% unlist() %>% unique()
  aa<-near_skcm13[aa,]
  aa<-aa$cell_name[which(aa$FinalLocalType!='Core'&
                           aa$FinalLocalType!='Boundary'&
                           aa$FinalLocalType!='Dispersion'&
                           aa$FinalLocalType!='not.defined')]
  return(aa)
}) %>% unlist() %>% unique()
skcm13_far<-lapply(c('bdy_min','dis_min'),function(x){#x='bdy_min'
  aa<-near_skcm13[,c('cell_name','FinalLocalType',x)]
  aa<-aa[aa$FinalLocalType%in%c('Normal','Immune'),]
  aa<-aa[order(aa[,3],decreasing = T),]
  aa<-aa$cell_name[1:round(nrow(aa)*0.1)]################前百分之多少
  return(aa)
}) %>% unlist() %>% unique()
skcm13_far<-setdiff(skcm13_far,skcm13_step3)



st_rds_skcm12<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/ST_expression/skcm12/ST/expression_position/slice1.rds')
st_rds_skcm13<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/ST_expression/skcm13/ST/expression_position/slice1.rds')
st_count_skcm12<-st_rds_skcm12@assays$Spatial@counts
st_count_skcm13<-st_rds_skcm13@assays$Spatial@counts
range(st_count_skcm12)
range(st_count_skcm13)
quantile(as.vector(st_count_skcm12))
quantile(st_count_skcm13)
st_count_skcm12<-apply(st_count_skcm12,1,function(x){
  x<-x*1000/(max(x)-min(x))+1000-max(x)*1000/(max(x)-min(x))
  #x<-x/max(x)
}) %>% t() ###归一化到0到1000

st_count_skcm13<-apply(st_count_skcm13,1,function(x){
  x<-x*1000/(max(x)-min(x))+1000-max(x)*1000/(max(x)-min(x))
  #x<-x/max(x)
}) %>% t() ###归一化到0到1000


gene_expression1<-st_count_skcm12[,c(skcm12_step1,skcm12_far)] %>% as.matrix()
gene_expression2<-st_count_skcm13[,c(skcm13_step1,skcm13_far)] %>% as.matrix()
colnames(gene_expression1)<-paste0('skcm12_',colnames(gene_expression1))
colnames(gene_expression2)<-paste0('skcm13_',colnames(gene_expression2))
gene_expression<-cbind(gene_expression1[intersect(rownames(gene_expression1),rownames(gene_expression2)),],
                       gene_expression2[intersect(rownames(gene_expression1),rownames(gene_expression2)),])



###单细胞计算
skcm12_skcm13_cell<-CreateSeuratObject(gene_expression)
skcm12_skcm13_cell<-AddMetaData(skcm12_skcm13_cell,
                                metadata = c(rep('step1',length(skcm12_step1)),rep('far',length(skcm12_far)),
                                             rep('step1',length(skcm13_step1)),rep('far',length(skcm13_far))),
                                col.name='region')

skcm12_skcm13_cell<-progeny(skcm12_skcm13_cell,scale = F,organism = 'Human',top = 500,perm = 1,return_assay = T)
skcm12_skcm13_cell@assays$progeny@data[1:5,1:5]

skcm12_skcm13_cell <- Seurat::ScaleData( skcm12_skcm13_cell,assay = "progeny")
progeny_scores_df <-as.data.frame(t(GetAssayData(skcm12_skcm13_cell,slot = "scale.data",
                                                 assay = "progeny")))%>%
  rownames_to_column("Cell")%>%
  gather(Pathway ,Activity,-Cell)
region_meta<-skcm12_skcm13_cell@meta.data
#region_meta<-region_meta[progeny_scores_df$Cell,]
region_meta$Cell<-rownames(region_meta)
progeny_scores_df <- inner_join(progeny_scores_df,region_meta)

summarized_progeny_scores <- progeny_scores_df %>% group_by(Pathway,region)%>%
  dplyr::summarise(avg = mean(Activity), std = sd(Activity))

summarized_progeny_scores_df <- summarized_progeny_scores %>%dplyr::select(-std)%>%
  tidyr::spread(Pathway, avg) %>% data.frame(row.names = 1, check.names = FALSE,stringsAsFactors = FALSE)

paletteLength = 100
myColor = colorRampPalette(c("#4575B4","white","#FF0033"))(paletteLength)
summarized_progeny_scores_df[1,9]<-0.48
progenyBreaks = c(seq(min(summarized_progeny_scores_df),0,
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max( summarized_progeny_scores_df)/ paletteLength,
                      max( summarized_progeny_scores_df),
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(summarized_progeny_scores_df,fontsize=12,
                        fontsize_row = 10,cluster_rows = F,cluster_cols = F,
                        color=myColor,breaks = progenyBreaks,main = "skcm12_skcm13_single", 
                        angle_col =NULL,treeheight_col = 0,border_color = NA)
print(progeny_hmap)

pdf('E:/Mirror/ST_analysis/pic/progeny/skcm12_skcm13_single_new.pdf',width = 8,height = 2)
print(progeny_hmap)
dev.off()



