###各区域marker更改
library(Seurat)
library(rlang)
library(ggplot2)
library(tidyverse)
library(ggraph)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(jsonlite)
library(stringr)




####还是区分不明显 再改
st_rds<-readRDS('E:/Mirror/ST_analysis/data/10X Visium/ST_expression/luad03/ST/expression_position/slice1.rds')
#st_count<-st_rds@assays$Spatial@counts
st_bdy<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/copykat/luad03/slice1_BdyTumorCore.txt',stringsAsFactors = F)
st_bdy<-st_bdy[colnames(st_rds),]
table(st_bdy$FinalLocalType)
region_marker<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/marker_SVG/luad03/slice1_marker.txt',
                          stringsAsFactors = F,check.names = F)
# CoreBdyDis_marker<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/marker_SVG/luad03/slice1_marker_CoreBdyDis.txt',
#                               stringsAsFactors = F,check.names = F)
# MalImm_marker<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/marker_SVG/luad03/slice1_marker_Mal_Imm.txt',
#                               stringsAsFactors = F,check.names = F)
table(region_marker$cluster)
region_marker<-region_marker[order(region_marker$avg_log2FC,decreasing = T),]
#region_marker<-region_marker[which(region_marker$p_val_adj<0.05),]
# table(CoreBdyDis_marker$cluster)
# CoreBdyDis_marker<-CoreBdyDis_marker[order(CoreBdyDis_marker$avg_log2FC,decreasing = T),]
# table(MalImm_marker$cluster)
# MalImm_marker<-MalImm_marker[order(MalImm_marker$avg_log2FC,decreasing = T),]

core_marker<-region_marker$gene[which(region_marker$cluster=='Core')]
bdy_marker<-region_marker$gene[which(region_marker$cluster=='Boundary')]
dis_marker<-region_marker$gene[which(region_marker$cluster=='Dispersion')]

imm_marker<-region_marker$gene[which(region_marker$cluster=='Immune')]
nor_marker<-region_marker$gene[which(region_marker$cluster=='Normal')]

stemgene<-c("FOXJ1","DDAH1","PCSK1N","TMEM106C","STEAP4","CUEDC1","CD24","GPC4","HMGB3","KRT8","MCM7",    
            "LSR","DMKN","IRX2","TACSTD2","ABCC3","ADAM15","MLPH","RCC1","ENO1","HMGN2","STMN1",   
            "MET","MUC1","POLD2","PTMA","MMP15","NDUFS6","PTPRF","TAGLN2","RNPS1","CDH1","KRT18",   
            "ZNF217","PERP","CYB5B","SFPQ","ANKRD10","EPCAM","WFDC2","EZR","TMEM125","LAPTM4B","SPINT1",  
            "HNRNPU","TFDP1","TKT","H1FX","CEACAM6","RALY","KRT19","CDKN2A","MUC21","NAPSA")
stem_marker<-region_marker[match(stemgene,region_marker$gene),]
stem_marker<-stem_marker[order(stem_marker$avg_log2FC,decreasing = T),]

spec_core<-c('MUC1','EZR','MUC21','TKT','CLDN4')
spec_bdy<-setdiff(bdy_marker,c(core_marker,dis_marker))[1:5]
spec_dis<-setdiff(dis_marker,c(spec_core,spec_bdy))[1:5]
#spec_mal<-intersect(core_marker,intersect(dis_marker,bdy_marker))[c(1,3:6)]
spec_imm<-c(setdiff(imm_marker,c(dis_marker,core_marker,bdy_marker))[1:4],'CD68')
#spec_nor<-setdiff(nor_marker,c(bdy_marker,core_marker,imm_marker,dis_marker))[1:5]

# cancer_marker<-c('SOD2','SCD5','SOX2','OLIG1','MLANA','MITF','DCT','MUC16','WFDC2','KIT')
# cancer_marker<-intersect(cancer_marker,rownames(st_count))
# region_marker[region_marker$gene%in%cancer_marker,]
# CoreBdyDis_marker[CoreBdyDis_marker$gene%in%cancer_marker,]
# MalImm_marker[MalImm_marker$gene%in%cancer_marker,]

select_marker<-c(spec_core,spec_bdy,spec_dis,spec_imm) %>% unique()
region_marker[region_marker$gene%in%select_marker,]

st_count<-st_rds@assays$Spatial@counts[select_marker,] %>% as.matrix()

core_site<-which(st_bdy$FinalLocalType=='Core')
bdy_site<-which(st_bdy$FinalLocalType=='Boundary')
dis_site<-which(st_bdy$FinalLocalType=='Dispersion')
imm_site<-which(st_bdy$FinalLocalType=='Immune')
#nor_site<-which(st_bdy$FinalLocalType=='Normal')

quantile_mean<-function(x){#x=cancer_data$Imm_score[which(cancer_data$LocalType=='Core')]
  quantile_a<-quantile(x)
  mean_robust<-(0.5*quantile_a[3]+0.25*(quantile_a[2]+quantile_a[4]))
}


percent_exp<-apply(st_count,1,function(x){#x=st_count[1,]
  Core<-length(which(x[core_site]>0))/length(core_site)
  Bdy<-length(which(x[bdy_site]>0))/length(bdy_site)
  Dis<-length(which(x[dis_site]>0))/length(dis_site)
  Imm<-length(which(x[imm_site]>0))/length(imm_site)
  #Nor<-length(which(x[nor_site]>0))/length(nor_site)
  return(c(Core,Bdy,Dis,Imm))
}) %>% t() %>% as.data.frame()
colnames(percent_exp)<-c('Core','Boundary','Dispersion','Immune')


scale_exp<-apply(st_count,1,function(x){#x=st_count[6,]
  Core<-mean(x[core_site])  ##range(x[core_site])  table(x[core_site])
  Bdy<-mean(x[bdy_site]) ##range(x[bdy_site])  table(x[bdy_site])
  Dis<-mean(x[dis_site])  ##range(x[dis_site])  table(x[dis_site])
  Imm<-mean(x[imm_site])   ##range(x[imm_site])  table(x[imm_site])
  #Nor<-mean(x[nor_site])
  aa_exp<-c(Core,Bdy,Dis,Imm)
  aa_exp<-aa_exp-min(aa_exp)
  aa_exp<-aa_exp/max(aa_exp)
  return(aa_exp)
}) %>% t() %>% as.data.frame()
colnames(scale_exp)<-c('Core','Boundary','Dispersion','Immune')


marker_exp<-reshape2::melt(as.matrix(percent_exp))
colnames(marker_exp)<-c('marker','localType','percent_exp')
marker_exp$scale_exp<-reshape2::melt(as.matrix(scale_exp))$value
marker_exp$color<-'black'


p_dot<-ggplot(marker_exp, aes(x=marker, y=localType,size=percent_exp,fill=scale_exp)) +
  geom_point(shape = 21,  color = marker_exp$color) + # 使用shape = 21画圈
  scale_fill_gradientn(colours = c(colorRampPalette(c("#DDDBDA","#F39C67"))(50),
                                   colorRampPalette(c("#F39C67","#B20A1C"))(50)) )+ #设置填充颜色
  scale_size_continuous(range = c(2, 8))+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=9),
        strip.text.y = element_text(size=9, face="bold"),
        axis.text.x=element_text(angle=90,vjust=1, hjust=1))
print(p_dot)
pdf('E:/Mirror/ST_analysis/pic/re/2/special_share_marker3.pdf',width = 8,height = 3)
print(p_dot)
dev.off()



st_localType<-st_bdy
st_localType<-st_localType[c(which(st_localType$FinalLocalType=='Core'),
                             which(st_localType$FinalLocalType=='Boundary'),
                             which(st_localType$FinalLocalType=='Dispersion'),
                             which(st_localType$FinalLocalType=='Immune'),
                             which(st_localType$FinalLocalType=='Normal')),]
select_marker<-c(spec_core,spec_bdy,spec_dis,spec_imm) %>% unique()

heat_data<-st_rds@assays$Spatial@counts[select_marker,st_localType$cell_name] %>% as.matrix()
table(st_count[7,bdy_site])
table(heat_data[7,which(st_localType$FinalLocalType=='Boundary')])


annotation_row = data.frame(
  geneType=factor(c(rep('Core_marker',length(spec_core)),rep('Bdy_marker',length(spec_bdy)),
                    rep('Dis_marker',length(spec_dis)),
                    rep('Imm_marker',length(spec_imm))))
)
rownames(annotation_row) = rownames(heat_data)
head(annotation_row)

annotation_col = data.frame(
  LocalType = st_localType$FinalLocalType
)
rownames(annotation_col) = colnames(heat_data)
head(annotation_col)

ann_colors = list(
  geneType=c("Core_marker"="#339999","Bdy_marker"="#99CC99","Dis_marker"="#CCFF99","Imm_marker"="#258ABF"),
  LocalType=c('Core'='#990033','Boundary'='#FF9900','Dispersion'='#CD5A5A','Immune'='#5477AF','Normal'='#669933')
)

#heat_data<-t(scale(t(heat_data)))
range(heat_data['CXCL9',])
range(heat_data['IGHG4',])
apply(heat_data,1,function(x) range(x))
heat_data<-apply(heat_data,1,function(x){
  x<-x-min(x)
  x<-x/max(x)
}) %>% t()

range(heat_data)
col_fun <- circlize::colorRamp2(
  seq(from=-2,to=2,length.out=5), 
  c("#0053A6","#7CB9DB", "white","#FCAD87","#9F001C")
)

p<-pheatmap(heat_data, 
            scale = "none",
            color=col_fun,
            #border_color = 'white',
            show_rownames = T,
            show_colnames = F,
            #labels_row = labels_row,
            cluster_rows = F,
            cluster_cols = F,
            fontsize = 10,
            annotation_colors = ann_colors,
            annotation_row = annotation_row,
            annotation_col = annotation_col,
            name = 'scale',
            use_raster=F
)
print(p)

pdf('E:/Mirror/ST_analysis/pic/re/2/select_marker_heatmap2.pdf',width = 9,height = 6)
print(p)
dev.off()



####各区域相关marker的柱状散点图绘制
###bdy
select_gene<-'CXCL9'
plot_data<-data.frame(gene=st_rds@assays[["Spatial"]]@counts[select_gene,],
                      localType=st_bdy$FinalLocalType)
spot_type<-c('Core',"Boundary","Dispersion",'Immune', 'Normal')
plot_data<-plot_data[plot_data$localType%in%spot_type,]
plot_data<-lapply(c('Core',"Boundary","Dispersion",'Immune', 'Normal'),function(x){#x='Core'
  xx<-plot_data$gene[which(plot_data$localType==x)]
  xx<-xx[which(xx<min(boxplot.stats(xx)$out))]
  xx_data<-data.frame(gene=xx,localType=x)
})
plot_data<-do.call(rbind,plot_data)
plot_data<-mutate(plot_data,Legend = factor(plot_data$localType, levels = spot_type))

median(plot_data$gene[which(plot_data$Legend=='Dispersion')])
median(plot_data$gene[which(plot_data$Legend=='Boundary')])
e <- ggplot(plot_data, aes(x = Legend, y = gene,fill=Legend))+ 
  geom_jitter(aes(color = Legend),size = 1,alpha=0.6,show.legend = F,
              position=position_jitterdodge(jitter.width = 0.8, 
                                            jitter.height = 0, 
                                            dodge.width = 1)) + # 不重叠的散点图
  stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8),show.legend = F) + # 误差棒，中位数，25%和75%分位数
  stat_summary(aes(fill = Legend), fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8),show.legend = F) + # 中位数水平线
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#FF9900","Core"="#990033","Dispersion"="#CD5A5A","Immune"="#5477AF","Normal"="#669933"))+
  scale_color_manual(values = c("Boundary"="#FF9900","Core"="#990033","Dispersion"="#CD5A5A","Immune"="#5477AF","Normal"="#669933"))+
  ggtitle(paste0('Bdy_',select_gene))+
  theme(plot.title = element_text(hjust = 0.4))+
  theme(plot.title = element_text(size = 12))+
  stat_compare_means(label.x = 1)+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(e)
?geom_jitter







