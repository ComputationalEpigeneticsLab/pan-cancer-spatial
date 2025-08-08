####所有干性基因集对癌症core区基因秩融合的GSEA结果
library(org.Hs.eg.db) #人类注释数据
# BiocManager::install('GO.db')
# BiocManager::install('clusterProfiler')
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(dplyr)
library(enrichplot)
library(fgsea)
library(stats)
library(Seurat)
library(stringr)
library(plyr)
#library(estimate)
library(RobustRankAggreg)
library(GseaVis)
library(stats)


###来源于125_各种基因集GSEA.R
GSEA_ESscore<-read.delim('E:/Mirror/ST_analysis/data/pic_data/GSEA_Stem&genesetScore_NESscore.txt',
                         stringsAsFactors = F,check.names = F)
GSEA_FDR<-read.delim('E:/Mirror/ST_analysis/data/pic_data/GSEA_Stem&genesetScore_NESFDR.txt',
                     stringsAsFactors = F,check.names = F)

GSEA_ESscore<-GSEA_ESscore[grep('Stem_of_',rownames(GSEA_ESscore)),]
GSEA_FDR<-GSEA_FDR[grep('Stem_of_',rownames(GSEA_FDR)),]

GSEA_ESscore<-GSEA_ESscore[setdiff(rownames(GSEA_ESscore),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_H3K27_bound')),]
GSEA_FDR<-GSEA_FDR[setdiff(rownames(GSEA_FDR),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_H3K27_bound')),]

# GSEA_ESscore<-GSEA_ESscore[setdiff(rownames(GSEA_ESscore),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_Eed_targets',
#                                                             'Stem_of_Ben_H3K27_bound','Stem_of_Ben_Suz12_targets',
#                                                             'Stem_of_Ben_NOS_TFs','Stem_of_Kim_et_al_core_m2h',
#                                                             'Stem_of_Kim_et_al_Myc_m2h')),]
# GSEA_FDR<-GSEA_FDR[setdiff(rownames(GSEA_FDR),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_Eed_targets',
#                                                 'Stem_of_Ben_H3K27_bound','Stem_of_Ben_Suz12_targets',
#                                                 'Stem_of_Ben_NOS_TFs','Stem_of_Kim_et_al_core_m2h',
#                                                 'Stem_of_Kim_et_al_Myc_m2h')),]

dis_dot<-ifelse(as.matrix(GSEA_FDR) < 0.05, "*", "")
dis_dot[is.na(dis_dot)]<-""
range(GSEA_ESscore,na.rm = T)
stem_heat_data<-GSEA_ESscore %>% as.matrix()
# stem_heat_data[which(stem_heat_data<0)]<-NA
# stem_heat_data<-scale(GSEA_ESscore)
# stem_heat_data<-apply(GSEA_ESscore,1,function(x){
#   x<-x-min(x)
#   x<-x/max(x)
# }) %>% t()
range(stem_heat_data)
bk<-seq(-1,2,length.out=100)
color_pheatmap<-c(colorRampPalette(c("#0669AD",'#89BDD9'))(13),
                  colorRampPalette(c("#89BDD9",'white'))(20),
                  colorRampPalette(c("white",'#E9C1C6'))(20),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(46)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(stem_heat_data), 
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
                      number_color = "black",
                      fontsize = 10,
                      cellwidth=15,
                      cellheight=15,
                      main = "Stem_GSEA",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/re/3/GSEA_stem/GSEA_stem_gai_NES.pdf',width = 7,height = 5)
print(p)
dev.off()

pdf('E:/Mirror/ST_analysis/pic/re/3/GSEA_stem/GSEA_stem_noScale_noFilter.pdf',width = 7,height = 5)
print(p)
dev.off()