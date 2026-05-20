####对每种癌症 使用所有切片core区marker的RRA排序结果和各种基因集计算GSEA
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
library(reshape2)


dir_out<-'/Fs3/'
GSEA_ESscore<-read.delim(paste0(dir_out,'GSEA_Stem&genesetScore_NESscore.txt'),
                         stringsAsFactors = F,check.names = F)
GSEA_FDR<-read.delim(paste0(dir_out,'GSEA_Stem&genesetScore_NESFDR.txt'),
                     stringsAsFactors = F,check.names = F)

GSEA_ESscore<-GSEA_ESscore[grep('Stem_of_',rownames(GSEA_ESscore)),]
GSEA_FDR<-GSEA_FDR[grep('Stem_of_',rownames(GSEA_FDR)),]

GSEA_ESscore<-GSEA_ESscore[setdiff(rownames(GSEA_ESscore),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_H3K27_bound')),]
GSEA_FDR<-GSEA_FDR[setdiff(rownames(GSEA_FDR),c('Stem_of_Ben_PRC2_targets','Stem_of_Ben_H3K27_bound')),]


dis_dot<-ifelse(as.matrix(GSEA_FDR) < 0.05, "*", "")
dis_dot[is.na(dis_dot)]<-""
range(GSEA_ESscore,na.rm = T)
stem_heat_data<-GSEA_ESscore %>% as.matrix()
range(stem_heat_data)
bk<-seq(-1,2,length.out=100)
color_pheatmap<-c(colorRampPalette(c("#0669AD",'#89BDD9'))(13),
                  colorRampPalette(c("#89BDD9",'white'))(20),
                  colorRampPalette(c("white",'#E9C1C6'))(20),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(46)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
column_order <- c('CRC','HNSCC','IPMN','CESC','CSCC','LIHC','DSRCT','EC','LUAD','LUSC','MIBC','GIST',
                  'HB','HGSC','PCNSL','PDAC','PN','PRAD','RCC','TC','OS','OSCC','OVCA','GBM',
                  'LGACC','LNC','HN-AS','NPC','GC','BRCA','DIPG','SKCM','LAM','PTCL','DLBCL')

# 确保数据框的列按照指定顺序排列
stem_heat_data <- stem_heat_data[, column_order]


p<-pheatmap::pheatmap(as.matrix(stem_heat_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = F,
                      cluster_cols = F,
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

