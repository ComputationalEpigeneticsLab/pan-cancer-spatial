####每种癌症类型取core的计算各个MP之间的相关性
library(corrplot)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(pheatmap)

color_pheatmap <- c(colorRampPalette(c("#0335FF","#809AFF"))(10),
                    colorRampPalette(c("#809AFF","white"))(10),
                    colorRampPalette(c("white","#DB878F"))(10),
                    colorRampPalette(c("#DB878F","#BE2331"))(10))

corr_re<-readRDS('/F3/F3c/corr_re.rds')
dis_dot<-as.matrix(corr_re[["p.adj"]])
dis_dot[as.matrix(corr_re[["p.adj"]]) >= 0.05 | abs(as.matrix(corr_re[["r"]]))<0.45]<-1
diag(dis_dot)<-1

dir_pic <- "/F3/F3c/"
me<-c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid")
for(j in 1:length(me)){#j=7
  pdf(paste0(dir_pic,'_','_',me[j],'_corr_new.pdf'),width = 8,height = 8)##selectSlice
  # par(mar = c(5, 4, 4, 2) + 0.1)  # 标准边距
  corrplot(corr_re$r, method = "circle", col = color_pheatmap, 
           order = "hclust",hclust.method = me[j],
           tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
           p.mat = dis_dot, ###一个显著性都没有时不能用
           diag = T, #type = 'upper',###显示一半及对角线
           addCoef.col = NULL,number.cex = 0.8,###相关性数字颜色及大小
           title = me[j],mar = c(0, 0, 2, 0),
           sig.level = c(0.0001,0.001), pch.cex = 1,
           insig = 'label_sig', pch.col = 'grey100')
  dev.off()
}



