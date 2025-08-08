#####ecosystem marker对免疫基因集做富集分析
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
library(stats)


geneset1<-read.delim('E:/Mirror/ST_analysis/data/geneset/genest_mp.txt',stringsAsFactors = F,check.names = F)
all_geneset1<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(all_geneset1)<-geneset1$pathway
pathway_data<-data.frame(pathway=rep(names(all_geneset1),unlist(lapply(all_geneset1,length))),
                         gene=unlist(all_geneset1))


ecosystem_cluster_maker<-read.delim('E:/Mirror/ST_analysis/other_data/ecosystem_cluster_maker.txt',stringsAsFactors = F,check.names = F)
ecosystem_cluster_maker<-ecosystem_cluster_maker[order(ecosystem_cluster_maker$scores,decreasing = T),]


GSEA_ESscore<-data.frame(pathway=names(all_geneset1))
GSEA_NESscore<-data.frame(pathway=names(all_geneset1))
GSEA_FDR<-data.frame(pathway=names(all_geneset1))
for(i in paste0('ecosystem',1:5)){#i='ecosystem1'
  marker1<-ecosystem_cluster_maker[which(ecosystem_cluster_maker$group==i),]
  marker1<-marker1[order(marker1$scores,decreasing = T),]
  geneList<-marker1$scores
  names(geneList)<-marker1$names
  
  GSEA_res <- GSEA(geneList,TERM2GENE = pathway_data,minGSSize = 2, maxGSSize = 1000, pvalueCutoff=1,eps = 0)
  
  slice_ESscore<-data.frame(pathway=GSEA_res@result$Description,
                            ESscore=GSEA_res@result$enrichmentScore)#enrichmentScore  NES
  colnames(slice_ESscore)[2]<-i
  slice_NESscore<-data.frame(pathway=GSEA_res@result$Description,
                             ESscore=GSEA_res@result$NES)#enrichmentScore  NES
  colnames(slice_NESscore)[2]<-i
  slice_FDR<-data.frame(pathway=GSEA_res@result$Description,
                        FDR=GSEA_res@result$p.adjust)
  colnames(slice_FDR)[2]<-i
  
  GSEA_ESscore<-merge(GSEA_ESscore,slice_ESscore,by='pathway',all=T)
  GSEA_NESscore<-merge(GSEA_NESscore,slice_NESscore,by='pathway',all=T)
  GSEA_FDR<-merge(GSEA_FDR,slice_FDR,by='pathway',all=T)
  print(i)
}
rownames(GSEA_ESscore)<-GSEA_ESscore$pathway
GSEA_ESscore<-GSEA_ESscore[,-1]
rownames(GSEA_NESscore)<-GSEA_NESscore$pathway
GSEA_NESscore<-GSEA_NESscore[,-1]
rownames(GSEA_FDR)<-GSEA_FDR$pathway
GSEA_FDR<-GSEA_FDR[,-1]

write.table(GSEA_ESscore,'E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_ESscore.txt',quote = F,sep = '\t')
write.table(GSEA_NESscore,'E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_NESscore.txt',quote = F,sep = '\t')
write.table(GSEA_FDR,'E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_FDR.txt',quote = F,sep = '\t')

GSEA_ESscore<-read.delim('E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_ESscore.txt')
GSEA_NESscore<-read.delim('E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_NESscore.txt')
GSEA_FDR<-read.delim('E:/Mirror/ST_analysis/data/pic_data/MP_eco_GSEA/Imm_GSEA_FDR.txt')

heat_data<-t(GSEA_NESscore)
range(heat_data,na.rm = T)
#heat_data<-t(scale(heat_data))###按列即按通路归一化
#heat_data<-t(scale(t(heat_data)))###按行即按MP类别归一化
range(heat_data)
dis_dot<-ifelse(as.matrix(t(GSEA_FDR)) < 0.05, "*", "")
dis_dot[is.na(dis_dot)]<-""

bk<-seq(-1,1,length.out=100)
color_pheatmap<-c(colorRampPalette(c("#0669AD",'#89BDD9'))(40),
                  colorRampPalette(c("#89BDD9",'white'))(10),
                  colorRampPalette(c("white",'#E9C1C6'))(10),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(40)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(heat_data),
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      #show_colnames = F,
                      cluster_rows = F,
                      cluster_cols = T,
                      treeheight_row = F,
                      # treeheight_col = T,
                      display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "black",
                      fontsize = 10,
                      # cellwidth=15,
                      # cellheight=15,
                      main = "ecosystem_ImmGSEA",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/NMF_module/ecosystem_ImmGSEA2.pdf',width = 8,height = 4.5)
print(p)
dev.off()
