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
library(reshape2)


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



####
library(enrichplot)
plotEnrichment(pathway = ref_down[[1]], stats = rank_FC, gseaParam = 1,ticksSize = 0.2)+ 
  labs(title="ATRA_down")

###干性基因集
geneset1<-read.delim('/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_StemGSEA/stem_geneset.txt',stringsAsFactors = F,check.names = F)
stemness_geneset<-lapply(1:nrow(geneset1),function(x) strsplit(geneset1[x,2],',')%>%unlist() )
names(stemness_geneset)<-paste0('Stem_of_',geneset1$geneset_name)
pathway_data<-data.frame(pathway=rep(names(stemness_geneset),unlist(lapply(stemness_geneset,length))),
                         gene=unlist(stemness_geneset))

# GSEA_ESscore<-read.delim('E:/Mirror/ST_analysis/data/pic_data/GSEA_Stem&genesetScore_ESscore.txt',
#                          stringsAsFactors = F,check.names = F)
# GSEA_FDR<-read.delim('E:/Mirror/ST_analysis/data/pic_data/GSEA_Stem&genesetScore_FDR.txt',
#                      stringsAsFactors = F,check.names = F)
# GSEA_ESscore<-GSEA_ESscore[names(stemness_geneset),]
# GSEA_FDR<-GSEA_FDR[names(stemness_geneset),]


####FC数据
dir_rds<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st/'
file_rds<-list.files(pattern = 'rds',path = dir_rds,recursive = T)
dataset_slice<-gsub('processed_|.rds','',file_rds)
cancer<-unlist(lapply(strsplit(dataset_slice,'/'),function(x)x[1])) %>% unique()
dataSet<-unlist(lapply(strsplit(dataset_slice,'/'),function(x)x[2]))
slice<-unlist(lapply(strsplit(dataset_slice,'/'),function(x)x[3]))
cancer[1:10]
dataSet[1:10]
slice[1:10]
dataset_slice[1:10]



dir_copykat<-'/data/zhouweiwei/ST_analysis/data/10X_Visium/new_st_copykat/'
file_bdy<-list.files(pattern = '_BdyCoreBud.txt',path = dir_copykat,recursive = T)
# file_bdy<-file_bdy[-grep('slice',file_bdy)]
length(file_bdy)
file_bdy[1:10]
table(dataset_slice==unlist(lapply(strsplit(file_bdy,'_BdyC'),function(x)x[1])))


for(i in 1:length(cancer)){
  #i=21 
  ###LUAD的
  cancer_slice_rds<-file_rds[grep(cancer[i],file_rds)]
  cancer_slice_bdy<-file_bdy[grep(cancer[i],file_bdy)]
  
  FC_list<-list()
  for(j in 1:length(cancer_slice_bdy)){
    st_rds<-readRDS(paste0(dir_rds,cancer_slice_rds[j]))
    st_count<-as.matrix(st_rds@assays[["Spatial"]]@counts)
    st_bdy<-read.delim(paste0(dir_bdy,cancer_slice_bdy[j]),stringsAsFactors = F,check.names = F)
    st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
    st_count<-st_count[,st_bdy$cell_name]
    core_site<-which(st_bdy$FinalLocalType=='Core')
    other_site<-which(st_bdy$FinalLocalType!='Core')
    FC_slice<-apply(st_count,1,function(x){#x=st_count[1,]
      aa<-mean(x[core_site])/mean(x[other_site])
    })
    FC_slice<-FC_slice[order(FC_slice,decreasing = T)]
    FC_list[[j]]<-names(FC_slice)
    print(cancer_slice_bdy[j])
  }
  
  
  ###秩融合
  RRA_score=aggregateRanks(FC_list)
  RRA_score<-RRA_score[order(RRA_score$Score,decreasing = F),]
  geneList<-1:nrow(RRA_score) %>% rev()
  names(geneList)<-RRA_score$Name
  
  GSEA_res <- GSEA(geneList,TERM2GENE = pathway_data,minGSSize = 2, maxGSSize = 1000, pvalueCutoff=1)
  
  library(GseaVis)
  #nn=which(GSEA_res@result[["Description"]]=='Stem_of_Zhang_Stem_Sig')
  nn=match(setdiff(rownames(GSEA_ESscore),'Stem_of_Kim_et_al_core_m2h'),GSEA_res@result[["Description"]])
  if(length(nn)==1){
    #?gseaNb
    for(nn in 1:length(GSEA_res@result[["Description"]])){#nn=1
    a_gsea<-gseaNb(object = GSEA_res,
                   geneSetID = GSEA_res@result[["Description"]][nn],
                   #curveCol = rep('#EB4747',18),
                   curveCol = c('#66CC00','#A78E41',"#CCCC99","#93C647","#CC3333","#ED703F","#D2AF83","#F3A383",'#7B3257',
                                "#8B964F","#FF9900","#EFA7A9","#EDDC6D","#FFFF00",'#CC9933','#FFCCCC','#996699','#A6B864'),
                   subPlot = 2,###只保留曲线
                   termWidth = 30,# 名称太长截断换行
                   addPval = T,
                   pvalX = 0.7,pvalY = 0.8,
                   pCol = 'black',
                   pHjust = 0
    )
    #pdf(paste0(dir_pic_GSEA,'EnrichPlot_',toupper(cancer[i]),'_',GSEA_res@result[["Description"]][nn],'.pdf'),height = 6,width = 8)
    pdf('/data/zhouweiwei/ST_analysis/data/10X_Visium/plot/GSEA_stem/LUAD_all_stemGeneSet2.pdf',height = 6,width = 8)
    print(a_gsea)
    dev.off()
    
    #}
  }
  
  
  print(cancer[i])
}

a_gsea<-gseaNb(object = GSEA_res,
                   geneSetID = GSEA_res@result[["Description"]][1:18],
                   #curveCol = rep('#EB4747',18),
                   curveCol = c('#66CC00','#A78E41',"#CCCC99","#93C647","#CC3333","#ED703F","#D2AF83","#F3A383",'#7B3257',
                                "#8B964F","#FF9900","#EFA7A9","#EDDC6D","#FFFF00",'#CC9933','#FFCCCC','#996699','#A6B864'),
                   subPlot = 2,###只保留曲线
                   #termWidth = 30,# 名称太长截断换行
                   #addPval = T,
                   #pvalX = 0.7,pvalY = 0.8,
                   #pCol = 'black',
                   #pHjust = 0
    )
    #pdf(paste0(dir_pic_GSEA,'EnrichPlot_',toupper(cancer[i]),'_',GSEA_res@result[["Description"]][nn],'.pdf'),height = 6,width = 8)
pdf('/data/zhouweiwei/ST_analysis/data/10X_Visium/plot/GSEA_stem/LUAD_all_stemGeneSet5.pdf',height = 12,width = 10)
print(a_gsea)
dev.off()  

p1<-gseaNb(object = gsea.re1,
           geneSetID = rownames(g1)[1:num1])

for(nn in 1:18){#nn=18
  a_GSEA<-plotEnrichment(pathway = stemness_geneset[[nn]], stats = geneList, gseaParam = 1,ticksSize = 0.2)+ 
    labs(title=paste0(toupper(cancer[i]),'_',names(stemness_geneset)[nn],'_ES:',round(GSEA_ESscore[nn,i],3),
                      '_FDR:',round(GSEA_FDR[nn,i],3)))
  pdf(paste0(dir_pic_GSEA,toupper(cancer[i]),'_',names(stemness_geneset)[nn],'_Enrichment','.pdf'),height = 6,width = 8)
  print(a_GSEA)
  dev.off()
  #?plotEnrichment
}

p <- gseaNb(
  object = GSEA_res,
  geneSetID = 1:5,  # 尝试多个ID
  curveCol = c('#66CC00','#A78E41',"#CCCC99","#93C647","#CC3333"),
  subPlot = 2,
  termWidth = 30,
  addPval = FALSE,  # 多个基因集时不显示p值
  pvalX = 0.7,
  pvalY = 0.8,
  pCol = 'black',
  pHjust = 0
)

pdf('/data/zhouweiwei/ST_analysis/data/10X_Visium/plot/GSEA_stem/LUAD_multi_gseaNb.pdf', 
    height = 6, width = 8)
print(p)
dev.off()



# 最简单有效的方法：使用gseaplot2绘制多个基因集
library(enrichplot)

show_n <- min(18, length(GSEA_res@result[["Description"]]))

# 准备颜色
my_colors <- c('#66CC00','#A78E41',"#CCCC99","#93C647","#CC3333","#ED703F","#D2AF83","#F3A383",'#7B3257',
               "#8B964F","#FF9900","#EFA7A9","#EDDC6D","#FFFF00",'#CC9933','#FFCCCC','#996699','#A6B864')

# 绘制多个基因集在同一曲线图上
p <- gseaplot2(
  x = GSEA_res,
  geneSetID = 1:show_n,
  color = my_colors[1:show_n],
  pvalue_table = FALSE,
  title = paste("Top", show_n, "Stem Gene Sets in LUAD"),
  ES_geom = "line",  # 使用线图
  subplots = 1:3     # 显示所有三个子图：ES曲线、排名、热图
)

# 调整图例
p <- p + theme(
  legend.position = "right",
  legend.text = element_text(size = 9),
  legend.title = element_blank()
)

# 保存
pdf('/data/zhouweiwei/ST_analysis/data/10X_Visium/plot/GSEA_stem/LUAD_all_stemGeneSets_combined1.pdf', 
    height = 8, width = 12)
print(p)
dev.off()

cat("已生成包含", show_n, "个基因集的组合GSEA图\n")



##########################挑LUAD_Stem_of_Ben_ES_exp1可视化富集曲线图
####使用其中靠前的基因进行热图展示

#i=11 ###LUAD多个切片基因秩融合结果
cancer_slice_rds<-file_rds[grep(cancer[i],file_rds)]
cancer_slice_bdy<-file_bdy[grep(cancer[i],file_bdy)]

FC_list<-list()
for(j in 1:length(cancer_slice_bdy)){
  st_rds<-readRDS(paste0(dir_rds,cancer_slice_rds[j]))
  st_count<-as.matrix(st_rds@assays[["Spatial"]]@counts)
  st_bdy<-read.delim(paste0(dir_bdy,cancer_slice_bdy[j]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType!='not.defined'),]
  st_count<-st_count[,st_bdy$cell_name]
  core_site<-which(st_bdy$FinalLocalType=='Core')
  other_site<-which(st_bdy$FinalLocalType!='Core')
  FC_slice<-apply(st_count,1,function(x){#x=st_count[1,]
    aa<-mean(x[core_site])/mean(x[other_site])
  })
  FC_slice<-FC_slice[order(FC_slice,decreasing = T)]
  FC_list[[j]]<-names(FC_slice)
  print(cancer_slice_bdy[j])
}


###秩融合
RRA_score=aggregateRanks(FC_list)
RRA_score<-RRA_score[order(RRA_score$Score,decreasing = F),]
geneList<-1:nrow(RRA_score) %>% rev()
names(geneList)<-RRA_score$Name

select_gene<-data.frame(RRA_gene=geneList)
select_gene$stem_gene<-names(geneList)
select_gene$stem_gene[!select_gene$stem_gene%in%stemness_geneset$Stem_of_Ben_ES_exp1]<-''


###找在所有切片里面都有的基因
all_slice_gene<-intersect(Reduce(intersect, FC_list),stemness_geneset$Stem_of_Ben_ES_exp1)
select_gene_gsea<-intersect(names(geneList)[1:740],all_slice_gene) ###45个基因  1800  740



cancer_slice_rds<-file_rds[grep(cancer[i],file_rds)]
cancer_slice_bdy<-file_bdy[grep(cancer[i],file_bdy)]
dataSlice<-unlist(lapply(strsplit(cancer_slice_bdy,'_'),function(x) x[1]))

core_mean<-c()
for(j in 1:length(cancer_slice_bdy)){#j=1
  st_rds<-readRDS(paste0(dir_rds,cancer_slice_rds[j]))
  st_bdy<-read.delim(paste0(dir_bdy,cancer_slice_bdy[j]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[which(st_bdy$FinalLocalType=='Core'),]
  
  st_count<-st_rds@assays[["Spatial"]]@counts[select_gene_gsea,st_bdy$cell_name] %>% as.matrix() %>% as.data.frame()
  slice_mean<-apply(st_count,1,mean)
  core_mean<-rbind(core_mean,slice_mean)
  rownames(core_mean)[j]<-dataSlice[j]
  
  print(cancer_slice_bdy[j])
}
write.table(core_mean,'/data/zhouweiwei/ST_analysis/data/10X_Visium/plot/GSEA_stem/selectGene_of_Ben_ES_exp1.txt',quote = F,sep = "\t")


select_gene_data<-read.table(file = "D:\\pan_cancer\\0_修稿\\F2\\GSEA\\selectGene_of_Ben_ES_exp1.txt",header = T,sep = "\t")
select_gene_data<-scale(select_gene_data)
select_gene_data<-apply(select_gene_data,1,function(x){
  x<-x-min(x)
  x<-x/max(x)
}) %>% t()
range(select_gene_data)
bk<-seq(0,1,length.out=100)
color_pheatmap<-c(colorRampPalette(c("white",'#E9C1C6'))(5),
                  colorRampPalette(c("#E9C1C6",'#BF404D'))(95)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(select_gene_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = T,
                      cluster_cols = F,
                      #treeheight_row = T,treeheight_col = T,
                      fontsize_number=15,
                      fontsize = 10,
                      # cellwidth=15,
                      # cellheight=15,
                      main = "selectGene_of_Ben_ES_exp1_LUAD",
                      breaks = bk,
                      name = 'scale_exp'
)
#print(p)
pdf('D:\\pan_cancer\\0_修稿\\F2\\GSEA\\selectGene_of_Ben_ES_exp1_2_scale.pdf',width = 6,height = 9)
print(p)
dev.off()


































