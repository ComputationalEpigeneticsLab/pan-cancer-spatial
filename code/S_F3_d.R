library("GSVA")
library(fgsea)
library(ggplot2)
library(ggpubr)

stem_geneset<-read.delim('E:/Mirror/ST_analysis/data/geneset/stem_geneset.txt',stringsAsFactors = F,check.names = F)
pathway_list_by <- by(stem_geneset,stem_geneset$geneset_name,function(x){return(unlist(strsplit(as.character(x[,2]),",")))})
pathway_list <- list()
for (i in 1:length(pathway_list_by)){
  pathway_list[[i]]=pathway_list_by[[i]]
  names(pathway_list)[i]=names(pathway_list_by[i])
}


dir_GSVAout<-'E:/Mirror/ST_analysis/data/TCGA/stem_GSVA/'
dir_TCGA<-'E:/Mirror/ST_analysis/data/TCGA/count/'
file_TCGA_count<-list.files(pattern = 'htseq_count.txt',path = dir_TCGA,recursive = T)
TCGA_cancer<-unlist(lapply(strsplit(file_TCGA_count,'/'),function(x)x[1]))

###mRNAsi和StemnessIndex
i=10 ###LUAD
TCGA_count<-read.delim(paste0(dir_TCGA,file_TCGA_count[i]),stringsAsFactors = F,check.names = F)
mRNAsi_model<-readRDS("E:/Mirror/cell_stem/coda/innerpath/model.RNA_test.Rds")
common <- intersect(names(mRNAsi_model$w), rownames(TCGA_count))
X <- TCGA_count[common, ]
w <- mRNAsi_model$w[common]
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
mRNAsi_score<-data.frame(sample=colnames(TCGA_count),score=score)


####LUAD的TCGA样本不同阶段的干性基因集GSVA富集程度
LUAD_stem<-read.delim('E:/Mirror/ST_analysis/data/TCGA/stem_GSVA/LUAD/stem_GSVA.txt',stringsAsFactors = F,check.names = F)
LUAD_stem<-t(LUAD_stem) %>% as.data.frame()
LUAD_stem$mRNAsi_score<-mRNAsi_score[rownames(LUAD_stem),'score']
LUAD_stem$StemnessIndex_score<-StemnessIndex_score[rownames(LUAD_stem),'score']
LUAD_stem<-t(LUAD_stem) %>% as.data.frame()
write.table(LUAD_stem,'E:/Mirror/ST_analysis/data/TCGA/stem_GSVA/LUAD/stem_GSVA.txt',quote = F,sep = '\t')

LUAD_pheno<-read.delim('E:/Mirror/ST_analysis/data/TCGA/count/LUAD/TCGA-LUAD.GDC_phenotype.tsv',stringsAsFactors = F,check.names = F)

length(intersect(colnames(LUAD_stem),LUAD_pheno$submitter_id.samples))
LUAD_stage<-LUAD_pheno[,c('submitter_id.samples','tumor_stage.diagnoses')]
rownames(LUAD_stage)<-LUAD_stage$submitter_id.samples
LUAD_stage<-LUAD_stage[colnames(LUAD_stem),]
table(LUAD_stage$submitter_id.samples==colnames(LUAD_stem))
LUAD_stage$stage<-LUAD_stage$tumor_stage.diagnoses
table(LUAD_stage$tumor_stage.diagnoses)
LUAD_stage$stage[grep('iv',LUAD_stage$stage)]<-'IV'
LUAD_stage$stage[grep('iii',LUAD_stage$stage)]<-'III'
LUAD_stage$stage[grep('ii',LUAD_stage$stage)]<-'II'
LUAD_stage$stage[grep('i',LUAD_stage$stage)]<-'I'
table(LUAD_stage$stage,LUAD_stage$tumor_stage.diagnoses)
LUAD_stage<-LUAD_stage[which(LUAD_stage$stage!='not reported'),]
LUAD_stage<-LUAD_stage[order(LUAD_stage$stage),]

table(unlist(lapply(strsplit(LUAD_stage$submitter_id.samples,'-'),function(x)x[4])))
LUAD_stage<-LUAD_stage[grep('01A',LUAD_stage$submitter_id.samples),]
table(LUAD_stage$stage)

LUAD_stem<-LUAD_stem[,rownames(LUAD_stage)]


####干性得分归一化到01之间
LUAD_stem<-apply(LUAD_stem,1,function(x){
  x<-x-min(x)
  x<-x/max(x)
}) %>% t()


####IV与I II III的FC值与秩和检验
IV_site<-which(LUAD_stage$stage=='IV')
ItoIII_site<-which(LUAD_stage$stage!='IV')

IIIandIV_site<-which(LUAD_stage$stage=='IV'|LUAD_stage$stage=='III')
IandII_site<-which(LUAD_stage$stage=='I'|LUAD_stage$stage=='II')

I_site<-which(LUAD_stage$stage=='I')
II_site<-which(LUAD_stage$stage=='II')
III_site<-which(LUAD_stage$stage=='III')
IV_site<-which(LUAD_stage$stage=='IV')

LUAD_stem_FC<-apply(LUAD_stem,1,function(x){
  #mean(x[IV_site])/mean(x[ItoIII_site])
  mean(x[c(II_site,III_site,IV_site)])/mean(x[I_site])
  #mean(x[IV_site])/mean(x[c(II_site,III_site,I_site)])
})

LUAD_stem_P<-apply(LUAD_stem,1,function(x){
  #wilcox.test(x[IV_site],x[ItoIII_site],alternative = 'greater')[["p.value"]]
  wilcox.test(x[c(II_site,III_site,IV_site)],x[I_site],alternative = 'greater')[["p.value"]]
  #wilcox.test(x[IV_site],x[c(II_site,III_site,I_site)],alternative = 'greater')[["p.value"]]
})



stem_name<-rownames(LUAD_stem)

for(i in 1:length(stem_name)){#i=21
  plot_data<-data.frame(stage=LUAD_stage$stage,score=as.vector(as.matrix(LUAD_stem[stem_name[i],])))
  #plot_data$stage[which(plot_data$stage=='II'|plot_data$stage=='III'|plot_data$stage=='I')]<-'ItoIII'
  plot_data$stage[which(plot_data$stage=='II'|plot_data$stage=='III'|plot_data$stage=='IV')]<-'IItoIV'
  
  
  plot_data<-lapply(unique(plot_data$stage),function(x){#x='II'
    xx<-plot_data$score[which(plot_data$stage==x)]
    xx<-setdiff(xx,boxplot.stats(xx)$out)
    xx_data<-data.frame(stage=x,score=xx)
    return(xx_data)
  })
  plot_data<-do.call(rbind,plot_data)
  # plot_data$stage<-"I"
  # plot_data$type<-c(rep('a',100),rep('b',396))
  
  e <- ggplot(plot_data, aes(x = stage, y = score,fill=stage))+ 
    # stat_boxplot(geom = 'errorbar',width=0.5,position = position_dodge(0.9))+
    # geom_boxplot(aes(fill = stage), color='black',width = 0.8,#lwd=0.3,fatten=0.9,
    #              position = position_dodge(0.9),alpha=1,outlier.alpha=0)+
    geom_jitter(aes(color = stage),size = 1,alpha=1,show.legend = F,
                position=position_jitterdodge(jitter.width = 0.8,
                                              jitter.height = 0,
                                              dodge.width = 1)) + # 不重叠的散点图
    stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8),show.legend = F) + # 误差棒，中位数，25%和75%分位数
    stat_summary(aes(fill = stage), fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8),show.legend = F) + # 中位数水平线
    theme_classic(base_size = 20)+
    theme(axis.text = element_text(color = 'black'))+
    scale_fill_manual(values = c('ItoIII'='#B8DCC2','I'='#B8DCC2','II'='#2A95B7','III'='#2A95B7','IV'='#295AB7','IItoIV'='#295AB7'))+
    scale_color_manual(values = c('ItoIII'='#B8DCC2','I'='#B8DCC2','II'='#2A95B7','III'='#2A95B7','IV'='#295AB7','IItoIV'='#295AB7'))+
    ggtitle(stem_name[i])+
    theme(plot.title = element_text(hjust = 0.4))+
    theme(plot.title = element_text(size = 12))+
    stat_compare_means(label.x = 1)+
    #theme_bw()+
    theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
          axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
  pdf(paste0('E:/Mirror/ST_analysis/pic/re/3/TCGA_stem/Box_IItoIVvsI_',stem_name[i],'.pdf'),width = 4,height = 3.5)
  print(e)
  dev.off()
  
}













