####所有TCGA样本做生存分析
library(dplyr)
library(ggrepel)
library(psych)


gene_anno<-read.delim('E:/Mirror/AML_DNAme/data/TCGA/allCancer/gencode.v36.annotation.gtf.gene.probemap')

#######所有TCGA癌症 count表达谱和生存数据 只要01A 11A样本 且生存天数大于等于30天
dir_out<-'E:/Mirror/ST_analysis/data/TCGA/allCancer/'
dir_TCGA<-'E:/Mirror/AML_DNAme/data/TCGA/allCancer/count/'
file_TCGA_count<-list.files(pattern = 'star_counts.tsv.gz',path = dir_TCGA,recursive = T)
dir_sur<-'E:/Mirror/AML_DNAme/data/TCGA/allCancer/survival/'
file_TCGA_pheno<-list.files(pattern = 'survival.tsv.gz',path = dir_sur,recursive = T)
TCGA_cancer<-unlist(lapply(strsplit(file_TCGA_count,'.star_counts'),function(x)x[1]))
table(TCGA_cancer==unlist(lapply(strsplit(TCGA_cancer,'.survival'),function(x)x[1])))

all_TCGA_count<-list()
all_TCGA_pheno<-c()
for(i in 1:length(file_TCGA_count)){
  #i=1
  #if(!dir.exists(paste0(dir_RRAout,TCGA_cancer[i]))) dir.create(paste0(dir_GSVAout,TCGA_cancer[i]))
  TCGA_count<-read.delim(paste0(dir_TCGA,file_TCGA_count[i]),stringsAsFactors = F,check.names = F)
  expr_matrix2<-merge(gene_anno,TCGA_count,by.x='id',by.y='Ensembl_ID',all=T)
  table(is.na(expr_matrix2$gene))
  expr_matrix2<-expr_matrix2[!is.na(expr_matrix2$gene),]
  expr_matrix2<-expr_matrix2[!duplicated(expr_matrix2$gene),]
  rownames(expr_matrix2)<-expr_matrix2$gene
  expr_matrix2<-expr_matrix2[,-c(1:6)]
  expr_matrix3<-round((2^expr_matrix2)-1)
  
  table(unlist(lapply(strsplit(colnames(expr_matrix3),'-'),function(x) x[4])))
  #expr_matrix3<-expr_matrix3[,grep('01A|01A|03A|04A|05A|06A|07A|08A|09A',colnames(expr_matrix3))]
  colnames(expr_matrix3)<-paste0(TCGA_cancer[i],'_',colnames(expr_matrix3))
  #TCGA_count[1:5,1:5]
  
  TCGA_pheno<-read.delim(paste0(dir_sur,file_TCGA_pheno[i]),stringsAsFactors = F,check.names = F)
  TCGA_pheno<-TCGA_pheno[,1:3]
  rownames(TCGA_pheno)<-paste0(TCGA_cancer[i],'_',TCGA_pheno$sample)
  #TCGA_pheno<-TCGA_pheno[which(TCGA_pheno$OS.time>=30),]
  
  expr_matrix3<-expr_matrix3[,intersect(colnames(expr_matrix3),rownames(TCGA_pheno))]
  TCGA_pheno<-TCGA_pheno[intersect(colnames(expr_matrix3),rownames(TCGA_pheno)),]
  TCGA_pheno$cancer<-TCGA_cancer[i]
  TCGA_pheno$sample_type<-unlist(lapply(strsplit(TCGA_pheno$sample,'-'),function(x) x[4]))
  
  all_TCGA_count[[i]]<-expr_matrix3
  names(all_TCGA_count)[i]<-TCGA_cancer[i]
  all_TCGA_pheno<-rbind(all_TCGA_pheno,TCGA_pheno)
  print(TCGA_cancer[i])
  print(table(unlist(lapply(strsplit(colnames(expr_matrix3),'-'),function(x) x[5]))))
}

table(rownames(all_TCGA_count[[1]])==rownames(all_TCGA_count[[2]]))
lapply(all_TCGA_count,function(x){
  table(rownames(all_TCGA_count[[1]])==rownames(x))
}) %>% unlist()
all_TCGA_count<-do.call(cbind,all_TCGA_count)
colnames(all_TCGA_count)<-unlist(lapply(strsplit(colnames(all_TCGA_count),'\\.'),function(x)x[2]))
table(colnames(all_TCGA_count)==rownames(all_TCGA_pheno))

sparsematrix <- as(all_TCGA_count, "sparseMatrix")
sparsematrix[1:5,1:5]
all_TCGA_count[1:5,1:5]
saveRDS(sparsematrix,file = 'E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_count.rds')
write.table(all_TCGA_pheno,'E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_pheno.txt',quote = F,sep = '\t')


all_TCGA_count<-as.matrix(readRDS('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_count.rds'))
all_TCGA_pheno<-read.delim('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_pheno.txt',stringsAsFactors = F,check.names = F)

all_TCGA_count<-as.matrix(readRDS('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/all_TCGA_count.rds'))
all_TCGA_pheno<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/all_TCGA_pheno.txt',stringsAsFactors = F,check.names = F)
###去批次
library(sva)
library(tinyarray)

table(all_TCGA_pheno$cancer)
table(all_TCGA_pheno$sample_type)
all_TCGA_pheno<-all_TCGA_pheno[grep('A',all_TCGA_pheno$sample_type),]
all_TCGA_count<-all_TCGA_count[,rownames(all_TCGA_pheno)]
all_TCGA_pheno$T_N<-'tumor'
all_TCGA_pheno$T_N[which(all_TCGA_pheno$sample_type=='11A')]<-'normal'
table(all_TCGA_pheno$sample_type,all_TCGA_pheno$T_N)

pdf('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/before_T_N.pdf',height = 8,width=8)
draw_pca(exp = all_TCGA_count, group_list = factor(all_TCGA_pheno$T_N))
dev.off()

pdf('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/before_cancer.pdf',height = 8,width=8)
draw_pca(exp = all_TCGA_count, group_list = factor(all_TCGA_pheno$cancer))
dev.off()


expr_count_combat <- ComBat_seq(counts = as.matrix(all_TCGA_count), 
                                batch = all_TCGA_pheno$cancer,
                                group = all_TCGA_pheno$T_N
)


pdf('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/after_T_N.pdf',height = 8,width=8)
draw_pca(exp = expr_count_combat, group_list = factor(all_TCGA_pheno$T_N))
dev.off()

pdf('/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/after_cancer.pdf',height = 8,width=8)
draw_pca(exp = expr_count_combat, group_list = factor(all_TCGA_pheno$cancer))
dev.off()

sparsematrix_combat <- as(expr_count_combat, "sparseMatrix")
sparsematrix_combat[1:5,1:5]
expr_count_combat[1:5,1:5]
saveRDS(sparsematrix_combat,file = '/data/zhouweiwei/ST_analysis/data/TCGA/allCancer/all_TCGA_count_combat.rds')



####重新生存分析
###得到的模型对所有TCGA癌症验证
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(dplyr)
library(stats)
library(Seurat)
library(stringr)
library(survminer)
#remove.packages("survival")
library(survival)
library(combinat)
library(meta) #加载包
library(stringr)
library(grid)


###使用空转对应的TCGA癌症样本进行生存模型构建
all_TCGA_count<-readRDS('E:/Mirror/ST_analysis/data/TCGA/ST_RRA/all_TCGA_count.rds')
all_TCGA_pheno<-read.delim('E:/Mirror/ST_analysis/data/TCGA/ST_RRA/all_TCGA_pheno.txt',stringsAsFactors = F,check.names = F)
all_TCGA_pheno$cancer<-unlist(lapply(strsplit(rownames(all_TCGA_pheno),'_TCGA'),function(x)x[1]))
TCGA_cancer<-unique(unlist(lapply(strsplit(rownames(all_TCGA_pheno),'_TCGA'),function(x) x[1])))###没有PCNSL HN-AS
table(unlist(lapply(strsplit(all_TCGA_pheno$submitter_id.samples,'-'),function(x) x[4])))###没有PCNSL HN-AS

allCancer_TCGA_count<-readRDS('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_count.rds')
allCancer_TCGA_pheno<-read.delim('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_pheno.txt',stringsAsFactors = F,check.names = F)
allCancer_TCGA_pheno$cancer<-unlist(lapply(strsplit(allCancer_TCGA_pheno$cancer,'-'),function(x)x[2]))
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[grep('A',allCancer_TCGA_pheno$sample_type),]
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[which(allCancer_TCGA_pheno$OS.time>=30),]
table(allCancer_TCGA_pheno$sample_type)
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[which(allCancer_TCGA_pheno$sample_type!='11A'),]
unique(allCancer_TCGA_pheno$cancer)
table(allCancer_TCGA_pheno$cancer)

allCancer_TCGA_count<-allCancer_TCGA_count[,rownames(allCancer_TCGA_pheno)]

match_cancer<-c('brca','cesc','gbm','lihc','luad','ovca','prad','skcm') %>% toupper()###名称完全对应的癌症
pan_cancer<-c("BRCA","CESC","CRC","CSCC","GBM","GIST","LIHC","LUAD","MIBC","OVCA","PDAC","PRAD","RCC","SKCM")
PartMatch_cancer<-c("BRCA","CESC","COAD","HNSC","GBM","STAD","LIHC","LUAD","BLCA","OV","PAAD","PRAD","KIRC","SKCM")

pan_data_pheno<-allCancer_TCGA_pheno[allCancer_TCGA_pheno$cancer%in%PartMatch_cancer,]
pan_data_count<-allCancer_TCGA_count[,rownames(pan_data_pheno)]
dim(pan_data_count)

####使用13个MP的基因 CAF TAM的marker cellchat的相关基因#####
MP13<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/MP_list_intersect_initia25_intersect_cluster24.txt',
                 stringsAsFactors = F,check.names = F)
# MP13<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/MP_list_intersect_initia25_intersect_cluster24.txt',
#                  stringsAsFactors = F,check.names = F)
MP13<-MP13[,-1]
MP13<-reshape2::melt(as.matrix(MP13)) %>% as.data.frame()
MP13<-MP13[,-1]
length(unique(MP13$value))


###CAF TAM marker基因####
CAF_marker<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/CAF_TAM_marker/CAF_marker.txt')
TAM_marker<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/CAF_TAM_marker/TAM_marker.txt')
# CAF_marker<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/CAF_marker.txt')
# TAM_marker<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/TAM_marker.txt')
###数据来自 272_CAFTAM主导marker.R

length(unique(CAF_marker$marker[which(CAF_marker$sum>0)]))
length(unique(TAM_marker$marker[which(TAM_marker$sum>0)]))

####cellchat基因######
bdy_step1CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF_diviLR.txt',stringsAsFactors=F,check.names=F)
bdy_step1TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM_diviLR.txt',stringsAsFactors=F,check.names=F)

# matrix_CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_diviLR.txt',stringsAsFactors = F,check.names = F)
# matrix_TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_diviLR.txt',stringsAsFactors = F,check.names = F)

# bdy_step1CAF_divi<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/bdy_step1CAF_diviLR.txt',stringsAsFactors=F,check.names=F)
# bdy_step1TAM_divi<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/bdy_step1TAM_diviLR.txt',stringsAsFactors=F,check.names=F)

# matrix_CAF_divi<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/matrix_CAF_diviLR.txt',stringsAsFactors = F,check.names = F)
# matrix_TAM_divi<-read.delim('/data/zhouweiwei/ST_analysis/data/TCGA/use_gene/matrix_TAM_diviLR.txt',stringsAsFactors = F,check.names = F)
#sum(matrix_CAF_divi$LR_num)

bdy_step1CAF_divi$L_R<-'source'
bdy_step1CAF_divi$L_R[which(bdy_step1CAF_divi$target=='CAF')]<-'target'
bdy_CAFpath<-as.data.frame.array(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$L_R))
bdy_CAFpath$sum<-bdy_CAFpath$source+bdy_CAFpath$target
bdy_CAF_data<-bdy_step1CAF_divi[bdy_step1CAF_divi$pathway_name%in%rownames(bdy_CAFpath)[order(bdy_CAFpath$sum,decreasing = T)][1:3],]
bdy_CAF_gene<-c(unlist(strsplit(bdy_CAF_data$ligand,'_')),unlist(strsplit(bdy_CAF_data$receptor,'_'))) %>% unique()

bdy_step1TAM_divi$L_R<-'source'
bdy_step1TAM_divi$L_R[which(bdy_step1TAM_divi$target=='TAM')]<-'target'
bdy_TAMpath<-as.data.frame.array(table(bdy_step1TAM_divi$pathway_name,bdy_step1TAM_divi$L_R))
bdy_TAMpath$sum<-bdy_TAMpath$source+bdy_TAMpath$target
bdy_TAM_data<-bdy_step1TAM_divi[bdy_step1TAM_divi$pathway_name%in%rownames(bdy_TAMpath)[order(bdy_TAMpath$sum,decreasing = T)][1:3],]
bdy_TAM_gene<-c(unlist(strsplit(bdy_TAM_data$ligand,'_')),unlist(strsplit(bdy_TAM_data$receptor,'_'))) %>% unique()

# matrix_CAF_divi$L_R<-'source'
# matrix_CAF_divi$L_R[which(matrix_CAF_divi$target=='CAF')]<-'target'
# matrix_CAFpath<-as.data.frame.array(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$L_R))
# matrix_CAFpath$sum<-matrix_CAFpath$source+matrix_CAFpath$target
# matrix_CAF_data<-matrix_CAF_divi[matrix_CAF_divi$pathway_name%in%rownames(matrix_CAFpath)[order(matrix_CAFpath$sum,decreasing = T)][1:3],]
# matrix_CAF_gene<-c(unlist(strsplit(matrix_CAF_data$ligand,'_')),unlist(strsplit(matrix_CAF_data$receptor,'_'))) %>% unique()
# 
# matrix_TAM_divi$L_R<-'source'
# matrix_TAM_divi$L_R[which(matrix_TAM_divi$target=='TAM')]<-'target'
# matrix_TAMpath<-as.data.frame.array(table(matrix_TAM_divi$pathway_name,matrix_TAM_divi$L_R))
# matrix_TAMpath$sum<-matrix_TAMpath$source+matrix_TAMpath$target
# matrix_TAM_data<-matrix_TAM_divi[matrix_TAM_divi$pathway_name%in%rownames(matrix_TAMpath)[order(matrix_TAMpath$sum,decreasing = T)][1:3],]
# matrix_TAM_gene<-c(unlist(strsplit(matrix_TAM_data$ligand,'_')),unlist(strsplit(matrix_TAM_data$receptor,'_|:'))) %>% unique()

# cellchatGene<-c(bdy_CAF_gene,bdy_TAM_gene,matrix_CAF_gene,matrix_TAM_gene) %>% unique()
cellchatGene<-c(bdy_CAF_gene,bdy_TAM_gene) %>% unique()
length(unique(cellchatGene))
# cellchatGene<-c('COL1A1','ITGB1','COL1A2','SDC4',
#                 'SPP1','ITGB1','ANGPTL4','ITGA5',
#                 'FN1','CD44','ICAM1','ITGB2',
#                 'THBS1','CD47','THBS1','SDC1')
# intersect(rownames(coef_re),MP13$value)

# lapply(rownames(coef_re),function(x){
#   c(x%in%MP13$value,x%in%CAF_marker$marker[which(CAF_marker$sum>0)],
#     x%in%TAM_marker$marker[which(TAM_marker$sum>0)])
# })
###合并基因####
use_gene<-unique(c(MP13$value,
                   CAF_marker$marker[which(CAF_marker$sum>0)],
                   TAM_marker$marker[which(TAM_marker$sum>0)],
                   cellchatGene))
length(intersect(use_gene,rownames(all_TCGA_count)))###与TCGA交集为468个基因

use_gene<-unique(MP13$value)

use_gene<-intersect(use_gene,rownames(all_TCGA_count))

use_gene1<-intersect(use_gene,rownames(all_TCGA_count))
sur_count<-all_TCGA_count[use_gene1,]


###lasso-cox###
library(glmnet)
library(meta) #加载包
library(stringr)

CoxExample<-list()
CoxExample$x<-t(as.matrix(sur_count))
CoxExample$y<-as.matrix(all_TCGA_pheno[,c('OS.time','OS')])
colnames(CoxExample$y) <- c('time', 'status')
head(CoxExample$x)[,1:5] # 行：患者/样本，列：特征，可以是特征基因的表达谱
head(CoxExample$y)

fit = glmnet(CoxExample$x,CoxExample$y,
             alpha = 1,family = "cox") 
#?glmnet
plot(fit,xvar= 'lambda',label = TRUE)
## 交叉验证
# 默认 nfolds = 10
set.seed(111)
cvfit <- cv.glmnet(CoxExample$x,CoxExample$y,
                   alpha = 1,family = "cox")
#?cv.glmnet
plot(cvfit)

###筛选特征
cvfit$lambda.min
cvfit$lambda.1se
coefficient <- coef(cvfit, s= cvfit$lambda.min)###29个基因用的是 lambda.min  28基因用的是 lambda.1se
selected_index <- which(as.numeric(coefficient) != 0)
selected_features <- names(coefficient[selected_index,])
coefficient<-as.matrix(coefficient)
coefficient<-coefficient[which(coefficient[,1]!=0),]


###多因素cox
dat<-cbind(CoxExample$y,CoxExample$x[,selected_features]) %>% as.data.frame()
colnames(dat)<-gsub('-','_',colnames(dat))
res.cox <- coxph(Surv(time, status) ~ ., data =  dat)
res.cox_sum<-summary(res.cox)
coef_re<-res.cox_sum[["coefficients"]] %>% as.data.frame()
coef_re<-coef_re[which(coef_re$`Pr(>|z|)`<0.01),]

coef_re$coef<-10000*coef_re$coef

write.table(coef_re,'E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/lasso_cox_coef_LR.txt',quote = F,sep = '\t')
coef_re<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/lasso_cox_coef_LR.txt',stringsAsFactors = F)
rownames(coef_re)<-gsub('_','-',rownames(coef_re))
#coef_re<-coef_re[which(coef_re$coef>0),]
coef_func<-function(x){
  score<-sum(x*coef_re$coef)
  return(score)
}


###33种TCGA癌症
allCancer_TCGA_count<-readRDS('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_count.rds')
allCancer_TCGA_pheno<-read.delim('E:/Mirror/ST_analysis/data/TCGA/allCancer/all_TCGA_pheno.txt',stringsAsFactors = F,check.names = F)
allCancer_TCGA_pheno$cancer<-unlist(lapply(strsplit(allCancer_TCGA_pheno$cancer,'-'),function(x)x[2]))
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[grep('A',allCancer_TCGA_pheno$sample_type),]
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[which(allCancer_TCGA_pheno$OS.time>=30),]
table(allCancer_TCGA_pheno$sample_type)
allCancer_TCGA_pheno<-allCancer_TCGA_pheno[which(allCancer_TCGA_pheno$sample_type!='11A'),]

allCancer_TCGA_count<-allCancer_TCGA_count[,rownames(allCancer_TCGA_pheno)]

length(intersect(rownames(coef_re),rownames(allCancer_TCGA_count)))
# pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_KM_LR.pdf',width = 4,height = 3)
sur_data<-data.frame(sample=allCancer_TCGA_pheno$sample,row.names = rownames(allCancer_TCGA_pheno),cancer=allCancer_TCGA_pheno$cancer,
                     days_to_death=allCancer_TCGA_pheno$OS.time,
                     status=allCancer_TCGA_pheno$OS,
                     score=apply(allCancer_TCGA_count[rownames(coef_re),],2,coef_func),
                     group='Not')
#sur_data$group[which(sur_data$score<median(sur_data$score))]<-'Low'

optimalCutoff=surv_cutpoint(sur_data, time = "days_to_death", event = "status",
                            variables = c("score"))
res.cut<- surv_categorize(optimalCutoff)
sur_data$group<-str_to_title(as.character(res.cut$score))
sur_data$group<-str_to_title(sur_data$group)
table(sur_data$group)
write.table(sur_data,'E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/all_sur_score_LR.txt',quote = F,sep = '\t')



if(length(which(sur_data$group=='High'))>10&length(which(sur_data$group=='Low'))>10){
  # fit <- survfit(Surv(days_to_death,status) ~ group,data = sur_data)
  # p1<-ggsurvplot(fit, data = sur_data,pval = TRUE,legend.title = 'pan_cancer',
  #                palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
  # print(p1)
  # pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/pan_cancer_sur.pdf',width = 8,height = 7)
  # print(p1)
  # dev.off()
  
  # fit <- survfit(Surv(days_to_death,status) ~ group,data = sur_data)
  # p1<-ggsurvplot(fit, data = sur_data,pval = TRUE,legend.title = 'pan_cancer',
  #                palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
  # print(p1)
  
  sdiff <- survdiff(Surv(days_to_death,status) ~ group,data = sur_data)
  p.val = 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  HR = (sdiff$obs[1]/sdiff$exp[1])/(sdiff$obs[2]/sdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1]))
}else{
  p.val<-1
  HR<-1
  up95 = 1
  low95 = 1
}

cancer<-unique(allCancer_TCGA_pheno$cancer)
for(j in 1:length(cancer)){#j=11
  sur_data_pheno<-allCancer_TCGA_pheno[allCancer_TCGA_pheno$cancer%in%cancer[j],]
  sur_data_count<-allCancer_TCGA_count[,rownames(sur_data_pheno)]
  
  sur_data<-data.frame(sample=sur_data_pheno$sample,row.names = rownames(sur_data_pheno),
                       days_to_death=sur_data_pheno$OS.time,
                       status=sur_data_pheno$OS,
                       score=apply(sur_data_count[rownames(coef_re),],2,coef_func),
                       group='Not')
  #sur_data$group[which(sur_data$score<median(sur_data$score))]<-'Low'
  optimalCutoff=surv_cutpoint(sur_data, time = "days_to_death", event = "status",
                              variables = c("score"))
  res.cut<- surv_categorize(optimalCutoff)
  sur_data$group<-str_to_title(as.character(res.cut$score))
  # sur_data$group<-str_to_title(sur_data$group)
  table(sur_data$group)
  
  # fit <- survfit(Surv(days_to_death,status) ~ group,data = sur_data)
  # p1<-ggsurvplot(fit, data = sur_data,pval = TRUE,legend.title = cancer[j],
  #                palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
  # print(p1)
  
  if(length(which(sur_data$group=='High'))>10&length(which(sur_data$group=='Low'))>10){
    sdiff <- survdiff(Surv(days_to_death,status) ~ group,data = sur_data)
    p.val = c(p.val,1 - pchisq(sdiff$chisq, length(sdiff$n) - 1))
    HR = c(HR,(sdiff$obs[1]/sdiff$exp[1])/(sdiff$obs[2]/sdiff$exp[2]))
    HR_i<-(sdiff$obs[1]/sdiff$exp[1])/(sdiff$obs[2]/sdiff$exp[2])
    up95 = c(up95,exp(log(HR_i) + qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1])))
    low95 = c(low95,exp(log(HR_i) - qnorm(0.975)*sqrt(1/sdiff$exp[2]+1/sdiff$exp[1])))
  }else{
    p.val = c(p.val,1)
    HR = c(HR,1)
    up95 = c(up95,1)
    low95 = c(low95,1)
  }
}
dev.off()

cancer<-unique(all_sur_score_LR$cancer)
PartMatch_cancer<-c("BRCA","CESC","COAD","HNSC","GBM","STAD","LIHC","LUAD","BLCA","OV","PAAD","PRAD","KIRC","SKCM")
length(intersect(cancer,PartMatch_cancer))
###绘制AUC曲线#############################################################################
library(survivalROC)
all_sur_score_LR<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/all_sur_score_LR.txt',stringsAsFactors = F,check.names = F)
all_sur_score_LR$years<-all_sur_score_LR$days_to_death/365
table(unlist(lapply(strsplit(all_sur_score_LR$sample,'-'),function(x) x[4])))
all_sur_score_LR$sample_type<-unlist(lapply(strsplit(all_sur_score_LR$sample,'-'),function(x) x[4]))
# all_sur_score_LR<-all_sur_score_LR[which(all_sur_score_LR$sample_type=='01A'),]
###得分按癌型归一化到01之间
all_sur_score_LR$scale_score<-lapply(unique(all_sur_score_LR$cancer),function(x){
  aaa<-all_sur_score_LR[all_sur_score_LR$cancer%in%x,'score']
  score<-aaa-min(aaa)
  score<-score/max(score)
  return(score)
}) %>% unlist()

library(timeROC)
library(survival)

all_sur_score_LR_match<-all_sur_score_LR[all_sur_score_LR$cancer%in%PartMatch_cancer,]
unique(all_sur_score_LR_match$cancer)
table(unlist(lapply(strsplit(all_sur_score_LR_match$sample,'-'),function(x) x[4])))

all_sur_score_LR_NOmatch<-all_sur_score_LR[!all_sur_score_LR$cancer%in%PartMatch_cancer,]
#运行timeROC函数，计算画图所需数据
p_data<-all_sur_score_LR
tROC <-timeROC(T=p_data$years,
               delta = p_data$status,
               marker = p_data$score,
               cause = 1,times = c(1,3,5),
               ROC=T)
#开始画图
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/TCGA_allCancer_ROC.pdf',width = 5.5,height = 5)
par(mar= c(5,5,1,1),cex.lab=1.5,cex.axis= 1.2) #设置图形边界
plot(tROC,time=1,col="#9EDAE5",title='TCGA',lwd=2) #1年ROC
plot(tROC,time=3,col="#1F77B4",add=T,title=F,lwd=2) #3年ROC
plot(tROC,time=5,col="#8C564B",add=T,title=F,lwd=2) #5年ROC
legend(0,1, # 这里设置legend的位置坐标：横坐标0,纵坐标1。也可以像前面一样设置为"bottomright"等
       c(paste0("AUC at 1 years  ", round(tROC$AUC[1], 2)),
         paste0("AUC at 3 years  ", round(tROC$AUC[2], 2)),
         paste0("AUC at 5 years  ", round(tROC$AUC[3], 2))),
       col=c("#9EDAE5","#1F77B4","#8C564B"),lwd=2,cex=1.2,bty="n")
dev.off()



##KM曲线
fit <- survfit(Surv(days_to_death,status) ~ group,data = all_sur_score_LR)
p1<-ggsurvplot(fit, data = all_sur_score_LR,pval = TRUE,legend.title = 'pan_cancer',
               palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
print(p1)

fit <- survfit(Surv(days_to_death,status) ~ group,data = all_sur_score_LR)
p1<-ggsurvplot(fit, data = all_sur_score_LR,pval = TRUE,legend.title = 'pan_cancer',xlab='Time(Day)',
               risk.table = TRUE,        # 增加risk table
               risk.table.col = "strata",# risk table根据分组使用不同颜色
               legend.labs = c("High", "Low"),    # 图例标签
               palette = c("#B42210", "#00448C"),linetype = "strata", surv.median.line = "hv",conf.int = F)
print(p1)
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allCancer_KM_LR_last.pdf',width = 5.5,height = 5.5)
print(p1)
dev.off()

KM_14<-all_sur_score_LR[all_sur_score_LR$cancer%in%PartMatch_cancer,]
fit <- survfit(Surv(days_to_death,status) ~ group,data = KM_14)
p1<-ggsurvplot(fit, data = KM_14,pval = TRUE,legend.title = 'pan_cancer',xlab='Time(Day)',
               risk.table = TRUE,        # 增加risk table
               risk.table.col = "strata",# risk table根据分组使用不同颜色
               legend.labs = c("High", "Low"),    # 图例标签
               palette = c("#B42210", "#00448C"),linetype = "strata", surv.median.line = "hv",conf.int = F)
print(p1)
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/TrainCancer_KM_LR_last.pdf',width = 5.5,height = 5.5)
print(p1)
dev.off()

KM_19<-all_sur_score_LR[!all_sur_score_LR$cancer%in%PartMatch_cancer,]
fit <- survfit(Surv(days_to_death,status) ~ group,data = KM_19)
p1<-ggsurvplot(fit, data = KM_19,pval = TRUE,legend.title = 'pan_cancer',xlab='Time(Day)',
               risk.table = TRUE,        # 增加risk table
               risk.table.col = "strata",# risk table根据分组使用不同颜色
               legend.labs = c("High", "Low"),    # 图例标签
               palette = c("#B42210", "#00448C"),linetype = "strata", surv.median.line = "hv",conf.int = F)
print(p1)
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/TextCancer_KM_LR_last.pdf',width = 5.5,height = 5.5)
print(p1)
dev.off()

sROC=survivalROC(Stime=all_sur_score_LR$years, # 生存时间
                 status=all_sur_score_LR$status, # 生存状态
                 marker = all_sur_score_LR$score, #选择gene87
                 predict.time =5, # 看5年的时间段
                 method="KM")
plot(sROC$FP, sROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col="red", 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.5, cex.axis=1.2, font=1.2)
abline(0,1)
aucText=paste0("5 years"," (AUC=",sprintf("%.3f",sROC$AUC),")")
legend("bottomright", aucText,
       lwd=2,bty="n",col=c("red","green","blue"),cex=1.2)





LR_sur_cox<-data.frame(cancer=c('pan_cancer',cancer),
                       logRankPvalue=p.val,
                       HR=HR,
                       HRup95=up95,
                       HRlow95=low95)
LR_sur_cox$sig<-'no_sig'
LR_sur_cox$sig[which(LR_sur_cox$logRankPvalue<0.05)]<-'sig'
write.table(LR_sur_cox,'E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/LR_sur_cox_allCancer_LR.txt',quote = F,sep = '\t',row.names = F)
LR_sur_cox<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/LR_sur_cox_allCancer_LR.txt',stringsAsFactors = F,check.names = F)
class(LR_sur_cox$logRankPvalue)

plot_LR_sur<-LR_sur_cox
plot_LR_sur$cancer<-factor(plot_LR_sur$cancer,levels = c(cancer,'pan_cancer'))
p_sur<-ggplot(plot_LR_sur, aes(x=cancer, y=HR, colour=sig)) + 
  geom_errorbar(aes(ymin=HRlow95, ymax=HRup95), width=0,linewidth=0.8)+
  geom_errorbar(aes(ymin=HR, ymax=HR), width=0.3,linewidth=0.8)+
  geom_hline(yintercept = 1,lty=2,col="grey30",lwd=0.6)+
  scale_color_manual(values = c("sig"="#c3262f","no_sig"="#5264b0"))+
  theme_bw() +#去掉背景灰色
  theme_classic()+
  labs(title = 'lasso_cox',y='HR (95% CI)',x='cancer')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))
print(p_sur)
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_HR.pdf',width = 8,height = 4)
print(p_sur)
dev.off()


pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_meta_All_cancer.pdf',width = 10,height = 9)
meta_p<-LR_sur_cox[-1,]
#meta_p<-meta_p[!meta_p$cancer%in%PartMatch_cancer,]
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(cancer))
forest(mg1, common = F,layout = "RevMan5")
grid.text(label = paste("P =", format.pval(meta_p$logRankPvalue, digits = 2)),
          x = 0.9, y = seq((1-0.15), 0.15, length.out = nrow(meta_p)),
          gp = gpar(fontsize = 12))
dev.off()

pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_meta_other_cancer.pdf',width = 10,height = 6)
meta_p<-LR_sur_cox[-1,]
meta_p<-meta_p[!meta_p$cancer%in%PartMatch_cancer,]
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(cancer))
forest(mg1, common = F,layout = "RevMan5")
grid.text(label = paste("P =", format.pval(meta_p$logRankPvalue, digits = 2)),
          x = 0.9, y = seq(0.8, 0.2, length.out = nrow(meta_p)),
          gp = gpar(fontsize = 12))
dev.off()


pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_meta_match_cancer.pdf',width = 10,height = 6)
meta_p<-LR_sur_cox[LR_sur_cox$cancer%in%PartMatch_cancer,]
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(cancer))
forest(mg1, common = F,layout = "RevMan5")
grid.text(label = paste("P =", format.pval(meta_p$logRankPvalue, digits = 2)),
          x = 0.9, y = seq(0.72, 0.28, length.out = nrow(meta_p)),
          gp = gpar(fontsize = 12))
dev.off()



cancer<-LR_sur_cox$cancer
PartMatch_cancer<-c("BRCA","CESC","COAD","HNSC","GBM","STAD","LIHC","LUAD","BLCA","OV","PAAD","PRAD","KIRC","SKCM")
match_cancer<-c('brca','cesc','gbm','lihc','luad','ov','prad','skcm') %>% toupper()###名称完全对应的癌症


plot_LR_sur<-LR_sur_cox[!LR_sur_cox$cancer%in%match_cancer,]
plot_LR_sur$cancer<-factor(plot_LR_sur$cancer,levels = match_cancer)
plot_LR_sur$cancer<-factor(plot_LR_sur$cancer,levels = c(plot_LR_sur$cancer[2:nrow(plot_LR_sur)],'pan_cancer'))
p_sur<-ggplot(plot_LR_sur, aes(x=cancer, y=HR, colour=sig)) + 
  geom_errorbar(aes(ymin=HRlow95, ymax=HRup95), width=0,linewidth=0.8)+
  geom_errorbar(aes(ymin=HR, ymax=HR), width=0.3,linewidth=0.8)+
  geom_hline(yintercept = 1,lty=2,col="grey30",lwd=0.6)+
  scale_color_manual(values = c("sig"="#c3262f","no_sig"="#5264b0"))+
  theme_bw() +#去掉背景灰色
  theme_classic()+
  labs(title = 'lasso_cox',y='HR (95% CI)',x='cancer')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))
print(p_sur)
pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_HR_NOmatch_cancer.pdf',width = 10,height = 3)
print(p_sur)
dev.off()


meta_p<-LR_sur_cox[LR_sur_cox$cancer%in%PartMatch_cancer,]
meta_p<-LR_sur_cox[-1,]
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(cancer))

pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/allCancer/allcancer_meta_match_cancer.pdf',width = 10,height = 6)
forest(mg1, common = F,layout = "RevMan5")
library(grid)
grid.text(label = paste("P =", format.pval(meta_p$logRankPvalue, digits = 2)),
          x = 0.9, y = seq(0.7, 0.3, length.out = nrow(meta_p)),
          gp = gpar(fontsize = 12))
dev.off()






