#######生存模型验证
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
library(Hmisc)


#coef_re<-read.delim('/data/zhouweiwei/ST_analysis/data/cellchatTCGA/lasso_cox_coef.txt',stringsAsFactors = F,check.names = F)
# coef_re<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/lasso_cox_coef.txt',stringsAsFactors = F,check.names = F)
# coef_re<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/lasso_cox_coef_allCancer.txt',stringsAsFactors = F)
coef_re<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/lasso_cox_coef_LR.txt',stringsAsFactors = F,check.names = F)
#coef_re<-coef_re[which(coef_re$coef>0),]
coef_re<-coef_re[which(coef_re$`Pr(>|z|)`<0.01),]
coef_func<-function(x){
  score<-sum(x*coef_re$coef)
  return(score)
}

dir_GEO<-'E:/Mirror/ST_analysis/data/sur_GEO/last/'
file_exp<-list.files(pattern = 'exp.txt',path = dir_GEO,recursive = T)
file_pheno<-list.files(pattern = 'pheno.txt',path = dir_GEO,recursive = T)
gseName<-unlist(lapply(strsplit(file_pheno,'/'),function(x)x[1]))

GEO_order<-read.csv(paste0(dir_GEO,'GEO_last.csv'),stringsAsFactors = F,check.names = F)

dir_pic<-'E:/Mirror/ST_analysis/pic/cellchatTCGA/GEO_sur/allCancer/last/'
p.val<-c()
HR<-c()
up95<-c()
low95<-c()
all_sur_data<-c()
pdf(paste0(dir_pic,'all_KM_LR.pdf'),width = 4,height = 3)
for(i in 1:length(file_pheno)){
  #i=19
  gse_exp<-read.delim(paste0(dir_GEO,file_exp[i]),stringsAsFactors = F,check.names = F)
  gse_pheno<-read.delim(paste0(dir_GEO,file_pheno[i]),stringsAsFactors = F,check.names = F)
  gse_pheno<-gse_pheno[which(gse_pheno$OS.time>=30),]
  inter_sample<-intersect(rownames(gse_pheno),colnames(gse_exp))
  gse_exp<-gse_exp[,inter_sample]
  gse_pheno<-gse_pheno[inter_sample,]
  
  gse_exp[,1]
  
  inter_gene<-intersect(rownames(coef_re),rownames(gse_exp))
  coef_func<-function(x){
    score<-sum(x*coef_re[inter_gene,]$coef,na.rm = T)
    return(score)
  }
  
  sur_data<-data.frame(sample=rownames(gse_pheno),row.names = rownames(gse_pheno),
                       days_to_death=gse_pheno$OS.time,
                       status=gse_pheno$OS,
                       score=apply(gse_exp[inter_gene,],2,coef_func),
                       group='Not')
  #sur_data$group[which(sur_data$score<median(sur_data$score))]<-'Low'
  optimalCutoff=surv_cutpoint(sur_data, time = "days_to_death", event = "status",
                              variables = c("score"))
  res.cut<- surv_categorize(optimalCutoff)
  sur_data$group<-str_to_title(as.character(res.cut$score))
  #sur_data$group<-str_to_title(sur_data$group)
  table(sur_data$group)
  #?capitalize
  sur_data$GSE<-gseName[i]
  print(nrow(sur_data))
  all_sur_data<-rbind(all_sur_data,sur_data)
  
  fit <- survfit(Surv(days_to_death,status) ~ group,data = sur_data)
  p1<-ggsurvplot(fit, data = sur_data,pval = TRUE,legend.title = gseName[i],
                 palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
  print(p1)
  
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

write.table(all_sur_data,'E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/GEO_allSample_sur_data_LR_last.txt',quote = F,sep = '\t')
all_sur_data<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/GEO_allSample_sur_data_LR_last.txt',
                         stringsAsFactors = F,check.names = F)
# table(all_sur_data$GSE)
# 
# all_sur_data<-all_sur_data[all_sur_data$GSE%in%c('GSE143985','GSE161158','GSE17536','GSE29621',
#                                                  'GSE31210','GSE39084','GSE78220'),]#,'GSE96058'
fit <- survfit(Surv(days_to_death,status) ~ group,data = all_sur_data)
p1<-ggsurvplot(fit, data = all_sur_data,pval = TRUE,legend.title = 'all_sample',
               palette = c("#E7B800", "#2E9FDF"),linetype = "strata", surv.median.line = "hv",conf.int = TRUE)
print(p1)
# pdf('E:/Mirror/ST_analysis/pic/cellchatTCGA/GEO_sur/allSample_KM.pdf',width = 8,height = 7)
# print(p1)
# dev.off()

if(length(which(all_sur_data$group=='High'))>10&length(which(all_sur_data$group=='Low'))>10){
  sdiff <- survdiff(Surv(days_to_death,status) ~ group,data = all_sur_data)
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

dev.off()

num_sample<-as.data.frame(table(all_sur_data$GSE))
num_sample<-num_sample[match(GEO_order$GEO,num_sample$Var1),]
GEO_order$sample_num<-num_sample$Freq
write.csv(GEO_order,paste0(dir_GEO,'GEO_last_sampleNUM.csv'),quote = F,row.names = F)

LR_sur_cox_GEO<-data.frame(data=c(gseName,'pan_cancer'),
                           logRankPvalue=p.val,
                           HR=HR,
                           HRup95=up95,
                           HRlow95=low95)
LR_sur_cox_GEO$sig<-'no_sig'
LR_sur_cox_GEO$sig[which(LR_sur_cox_GEO$logRankPvalue<0.05)]<-'sig'
write.table(LR_sur_cox_GEO,'E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/LR_sur_cox_GEO_LR_last.txt',quote = F,sep = '\t',row.names = F)
LR_sur_cox_GEO<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchatTCGA/LR_sur_cox_GEO_LR_last.txt',stringsAsFactors = F)



library(meta) #加载包
meta_p<-LR_sur_cox_GEO[-nrow(LR_sur_cox_GEO),]
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(data))
#forest(mg1)

forest(mg1, common = F,layout = "RevMan5") # layout = "meta", "RevMan5", "JAMA", or "subgroup"

pdf(paste0(dir_pic,'meta_HR_LR_last.pdf'),width = 10,height = 7.2)
meta_p<-LR_sur_cox_GEO[-nrow(LR_sur_cox_GEO),]
meta_p<-meta_p[match(GEO_order$GEO,meta_p$data),]
meta_p$data<-paste0(1:nrow(meta_p),'_',meta_p$data)
mg1 <- metagen(log(HR),lower = log(HRlow95), upper = log(HRup95),
               data = meta_p, sm="HR", studlab = paste(data))
forest(mg1, common = F,layout = "RevMan5",
       text.pvalue = format.pval(meta_p$logRankPvalue, digits = 3, eps = 0.001))
grid.text(label = paste("P =", format.pval(meta_p$logRankPvalue, digits = 2)),
          x = 0.85, y = seq((1-0.29), 0.29, length.out = nrow(meta_p)),
          gp = gpar(fontsize = 12))
dev.off()


