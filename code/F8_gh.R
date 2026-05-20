
##ROC一起可视化################
color=read.delim('E:/panCancerST/35cancer_color.txt',check.names = F)
GEO_color=color$color
names(GEO_color)=color$cancer

TCGA_color=c("BLCA"='#D2BE18',"ACC"='#AA9C1F',"BRCA"='#b7996d', "CESC"='#B51672', "CHOL"='#DA9B8B',
       "COAD"='#8bc96d', "DLBC"='#3B62A3', "ESCA"='#222F63', "GBM"='#a4cde1',
       "HNSC"='#ECD0DC', "KICH"='#AF4F91', "KIRC"='#af93c4', "KIRP"='#9C80C7', "LAML"='#DCB059',
       "LGG"='#B26A1E',  "LIHC"='#c6b598', "LUAD"='#7A9AC5', "LUSC"='#85B3AB', "MESO"='#A4C5C1',
       "OV"='#BDA684',   "PAAD"='#8B3580', "PCPG"='#C194BF', "PRAD"='#f9b769', "READ"='#649657',
       "SARC"='#3D7B2E', "SKCM"='#31307F', "STAD"='#FFB978', "TGCT"='#96DD88', "THCA"='#325FA2',
       "THYM"='#FF7F0E', "UCEC"='#FF9694', "UCS" ='#D42728', "UVM"='#A840FA')

###
result$type = unlist(lapply(result$sample,function(x){strsplit(x,'-')[[1]][4]}))
result = result[result$type%in%c('01A','11A','11B','11C'),]
allCancer_data=result
###
cancerlist=GEO$cancer
allCancer_data=all_cancer
aLL_AUC_score<-c()
pdf("E:/panCancerST/bud_model/GEO_classifysample_ROC1.pdf",width = 6,height = 6)
i=1
all_score_cancer<-allCancer_data[allCancer_data$cancer%in%cancerlist[i],]
AUC_data<-data.frame(score=all_score_cancer$score,
                     lab=all_score_cancer$sampletype)
if(length(which(AUC_data$lab=='normal'))>2){
  AUC_data$lab<-factor(AUC_data$lab,levels = c('normal','tumor'))
  AUC_score<-as.numeric(roc(AUC_data$lab,AUC_data$score)$auc)
  names(AUC_score)<-cancerlist[i]
  aLL_AUC_score<-c(aLL_AUC_score,AUC_score)
  plot(roc(AUC_data$lab,AUC_data$score),
       col=GEO_color[cancerlist[i]],
       legacy.axes = TRUE,
       xlim=c(1,0),
       xlab = "1-Specificity", ylab = "Sensitivity")
}  

for(i in 2:length(cancerlist)){
  all_score_cancer<-allCancer_data[allCancer_data$cancer%in%cancerlist[i],]
  AUC_data<-data.frame(score=all_score_cancer$score,
                       lab=all_score_cancer$sampletype)
  if(length(which(AUC_data$lab=='normal'))>2){
    AUC_data$lab<-factor(AUC_data$lab,levels = c('normal','tumor'))
    AUC_score<-as.numeric(roc(AUC_data$lab,AUC_data$score)$auc)
    names(AUC_score)<-cancerlist[i]
    aLL_AUC_score<-c(aLL_AUC_score,AUC_score)
    plot(roc(AUC_data$lab,AUC_data$score),
         col=GEO_color[cancerlist[i]],
         legacy.axes = TRUE,
         xlim=c(1,0),
         xlab = "1-Specificity", ylab = "Sensitivity",add=T)
  }
}

legend_lab<-paste0(names(aLL_AUC_score),' (',round(aLL_AUC_score,3),')')
legend("bottomright",
       legend=legend_lab,
       col=GEO_color[names(aLL_AUC_score)],
       lwd=1,bty='n',inset=c(0.28,0))
dev.off()

write.csv(as.data.frame(aLL_AUC_score),'E:/panCancerST/bud_model/TCGA_ROC.csv')

##ROC值箱式图####
plotdata=GEO
plotdata=read.csv('E:/panCancerST/bud_model/TCGA_ROC.csv')
plotdata$database="GEO"

# 绘制图形（保留统计线，移除散点中线）
e=ggplot(plotdata, aes(x = database, y = ROC)) +
  stat_boxplot(geom = "errorbar",width = 0.3,position = position_dodge(0.9)) +#绘制统计参考线（误差棒+分位数线）
  geom_jitter(aes(color = cancer),size = 1,alpha = 1,show.legend = TRUE,width = 0.15,height = 0) +
  stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8),show.legend = F) + # 误差棒，中位数，25%和75%分位数
  stat_summary(fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8),show.legend = F) + # 中位数水平线
  scale_color_manual(values = GEO_color) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x="" ,y="ROC Value") +
  theme_classic(base_size = 20)+
  theme(axis.title.x = element_text(size=12,colour = "black"),
        axis.text.x = element_text(size=10,colour = "black"),
        axis.title.y =element_text(size=12,colour = "black"),
        axis.text.y = element_text(size = 10,colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.3,"cm"))
e
ggsave( "E:/panCancerST/bud_model/GEO_ROC箱式图.pdf",e,width = 4.8,height = 4)



##整合所有的表达值###############################################################
GEO = read.xlsx('E:/panCancerST/bud_model/GEOBulk_ROC.xlsx',sheet = 3)
gene=read.delim("E:/panCancerST/bud_model/50_lasso.txt",check.names = FALSE)

all_cancer = data.frame()
for (i in 1:nrow(GEO)) {
  meta = read.csv(paste0('E:/panCancerST/GEO/',GEO$cancer[i],'/',GEO$dataset[i],'/',GEO$dataset[i],'_cli.csv'))
  ExpMatrix = read.delim(paste0('E:/panCancerST/GEO/',GEO$cancer[i],'/',GEO$dataset[i],'/',GEO$dataset[i],'_symbol.txt'),check.names = FALSE)
  ExpMatrix = ExpMatrix[,meta$sampleID]
  
  matchgene=intersect(rownames(ExpMatrix),rownames(gene))
  exp_filter=ExpMatrix[matchgene,] %>% t() %>% as.data.frame()
  gene_filter=gene[matchgene,]
  
  coef_vector=gene_filter$s0
  result <- data.frame(
    sample = rownames(exp_filter),
    score = as.matrix(exp_filter) %*% coef_vector
  )
  result=result[match(meta$sampleID,rownames(result)),]
  result$sampletype=meta$sampletype
  result$cancer=GEO$cancer[i]
  all_cancer=rbind(all_cancer,result)
}
write.table(all_cancer,file = 'E:/panCancerST/bud_model/GEO_score.txt',quote = F,sep = '\t',row.names = F)









