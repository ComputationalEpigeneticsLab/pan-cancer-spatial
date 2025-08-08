#####对MP score一致性聚类结果
###找每个类中主导的MP
###类似于DEG 计算每个类的MP与其他类的FC和p值(秩和检验或T检验)







cluster_re<-read.csv('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/pam_pearson/pam_pearson.k=5.consensusClass.csv',
                     header = F,stringsAsFactors = F,check.names = F)
cluster_re$V1<-gsub('hn.as','hn-as',cluster_re$V1)
cluster_re$V1<-gsub('\\.','/',cluster_re$V1)
table(cluster_re$V2)

MP_score<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/merge_by_slice_sum&mean_core.txt',
                     stringsAsFactors = F,check.names = F)
length(intersect(cluster_re$V1,colnames(MP_score)))
MP_score<-MP_score[,cluster_re$V1]
MP_score<-as.matrix(MP_score)
cluster_num<-unique(cluster_re$V2)

all_C_DE_MP<-data.frame(MP=rownames(MP_score))
for(i in 1:length(cluster_num)){
  #i=1
  is_site<-which(cluster_re$V2==cluster_num[i])
  no_site<-which(cluster_re$V2!=cluster_num[i])
  
  DE_MP<-apply(MP_score,1,function(x){##x=MP_score[1,]
    #x<-as.vector(as.matrix(x))
    FC<-mean(x[is_site])/mean(x[no_site])
    pValue_W<-wilcox.test(x[is_site],x[no_site])[["p.value"]]
    pValue_T<-t.test(x[is_site],x[no_site])[["p.value"]]
    return(c(FC,pValue_W,pValue_T))
  }) %>% t()
  colnames(DE_MP)<-paste0('C',cluster_num[i],c('_FC','_W_P','_T_P'))
  all_C_DE_MP<-cbind(all_C_DE_MP,DE_MP)
}

write.table(all_C_DE_MP,'E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/all_DE_MP.txt',quote = F,sep = '\t',
            row.names = F)


#####MP聚类的5簇的差异MP
library(Seurat)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(ggplot2)
library(ggrepel)


###文件来自181_MP_cluster主导.R
all_C_DE_MP<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/all_DE_MP.txt',
                        stringsAsFactors = F,check.names = F)

cluster<-paste0('C',1:5)


####火山图
###bdy


for(i in 1:length(cluster)){#i=1
  v_data<-all_C_DE_MP[,c(1,grep(cluster[i],colnames(all_C_DE_MP)))]
  colnames(v_data)<-c('MP','FC','W_P','T_P')
  v_data$LogFC<-as.numeric(log2(v_data$FC))
  v_data$FDR<-as.numeric(p.adjust(v_data$W_P))
  v_data$type = as.factor(ifelse(v_data$FDR< 0.05 & abs(v_data$LogFC) >= log2(1.5), 
                                 ifelse(v_data$LogFC> log2(1.5) ,'Up','Down'),'Not'))
  table(v_data$type)
  v_data$label <- v_data$MP
  v_data$label[!v_data$type%in%c('Up','Down')]<-""
  length(which(v_data$label!=""))
  
  a<-ggplot(v_data, aes(x = LogFC, y = -log10(FDR), colour=type)) +
    geom_point(alpha=0.8, size=2) +
    scale_color_manual(values=c('Down'="#4393C3",'Not'="#808080",'Up'="#D6604D"))+
    geom_text_repel(aes(x = LogFC,y = -log10(FDR),label=label))+
    geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=2,col="black",lwd=0.6) +
    geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)+
    xlab('Log2FC')+
    ylab('-Log10(P-Value)')+
    ggtitle(cluster[i])+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position="right",
          legend.title = element_blank())
  #print(a)
  pdf(paste0('E:/Mirror/ST_analysis/pic/NMF_module/DE_MP_volcano/',cluster[i],'_DE_MP_volcano.pdf'),width = 5.5,height = 5)
  print(a)
  dev.off()
}



#?geom_label_repel


###火山图整成一个
rownames(all_C_DE_MP)<-all_C_DE_MP$MP
all_C_DE_MP<-all_C_DE_MP[,-1]

FC_data<-all_C_DE_MP[,grep('FC',colnames(all_C_DE_MP))]
colnames(FC_data)<-unlist(lapply(strsplit(colnames(FC_data),'_'),function(x)x[1]))
P_data<-all_C_DE_MP[,grep('W_P',colnames(all_C_DE_MP))]
P_data<-apply(P_data,2,p.adjust)
FC_data<-reshape2::melt(as.matrix(FC_data))
P_data<-reshape2::melt(as.matrix(P_data))
colnames(FC_data)<-c('MP','cluster','FC')
FC_data$p_value<-P_data$value
FC_data$logFC<-log2(FC_data$FC)
FC_data$logP<-(-log10(FC_data$p_value))
FC_data$lab<-paste0(FC_data$cluster,'_',FC_data$MP)
FC_data$lab[which(FC_data$logFC<log2(1.5)|FC_data$p_value>=0.05)]<-''
FC_data$cluster<-as.vector(FC_data$cluster)
FC_data$cluster[which(abs(FC_data$logFC)<log2(1.5)|FC_data$p_value>=0.05)]<-'no_sig'


a<-ggplot(FC_data, aes(x = logFC, y = logP, colour=cluster)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values=c('C1'="#469AD4",'C2'="#8DD1C6",'C3'="#F9E30E",'C4'='#968BD7','C5'='#E2806F','no_sig'='grey80'))+
  geom_text_repel(aes(x = logFC,y = logP,label=lab),size=4,force=T)+
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=2,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)+
  xlab('Log2FC')+
  ylab('-Log10(P-Value)')+
  ggtitle('all_MP_cluster')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",
        legend.title = element_blank())
print(a)
pdf(paste0('E:/Mirror/ST_analysis/pic/NMF_module/DE_MP_volcano/','all_MP','_DE_MP_volcano.pdf'),width = 5.5,height = 5)
print(a)
dev.off()






