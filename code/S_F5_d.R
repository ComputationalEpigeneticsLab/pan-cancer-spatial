#####使用MP在core区的均值对样本进行pca降维 tsne聚类
library(Rtsne)

cluster_re<-read.csv('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/pam_pearson/pam_pearson.k=5.consensusClass.csv',
                     header = F,stringsAsFactors = F,check.names = F)
cluster_re$V1<-gsub('hn.as','hn-as',cluster_re$V1)
cluster_re$V1<-gsub('\\.','/',cluster_re$V1)
cluster_re$V2<-paste0("Cluster",cluster_re$V2)
table(cluster_re$V2)

mean_MPscore<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/ConsensusCluster/merge_by_slice_sum&mean_core.txt',
                         stringsAsFactors = F,check.names = F)
mean_MPscore<-mean_MPscore[,cluster_re$V1]



aa<-as.matrix(t(mean_MPscore[,cluster_re$V1]))
###pca
pca_out=prcomp(as.matrix(aa))
#pca_out$x[1:3,1:3]
pc_s_n=ceiling(dim(pca_out$x)[2]/3)
# 使用Rtsne函数进行tSNE降维分析
set.seed(1181)
tsne_out <- Rtsne(pca_out$x[,1:pc_s_n],
                  pca=FALSE,
                  #dims=3,
                  check_duplicates = F,
                  perplexity=20,### perplexity<(nrow(tsne_exp)-1)/3
                  theta=0.0)
jiangwei <- as.data.frame(tsne_out$Y)
jiangwei$cluster<-cluster_re$V2

p_E<-ggplot(jiangwei, aes(x=V1, y=V2,color=cluster)) +
  geom_hline(aes(yintercept=0), colour="grey", linetype="dashed") +
  geom_vline(aes(xintercept=0), colour="grey", linetype="dashed")+
  geom_point(aes(colour=cluster),size=1.5) +
  scale_color_manual(values = c('Cluster1'="#469AD4",'Cluster2'="#8DD1C6",'Cluster3'="#F9E30E",
                                'Cluster4'='#968BD7','Cluster5'='#E2806F'))+ #设置填充颜色
  labs(title = "Cluster_MP",
       x = "tSNE1",
       y = "tSNE2") +
  theme_classic()+
  guides(colour = guide_legend(override.aes = list(size=3)))
print(p_E)
pdf('E:/Mirror/ST_analysis/pic/NMF_module/tSNE_cluster_MP.pdf',width = 5,height = 4)
print(p_E)
dev.off()