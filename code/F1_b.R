####地球图
library(plot1cell)
library(Seurat)
library(tidyverse)
library(stringr)
library(RColorBrewer)

all_umap<-read.delim('/pan_cancer/ST_merge_res0.8_metadatainfor_BdyCoreBud.txt',
                     stringsAsFactors = F,check.names = F)
integration_umap<-all_umap[,c(5,6,4,8,10)]


library(entropy)
###熵值计算
re_entropy<-lapply(unique(all_umap$leiden), function(x){###x=0
  aa<-all_umap[all_umap$leiden%in%x,]
  bb<-as.data.frame(table(aa$cancer))
  return(entropy(bb$Freq))
}) %>% unlist()


re_entropy<-data.frame(leiden=unique(all_umap$leiden),
                       entropy_score=re_entropy)

all_umap_a<-merge(all_umap,re_entropy,by='leiden')
rownames(all_umap_a)<-paste0(all_umap_a$batch,'_',all_umap_a$cell_name)
all_umap_a<-all_umap_a[rownames(all_umap),]

cancer_col<-read.delim('/pan_cancer/35cancer_color.txt',stringsAsFactors = F,check.names = F)
cancer_col2<-cancer_col$color
names(cancer_col2)<-cancer_col$cancer

####区域
p_spot<-ggplot(data = all_umap_a,aes(x = UMAP_1, y = UMAP_2,color=FinalLocalType)) + 
  geom_point(size=1)+
  scale_color_manual(values =  c("Boundary"="#F6B86D","Core"="#D62D28","Budding"="#EE762D",'Immune'='#a9d38a','Normal'='#2375ae'))+
  theme_classic()+
  ggtitle('cancer')+
  guides(colour = guide_legend(override.aes = list(size=2)))
print(p_spot)
png('/pic/FinalLocalType.png',width = 700,height = 600)
print(p_spot)
dev.off()

###癌症类型
p_spot2<-ggplot(data = all_umap_a,aes(x = UMAP_1, y = UMAP_2,color=cancer)) + 
  geom_point(size=1)+
  scale_color_manual(values =  cancer_col2)+
  theme_classic()+
  ggtitle('cancer')+
  guides(colour = guide_legend(override.aes = list(size=2)))
print(p_spot2)
png('/pic/cancer_type.png',width = 700,height = 600)
print(p_spot2)
dev.off()

###恶性与正常
table(all_umap_a$FinalLocalType)
all_umap_a$mal<-'Malignant'
all_umap_a$mal[all_umap_a$FinalLocalType%in%c('Immune',   'Normal')]<-'Non-Malignant'
table(all_umap_a$mal)
p_spot2<-ggplot(data = all_umap_a,aes(x = UMAP_1, y = UMAP_2,color=mal)) + 
  geom_point(size=1)+
  scale_color_manual(values =  c('Malignant'='#C92C2A','Non-Malignant'='#559499'))+
  theme_classic()+
  ggtitle('mal')+
  guides(colour = guide_legend(override.aes = list(size=2)))
print(p_spot2)
png('/pic/Mal_type.png',width = 680,height = 600)
print(p_spot2)
dev.off()


###熵值
p_spot2<-ggplot(data = all_umap_a,aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(colour=entropy_score),size=1) +
  scale_color_gradientn(colours = c(colorRampPalette(c("#E1CFC4","#F39C67"))(50),###数字可以改
                                    colorRampPalette(c("#F39C67","#B20A1C"))(50))
  )+
  theme_classic()+
  labs(title = 'entropy',x = "",y = "",col='entropy_score')
#print(p_spot2)
png('/pic/entropy.png',width = 680,height = 600)
print(p_spot2)
dev.off()




