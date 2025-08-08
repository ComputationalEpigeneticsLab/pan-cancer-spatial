#####不同ecosysytem的spot间的免疫通路秩和检验
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

integration_umap_MP<-read.delim('E:/Mirror/ST_analysis/data/pic_data/integration_umap_MP.txt',stringsAsFactors = F,check.names = F)

dir_EMT<-'E:/Mirror/ST_analysis/data/10X Visium/EMT/'
file_Imm<-list.files(pattern = 'MP5_score.txt',path = dir_EMT,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_Imm,'_MP5'),function(x)x[1]))
###数据来源191_MP聚类基因基因得分.R


all_FC_re<-c()
all_P_re<-c()
for(i in 1:length(file_Imm)){
  #i=1
  Imm_data<-read.delim(paste0(dir_EMT,file_Imm[i]),stringsAsFactors = F,check.names = F)
  slice_MP<-integration_umap_MP[integration_umap_MP$slice%in%dataSlice[i],]
  rownames(slice_MP)<-slice_MP$cell_name
  slice_MP<-slice_MP[rownames(Imm_data),]
  Imm_data$top_MP<-slice_MP$MP_top1
  Imm_data<-Imm_data[which(Imm_data$LocalType!='not.defined'),]
  Imm_data$ecosystem<-''
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_8','MP_10','MP_11')]<-'ecosystem1'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_2','MP_14')]<-'ecosystem2'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_3','MP_6','MP_7','MP_9','MP_12')]<-'ecosystem3'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_4')]<-'ecosystem4'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_5','MP_13')]<-'ecosystem5'
  #table(Imm_data$top_MP,Imm_data$ecosystem)
  
  all_DE_Imm<-data.frame(ImmPath=colnames(Imm_data)[1:18])
  for(j in paste0('ecosystem',1:5)){
    #j='ecosystem1'
    is_site<-which(Imm_data$ecosystem==j)
    no_site<-which(Imm_data$ecosystem!=j)
    
    DE_Imm<-apply(Imm_data[,1:18],2,function(x){##x=Imm_data[,1]
      #x<-as.vector(as.matrix(x))
      if(length(which(x%in%NA==T))>0) x[x%in%NA]<-0
      FC<-2^(mean(x[is_site])-mean(x[no_site]))
      pValue_W<-wilcox.test(x[is_site],x[no_site])[["p.value"]]
      if(pValue_W%in%NaN) pValue_W<-1
      #pValue_T<-t.test(x[is_site],x[no_site])[["p.value"]]
      return(c(FC,pValue_W))
    }) %>% t()
    colnames(DE_Imm)<-paste0(j,c('_FC','_W_P'))
    all_DE_Imm<-cbind(all_DE_Imm,DE_Imm)
  }
  all_DE_Imm<-all_DE_Imm[,-1]
  
  FC_re<-all_DE_Imm[,grep('FC',colnames(all_DE_Imm))]
  P_re<-all_DE_Imm[,grep('W_P',colnames(all_DE_Imm))]
  FC_re<-melt(as.matrix(FC_re))
  P_re<-melt(as.matrix(P_re))
  FC_re$slice<-dataSlice[i]
  P_re$slice<-dataSlice[i]
  
  all_FC_re<-rbind(all_FC_re,FC_re)
  all_P_re<-rbind(all_P_re,P_re)
  
  print(dataSlice[i])
}

table(paste0(all_FC_re$Var1,all_FC_re$slice)==paste0(all_P_re$Var1,all_P_re$slice))
all_P_re$FC<-all_FC_re$value
colnames(all_P_re)[1:3]<-c('ImmPath','ecosystem','P_value')
all_P_re$ecosystem<-unlist(lapply(strsplit(as.vector(all_P_re$ecosystem),'_W'),function(x)x[1]))
write.table(all_P_re,'E:/Mirror/ST_analysis/data/pic_data/ecosystem_DE_ImmPath.txt',quote = F,sep = '\t',row.names = F)

all_P_re<-read.delim('E:/Mirror/ST_analysis/data/pic_data/ecosystem_DE_ImmPath.txt',stringsAsFactors = F,check.names = F)


dfcol<-data.frame(ImmPath=unique(all_P_re$ecosystem),
                  logFC=0,
                  label=rep('',5))
mycol <- c('#4b2d80','#136599','#1c959d','#8e2823','#d1e181')

all_P_re$col<-as.vector(all_P_re$ImmPath)
all_P_re$logFC<-log2(all_P_re$FC)
all_P_re$col[which(all_P_re$P_value>=0.05|abs(all_P_re$logFC)<log2(1.2))]<-'no_sig'


all_P_re<-all_P_re[order(all_P_re$ecosystem),]
all_P_re<-lapply(unique(all_P_re$ecosystem),function(x){#x=unique(all_P_re$ecosystem)[5]
  aa<-all_P_re[all_P_re$ecosystem%in%x,]
  aa$lab_cand<-paste0(aa$slice,'of',aa$ImmPath)
  # aa_down<-aa[aa$col%in%'down',]
  # aa_down<-aa_down[order(aa_down$logFC),]
  # down_slice<-aa_down$slice[1:ifelse(nrow(aa_down)>2,2,nrow(aa_down))]
  # aa_up<-aa[aa$col%in%'up',]
  # aa_up<-aa_up[order(aa_up$logFC,decreasing = T),]
  # up_slice<-aa_up$slice[1:ifelse(nrow(aa_up)>2,2,nrow(aa_up))]
  aa<-aa[order(aa$FC,decreasing = T),]
  lab_slice<-aa$lab_cand[which(aa$P_value<0.05)]
  lab_slice<-lab_slice[1:ifelse(length(lab_slice)>2,2,length(lab_slice))]
  
  aa$lab<-''
  aa$lab[match(lab_slice,aa$lab_cand)]<-lab_slice
  return(aa)
})
all_P_re<-do.call(rbind,all_P_re)
all_P_re$lab<-gsub('slice','',all_P_re$lab)

#根据图p中log2FC区间确定背景柱长度：
max_bar<-lapply(unique(all_P_re$ecosystem),function(x){
  aa<-all_P_re$logFC[all_P_re$ecosystem%in%x]
  return(max(aa))
}) %>% unlist()
min_bar<-lapply(unique(all_P_re$ecosystem),function(x){
  aa<-all_P_re$logFC[all_P_re$ecosystem%in%x]
  return(min(aa))
}) %>% unlist()
dfbar<-data.frame(ecosystem=unique(all_P_re$ecosystem),
                  logFC=max_bar)
dfbar1<-data.frame(ecosystem=unique(all_P_re$ecosystem),
                   logFC=min_bar)
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = ecosystem,y = logFC),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = ecosystem,y = logFC),
           fill = "#dcdcdc",alpha = 0.6)

path_col<-c('#a4cde1','#277fb8','#96cb8f','#4dae47','#b79973','#f38989','#e32427','#003366','#f9b769',
            '#d4a6a8','#815e99','#b05a28','#663399','#FFCC00','#EEF0A5','#336633','#FC733D','#CCCC99')
names(path_col)<-unique(as.vector(all_P_re$ImmPath))
###加散点图和标签
p2<-p1+geom_jitter(data = all_P_re,aes(x = ecosystem,y = logFC,color = col),
                   size=0.8,position = position_jitter(seed = 1))+
  # scale_fill_manual(values = c("#f79f1f","#a3cb38","#1289a7"))+
  geom_text_repel(data = all_P_re,aes(x = ecosystem,y = logFC,label = ifelse(lab!='', lab, "")), position = position_jitter(seed = 1),
                  size=2.5,point.padding = 0,direction ='both')+
  scale_color_manual(values = path_col)+
  theme_bw()+
  labs(title = 'ecosystem_ImmPath',x = "",y = "logFC",col='ImmPath')+
  theme(#axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
    panel.background = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 10,colour = "black"),
    axis.title = element_text(size = 15))+
  theme(panel.grid = element_blank(),legend.position = 'right')+
  guides(colour = guide_legend(override.aes = list(size=3.5)))

###加x轴色块标签
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=ImmPath,y=logFC),
                     height=0.1,
                     color = "black",
                     fill = mycol,
                     alpha = 0.95,
                     show.legend = F)
print(p3)
pdf('E:/Mirror/ST_analysis/pic/NMF_module/ecosystem_DE_Imm/ecosystem_path_point_bar.pdf',width = 12,height = 6)
print(p3)
dev.off()

