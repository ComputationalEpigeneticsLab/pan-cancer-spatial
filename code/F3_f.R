###按MP为spot分ecosystem 并计算不同eco的差异免疫通路
library(dplyr) #进行数据转换
library(ggplot2)#绘图
library(enrichplot)
library(fgsea)
library(stats)
library(Seurat)
library(stringr)
library(plyr)
library(ggrepel)

integration_umap_MP<-read.delim('/10X_Visium/new_st_NMF/MPscore_new/all_MP.txt',
                                stringsAsFactors = F,check.names = F)
integration_umap_MP[1:3,]
eco_MP<-integration_umap_MP[,c(1:17)]
eco_MP[1:3,]
eco_MP$top_MP<-apply(eco_MP,1,function(x){
  names(x)[which.max(x)]
})
table(eco_MP$top_MP)
eco_MP[1:3,]
dim(eco_MP)
eco_MP$ecosystem<-''
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_8')]<-'ecosystem1'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_6')]<-'ecosystem2'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem3'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_7','MP_17')]<-'ecosystem4'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_1','MP_4')]<-'ecosystem5'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem6'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_15','MP_11')]<-'ecosystem7'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_8')]<-'ecosystem1'
# eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_6')]<-'ecosystem2'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_7','MP_17')]<-'ecosystem3'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_1','MP_4')]<-'ecosystem4'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
eco_MP$ecosystem[eco_MP$top_MP%in%c('MP_15','MP_11')]<-'ecosystem6'
table(eco_MP$top_MP,eco_MP$ecosystem)
integration_umap_MP$ecosystem<-eco_MP$ecosystem

all_mp<-read.delim('/10X_Visium/all_umap_order.txt',stringsAsFactors = F,check.names = F)
all_mp[1:3,]
table(rownames(integration_umap_MP)==rownames(all_mp))

integration_umap_MP$FinalLocalType<-all_mp$FinalLocalType

write.table(integration_umap_MP,'/10X_Visium/new_st_NMF/MPscore_new/all_MP_ecoSystem_6.txt',
            quote = F,sep = '\t')


dir_EMT<-'/10X_Visium/new_st_18Imm/'
file_Imm<-list.files(pattern = '_18Imm.txt',path = dir_EMT,recursive = T)
length(file_Imm)
file_Imm[1:10]
dataSlice<-unlist(lapply(strsplit(file_Imm,'_18I'),function(x)x[1]))
###数据来源360_newST_18Imm.R
dataSlice<-gsub('/','_',dataSlice)
dataSlice[1:10]


all_FC_re<-c()
all_P_re<-c()
for(i in 1:length(file_Imm)){
  #i=1
  Imm_data<-read.delim(paste0(dir_EMT,file_Imm[i]),stringsAsFactors = F,check.names = F)
  Imm_data[1:3,1:3]
  dim(Imm_data)
  slice_MP<-integration_umap_MP[integration_umap_MP$dataSlice%in%dataSlice[i],]
  rownames(slice_MP)<-slice_MP$cell_name
  slice_MP<-slice_MP[rownames(Imm_data),]
  slice_MP<-slice_MP[,c(1:17)]
  slice_MP$MP_top1<-apply(slice_MP,1,function(x){
    names(x)[which.max(x)]
  })
  table(slice_MP$MP_top1)
  slice_MP[1:3,]
  dim(slice_MP)
  
  Imm_data$top_MP<-slice_MP$MP_top1
  # Imm_data<-Imm_data[which(Imm_data$LocalType!='not.defined'),]
  Imm_data$ecosystem<-''
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_8')]<-'ecosystem1'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_6','MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem2'
  # Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_3','MP_12','MP_13','MP_14','MP_16')]<-'ecosystem3'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_7','MP_17')]<-'ecosystem3'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_1','MP_4')]<-'ecosystem4'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_9','MP_10','MP_2','MP_5')]<-'ecosystem5'
  Imm_data$ecosystem[Imm_data$top_MP%in%c('MP_15','MP_11')]<-'ecosystem6'
  #table(Imm_data$top_MP,Imm_data$ecosystem)
  
  all_DE_Imm<-data.frame(ImmPath=colnames(Imm_data)[1:18])
  for(j in unique(Imm_data$ecosystem)){
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
  FC_re<-reshape2::melt(as.matrix(FC_re))
  P_re<-reshape2::melt(as.matrix(P_re))
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
write.table(all_P_re,'/10X_Visium/new_st_NMF/eco_DE_Imm/ecosystem_DE_ImmPath_6.txt',quote = F,sep = '\t',row.names = F)

all_P_re<-read.delim('/10X Visium/new_st_NMF/eco_DE_Imm/ecosystem_DE_ImmPath_6.txt',stringsAsFactors = F,check.names = F)


####将所有ecosystem和一起画
col_data<-data.frame(ImmPath=unique(all_P_re$ImmPath),
                     col=c('#a4cde1','#277fb8','#96cb8f','#4dae47','#b79973','#f38989','#e32427','#003366','#f9b769',
                           '#d4a6a8','#815e99','#b05a28','#663399','#FFCC00','#EEF0A5','#336633','#FC733D','#CCCC99'))

dfcol<-data.frame(ImmPath=paste0('ecosystem',1:6),
                  logFC=0,
                  label=rep('',6))
mycol <- c('#4b2d80','#136599','#1c959d','#8e2823','#d1e181','#FBEB6B')#,'#AEEBE3'

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
pdf('/10X Visium/new_st_NMF/eco_DE_Imm/00pic/ecosystem_path_point_bar_6.pdf',width = 12,height = 6)
print(p3)
dev.off()













