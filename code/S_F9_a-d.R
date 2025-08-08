###cellchat分开统计
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(reshape2)


####
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_CAFstep1/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_Bdy_TAMstep1/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_TAM/'
dir_cellchat<-'E:/Mirror/ST_analysis/data/10X Visium/cellchat_matrix_CAF/'
file_LR<-list.files(pattern = 'net_lr.csv',path = dir_cellchat,recursive = T)
dataSlice<-unlist(lapply(strsplit(file_LR,'_'),function(x)x[1]))

all_path_LR<-c()
for(i in 1:length(file_LR)){
  #i=1
  LR_data<-read.csv(paste0(dir_cellchat,file_LR[i]),stringsAsFactors = F,check.names = F,row.names = 1)
  LR_data$slice<-dataSlice[i]
  colnames(LR_data)
  LR_data<-LR_data[which(LR_data$pval<0.05),]
  LR_data<-LR_data[,c('source','target','pathway_name','ligand','receptor','annotation','evidence','slice')]
  
  need_site<-lapply(1:nrow(LR_data),function(x){#x=1
    a<-c(LR_data$source[x],LR_data$target[x])
    other<-lapply(setdiff(unique(c(LR_data$source,LR_data$target)),'CAF'),function(y){
      return(y%in%a)
    }) %>% unlist()
    return('CAF'%in%a&&TRUE%in%other)
  }) %>% unlist()
  table(need_site)
  LR_data<-LR_data[need_site,]
  if(nrow(LR_data)>0){
    all_path_LR<-rbind(all_path_LR,LR_data)
  }
  print(dataSlice[i])
}
all_path_LR$pathway_name<-as.vector(all_path_LR$pathway_name)

write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_diviLR.txt',quote = F,sep = '\t',row.names = F)
write.table(all_path_LR,'E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_diviLR.txt',quote = F,sep = '\t',row.names = F)

# bdy_step1CAF<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF.txt',stringsAsFactors = F,check.names = F)
# sum(bdy_step1CAF$LR_num)
# 
# class(all_path_LR$pathway_name)
# a_all_path_LR<-all_path_LR[all_path_LR$slice%in%'brca01/slice1',]
# b_bdy_step1CAF<-bdy_step1CAF[bdy_step1CAF$slice%in%'brca01/slice1',]
#############################################################################要注意 不同细胞对之间可能有相同的通路和LR对
#########################################################################不能用table计算频率 重复出现的LR会导致pathway少

bdy_step1CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1CAF_diviLR.txt',stringsAsFactors=F,check.names=F)
bdy_step1TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/bdy_step1TAM_diviLR.txt',stringsAsFactors=F,check.names=F)

matrix_CAF_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_CAF_diviLR.txt',stringsAsFactors = F,check.names = F)
matrix_TAM_divi<-read.delim('E:/Mirror/ST_analysis/data/pic_data/cellchat/matrix_TAM_diviLR.txt',stringsAsFactors = F,check.names = F)
#sum(matrix_CAF_divi$LR_num)

bdy_step1CAF_divi$L_R<-'source'
bdy_step1CAF_divi$L_R[which(bdy_step1CAF_divi$target=='CAF')]<-'target'
bdy_CAFpath<-as.data.frame.array(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$L_R))
bdy_CAFpath$sum<-bdy_CAFpath$source+bdy_CAFpath$target

bdy_step1TAM_divi$L_R<-'source'
bdy_step1TAM_divi$L_R[which(bdy_step1TAM_divi$target=='TAM')]<-'target'
bdy_TAMpath<-as.data.frame.array(table(bdy_step1TAM_divi$pathway_name,bdy_step1TAM_divi$L_R))
bdy_TAMpath$sum<-bdy_TAMpath$source+bdy_TAMpath$target

matrix_CAF_divi$L_R<-'source'
matrix_CAF_divi$L_R[which(matrix_CAF_divi$target=='CAF')]<-'target'
matrix_CAFpath<-as.data.frame.array(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$L_R))
matrix_CAFpath$sum<-matrix_CAFpath$source+matrix_CAFpath$target

matrix_TAM_divi$L_R<-'source'
matrix_TAM_divi$L_R[which(matrix_TAM_divi$target=='TAM')]<-'target'
matrix_TAMpath<-as.data.frame.array(table(matrix_TAM_divi$pathway_name,matrix_TAM_divi$L_R))
matrix_TAMpath$sum<-matrix_TAMpath$source+matrix_TAMpath$target


####################################################################################
###river plot####################################################################################
cell_col<-c("B lymphocytes"="#66CC00",'CAF'='#A78E41',"Endothelial"="#CCCC99","Epithelial"="#93C647",
            "Macrophage"="#CC3333","Myeloid cell"="#ED703F","NK cell"="#D2AF83","T lymphocytes"="#F3A383",
            'TAM'='#7B3257',"Fibroblasts"="#8B964F","MAST cell"="#FF9900","Monocyte"="#EFA7A9",
            "Neutrophils"="#EDDC6D","Dendritic"="#FFFF00",'unknown'='grey85',
            
            "GC B cells in the DZ"='#CC9933',
            "Plasma cells"='#FFCCCC',"CD8+ T Memory"='#996699',"follicular B cells"='#CCCC99',
            "Treg"='#FFCC33',"Cytotoxic"='#FF6666',
            "TAM_C0"='#336699',"TAM_C1"='#99CCCC',"TAM_C2"='#CCFFFF',"TAM_C3"='#99CC33',
            "Naive"='#FF9900',"B cell Regulatory"='#990033',"Naive B cell"='#990066',
            'Core'='#d62d28','bdy'='#f6b86d','Dispersion'='#ee762d')

library(ggalluvial)
library(jsonlite)
library(plyr)
library(ggplot2)


bdy_step1CAF_divi
bdy_CAFpath

bdy_step1TAM_divi
bdy_TAMpath

matrix_CAF_divi
matrix_CAFpath

matrix_TAM_divi
matrix_TAMpath

river_data<-bdy_step1CAF_divi[bdy_step1CAF_divi$pathway_name%in%rownames(bdy_CAFpath)[order(bdy_CAFpath$sum,decreasing = T)][1:10],]

data_del <- river_data
data_del$cell_pairs <- paste(data_del[,"source"],data_del[,"target"],sep='_')
data_del$LR <- paste(data_del[,"ligand"],data_del[,"receptor"],sep='_of_')
# gene_count <- ddply(data_del,.(ligand,receptor),function(x) length(x$cell_pairs))
#gene_count_sort <- gene_count[order(gene_count[,3],decreasing = T),,drop=F]
#gene_count<-aggregate(data_del$LR,by=list(data_del$cell_pairs,data_del$pathway_name),function(x){length(x)})
gene_count<-as.data.frame.array(table(data_del$LR,data_del$pathway_name))
gene_count_sort<-lapply(1:ncol(gene_count),function(x){#x=1
  a<-gene_count[,x]
  names(a)<-rownames(gene_count)
  a<-a[order(a,decreasing = T)]
  return(names(a)[1:5])
}) %>% unlist() %>% unique()

gene_plot <- data.frame(ligand=unlist(lapply(strsplit(gene_count_sort,'_of_'),function(x)x[1])),
                        receptor=unlist(lapply(strsplit(gene_count_sort,'_of_'),function(x)x[2])))
plot_data <- data_del[data_del$ligand %in% gene_plot[,1]&data_del$receptor %in% gene_plot[,2],
                      c('source',"ligand","pathway_name","receptor","target")]
rownames(plot_data)<-1:nrow(plot_data)
plot_data<-plot_data[c(which(plot_data$source=="CAF"),which(plot_data$source!="CAF")),]
plot_data<-plot_data[c(which(plot_data$source=="TAM"),which(plot_data$source!="TAM")),]

#plot_data2 <- plot_data[order(plot_data$cell_from),]
plot_data2 <- plot_data
plot_data2$pathway_name <- paste0('FUN_',plot_data2$pathway_name) 
plot_data2 <- as.matrix(plot_data2)
dataP <- data.frame(frenq = 1,
                    Cohort = rep(c(1:nrow(plot_data2)),times = 5),
                    x = rep(c('source',"ligand","pathway_name","receptor","target"),each = nrow(plot_data2)),
                    stratum = c(plot_data2[,1],plot_data2[,2],plot_data2[,3],plot_data2[,4],plot_data2[,5])) 
dataP$x <- factor(dataP$x,levels = c('source',"ligand","pathway_name","receptor","target"))
dataP$stratum <- factor(dataP$stratum,levels = c(unique(plot_data2[,1]),union(plot_data2[,2],plot_data2[,4]),unique(plot_data2[,3])))

setdiff(c(unique(plot_data2[,1]),union(plot_data2[,2],plot_data2[,4]),unique(plot_data2[,3])),unique(dataP$stratum))
setdiff(unique(dataP$stratum),c(unique(plot_data2[,1]),union(plot_data2[,2],plot_data2[,4]),unique(plot_data2[,3])))
cellType<-unique(plot_data2[,1])
color<-cell_col[cellType]
gene <- union(plot_data2[,2],plot_data2[,4])
func <- unique(plot_data2[,3])
#color3 <- c(color[1:9,2],colorRampPalette(c("#F28E2B","#8CD17D","#4E79A7"))(52),c("#FABFD2","#F1CE63","#86BCB6","#FF9D9A"))
color3 <- c(color,
            colorRampPalette(c("#F2A75D","#8CD17D","#6485A7"))(length(gene)),
            colorRampPalette(c("#FFCC99","#FF9966"))(length(func)))
names(color3) <- c(cellType,gene,func)

p_river<-ggplot(dataP,
                aes(x=x,y=frenq,stratum=stratum,alluvium=Cohort,fill=stratum,label=stratum)) +
  geom_flow(width=1/9) +
  geom_stratum(width=1/9,linetype=1.5,size=0.5,alpha=1,color="white") +
  geom_text(stat="stratum",size=2.5,nudge_x=0.2) +
  scale_x_discrete(limits=c()) +
  theme_bw()+
  theme(legend.position="none",
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank()) +
  scale_fill_manual(values=color3)
print(p_river)
pdf(paste0(dir_pic,'river_matrix_TAM.pdf'),width = 9,height = 7)##饼图 width = 6,height = 5
print(p_river)
dev.off()





