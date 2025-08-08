####MP14基因交集热图
library(Seurat)
library(rlang)
library(ggplot2)
#library(tidyverse)
library(ggraph)
library(ggpubr)
library(dplyr)
library(NMF)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(scales)




Program_713_org<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/program_top50gene.txt',
                            stringsAsFactors = F,check.names = F)
Program_713_org$cluster<-gsub("Cluster","MP",Program_713_org$cluster)
table(Program_713_org$cluster)
Program_713<-Program_713_org

order_MP<-c('MP_3','MP_2','MP_5','MP_4','MP_7','MP_6','MP_9','MP_11','MP_10','MP_8','MP_13','MP_12','MP_14')

Program_713<-Program_713[lapply(order_MP,function(x) which(Program_713$cluster==x))%>%unlist(),]
#lapply(order_MP,function(x) which(Program_713$cluster==x))%>%unlist()

Program_713_gene<-lapply(Program_713$gene,function(x) unlist(strsplit(x,',')))
Program_713_gene<-do.call(cbind,Program_713_gene)
colnames(Program_713_gene)<-Program_713$slice
plot_data <- apply(Program_713_gene , 2, function(x) apply(Program_713_gene , 2, function(y) length(intersect(x,y)))) 


plot_data<-apply(plot_data,2,function(x){
  100*x/(100-x)
})

###每个MP内部排序
all_order<-c()
for(MP_i in order_MP){#MP_i=order_MP[1]
  order_data<-plot_data[Program_713$slice[Program_713$cluster%in%MP_i],Program_713$slice[Program_713$cluster%in%MP_i]]
  program_order<-apply(order_data,1,function(x){sum(x,na.rm = T)})
  program_order<-program_order[order(program_order,decreasing = T)]
  all_order<-c(all_order,names(program_order))
}

plot_data<-plot_data[all_order,all_order]
Program_713<-Program_713[match(all_order,Program_713$slice),]

la <- rowAnnotation(df = data.frame(MP_list=Program_713$cluster),
                    col = list(MP_list=c(#"MP_1"="#A7CBEF",
                      "MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
                      "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
                      "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
                      'MP_14'='#FF9694'
                    )
                    ))
# "MP_8"="#99CC66","MP_9"="#99CC66",
# "MP_2"="#5F91C7",
# "MP_3"="#EFA7A9","MP_6"="#EFA7A9","MP_7"="#EFA7A9","MP_14"="#EFA7A9",
# "MP_4"="#F9D977",
# "MP_5"="#ECBEAA","MP_12"="#ECBEAA","MP_13"="#ECBEAA",
# "MP_10"="#B196C1",
# "MP_11"="#7CBFB6"
range(plot_data)
quantile(plot_data)
col_fun <- colorRamp2(
  c(0, 50, 100), 
  c("#1f78b4", "white","#ff7f00")
)
col_fun <- colorRamp2(
  c(0, 5, 15, 25, 100), 
  c("white","#FBF4C3",'#FBF4C3','#CC4476','#251F49')
)
col_fun <- colorRamp2(
  c(0, 5, 15, 20, 100), 
  c("white","#FBF4C3",'#FBF4C3','#CC4476','#251F49')
)
# col_fun <- colorRamp2(
#   c(0, 5, 15, 20, 100), 
#   c("white","#FEF6C5",'#FBEE92','#CC4476','#251F49')
# )
p1<-Heatmap(plot_data,
            col = col_fun,
            #top_annotation =ha,#####顶部注释
            left_annotation=la,#####左侧注释
            show_row_names = F,#####不显示行名
            show_column_names =F,#####不显示列名
            cluster_rows =F,
            cluster_columns = F,
            #name="TCGA_Immune_infiltration",
            column_title='MP',
            row_names_gp = gpar(fontsize = 10))
print(p1)

pdf('E:/Mirror/ST_analysis/pic/re/5/core_MP13_intersect_order.pdf',height = 8,width = 8)
print(p1)
dev.off()

# pdf('E:/Mirror/ST_analysis/pic/NMF_module/core_re/core_MP14_intersect.pdf',height = 8,width = 8)
# print(p1)
# dev.off()
