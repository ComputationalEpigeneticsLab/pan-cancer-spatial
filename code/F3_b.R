######MP与cancer program的关系
###桑基图展示
library(ggalluvial)
library(jsonlite)
library(plyr)
library(ggplot2)



new_MP<-readRDS('/10X Visium/new_st_NMF/intra23_inter23Cluster_list(20-300).rds')
new_MP2<-unlist(new_MP)
new_MP2<-data.frame(cluster=rep(names(new_MP),as.numeric(unlist(lapply(new_MP,length)))),
                    slice=new_MP2)
table(new_MP2$cluster)


plot_data<-new_MP2
plot_data$cluster<-gsub('Cluster','MP',plot_data$cluster)
plot_data$cancer<-unlist(lapply(strsplit(plot_data$slice,'_'),function(x)x[1]))
table(plot_data$cancer)
length(unique(plot_data$cancer))
plot_data$program<-paste0(plot_data$cancer,'_program')
colnames(plot_data)[1]<-'MP'


plot_data2 <- plot_data[,c('MP','program')]
plot_data2 <- as.matrix(plot_data2)
dataP <- data.frame(frenq = 1,
                    Cohort = rep(c(1:nrow(plot_data2)),times = 2),
                    x = rep(c("MP","program"),each = nrow(plot_data2)),
                    stratum = c(plot_data2[,1],plot_data2[,2])) 
dataP$x <- factor(dataP$x,levels = c("MP","program"))
dataP$stratum <- factor(dataP$stratum,levels = c(unique(plot_data2[,1]),unique(plot_data2[,2])))
unique(dataP$stratum)

cancer_col<-read.delim('/35cancer_color.txt')
cancer_col<-cancer_col[cancer_col$cancer%in%unique(plot_data$cancer),]
col2<-cancer_col$color
names(col2)<-paste0(cancer_col$cancer,'_program')

color3<-c(c("MP_1"="#fb6a4b","MP_2"="#fe9376","MP_3"="#008B8B","MP_4"="#41b9C1","MP_5"="#6A8EC9",
            "MP_6"="#817cb9","MP_7"="#cb78a6","MP_8"="#c65861","MP_9"="#652884","MP_10"="#444577",
            "MP_11"="#8A7355","MP_12"="#B3BB61","MP_13"="#9d5c39","MP_14"="#fcb93e","MP_15"="#FFB978",
            "MP_16"="#399335",'MP_17'='#96DD88'),
          col2)

#pdf("riverplot_23_2_28.pdf",width = 5)
p_river<-ggplot(dataP,
                aes(x=x,y=frenq,stratum=stratum,alluvium=Cohort,fill=stratum,label=stratum)) +
  geom_flow(width=1/9) +
  geom_stratum(width=1/9,linetype=1.5,size=0.5,alpha=1,color="white") +
  geom_text(stat="stratum",size=3,nudge_x=0.2) +
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
#dev.off()
##5*8
print(p_river)
pdf(paste0('/10X Visium/new_st_NMF/','MP_program_river.pdf'),width = 5,height = 7)
print(p_river)
dev.off()










