######MP与cancer program的关系
###桑基图展示
library(ggalluvial)
library(jsonlite)
library(plyr)
library(ggplot2)





program_MP<-read.delim('E:/Mirror/ST_analysis/data/10X Visium/NMF_module/program_top50gene.txt',
                       stringsAsFactors = F,check.names = F)
plot_data<-program_MP
plot_data$cluster<-gsub('Cluster','MP',plot_data$cluster)
plot_data<-plot_data[which(plot_data$cluster!='MP_1'),]
plot_data$cancer<-unlist(lapply(strsplit(plot_data$slice,'_'),function(x)x[1]))
plot_data$cancer<-substr(plot_data$cancer,1,nchar(plot_data$cancer)-2) %>% toupper()
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

color3<-c(c("MP_2"="#1F77B2","MP_3"="#FF7F0E","MP_4"="#279C68","MP_5"="#D42728",
            "MP_6"="#A840FA","MP_7"="#8A564B","MP_8"="#E177C0","MP_9"="#B3BB61",
            "MP_10"="#17BCCD","MP_11"="#ACC5E6","MP_12"="#FFB978","MP_13"="#96DD88",
            'MP_14'='#FF9694','unresolved'='#7F7F7F'),
          c('BRCA_program'='#b7996d','CESC_program'='#e32427','CRC_program'='#8bc96d','CSCC_program'='#b05a28','GBM_program'='#a4cde1',
            'GIST_program'='#96cb8f','HGSC_program'='#277fb8','HN-AS_program'='#f38989','IPMN_program'='#5c9e43','LIHC_program'='#c6b598',
            'LUAD_program'='#7A9AC5','MIBC_program'='#60592E','OSCC_program'='#C5BE97','OVCA_program'='#C89192','PCNSL_program'='#44637F',
            'PDAC_program'='#549da3','PRAD_program'='#f9b769','RCC_program'='#af93c4','SKCM_program'='#d4a55b'))

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
pdf(paste0('E:/Mirror/ST_analysis/pic/re/MP_score/','MP_program_river.pdf'),width = 5,height = 6)
print(p_river)
dev.off()










