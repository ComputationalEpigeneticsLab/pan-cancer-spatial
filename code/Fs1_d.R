library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(aplot)
library(openxlsx)
library(ggpubr)
library(patchwork)

setwd('/ST/')
location_data = read.table('ST_merge_res0.8_metadatainfor_BdyCoreBud.txt',sep = '\t',row.names = 1,header = T);gc()
location_data[location_data$FinalLocalType%in%c('Budding'),'FinalLocalType']='Dispersion'

cancer_name <- unique(location_data$cancer)
location_info_all <- data.frame()
for (i in 1:length(unique(location_data$cancer))) {
  #i <- 1
  cancer0 <- cancer_name[i]
  cancer_data <- location_data[location_data$cancer==cancer0,]
  location_info <- as.data.frame(table(cancer_data$FinalLocalType))
  location_info$cancer <- cancer0
  #location_info<-location_info[which(location_info$Var1!='not.defined'),]
  location_info$propotion <- location_info$Freq/sum(location_info$Freq)
  location_info <- location_info %>%
    arrange(match(Var1, c("Core","Boundary","Dispersion","Immune","Normal"))) %>% as.data.frame()
  
  location_info_all <- rbind(location_info_all,location_info)
}
order = location_info_all[location_info_all$Var1%in%'Core' & location_info_all$cancer != 'all',] %>% arrange(desc(propotion))
cancer_order = c('all',order$cancer)

location_info_all_info<- as.data.frame(table(location_data$FinalLocalType))
location_info_all_info$cancer <- "all"
location_info_all_info$propotion <- location_info_all_info$Freq/sum(location_info_all_info$Freq)
location_info_all <- rbind(location_info_all,location_info_all_info)

location_info_all <- mutate(location_info_all,
                            Var1 = factor(location_info_all$Var1,
                                          levels = rev(c("Core","Boundary","Dispersion","Immune","Normal"))))

cancer_types <- c('all',"GBM", "HNSCC", "HN-AS", "DIPG", "PN", "PCNSL", "DLBCL", "LGACC", "OSCC", 
                  "NPC", "TC", "DSRCT", "PTCL", "BRCA", "IPMN", "LUSC", "LUAD", "LNC", 
                  "LAM", "LIHC", "HB", "GC", "GIST", "CRC", "RCC", "PDAC", "EC", "OVCA", 
                  "HGSC", "CESC", "PRAD", "MIBC", "OS", "CSCC", "SKCM")
cancer_colors <- c('#97989B',"#915824", "#AD8031", "#C9A063", "#C5A97F", "#CBBE82", "#BAB64D", "#C7C458", 
                   "#907F28", "#D8E7A3", "#E4E653", "#D79E63", "#EDC824", "#DF7561", "#91302C", 
                   "#B66C6A", "#C693C1", "#C7A1C9", "#B8A3B9", "#B8AECA", "#55A7AE", "#89CDC9", 
                   "#649DD4", "#C1E4F8", "#D79E63", "#81B464", "#506248", "#BF5099", "#D89AB1", 
                   "#F5B8C9", "#F5B2B1", "#488B45", "#EA9D19", "#9F8A83", "#C4CDB5", "#54626A")
names(cancer_colors)=cancer_types

location_info_all <- mutate(location_info_all,
                            cancer = factor(location_info_all$cancer,
                                            levels = cancer_order))

region_types <- c("Core", "Boundary", "Dispersion", "Stormal", "Immune", "Normal")
region_colors <- c("#ca3028", "#eeb46c", "#e37330", "#a0d5e8", "#a4cc87", "#5272a6")
names(region_colors)=region_types

p1 <- ggplot(location_info_all,aes(x=cancer,y=propotion,fill=Var1))+
  geom_col(position = "fill", width = 0.6)+
  scale_fill_manual(values =  region_colors)+
  labs(x = '', fill = "",title = "")+
  theme_test()+
  theme(plot.title = element_text(size = 10,colour = "black",face="bold",hjust = 0.5),
        axis.title.y = element_text(size = 9,color = "black",face = "bold", vjust = 0.5, hjust = 0.5,angle = 90),
        legend.title = element_text(color="black",size=9,face="bold"),
        legend.text = element_text(color="black",size = 8),
        #axis.text.x = element_blank(),
        axis.text.x = element_text(size = 8,color = "black",vjust = 0.5,hjust = 1,angle = 90),
        axis.text.y = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5, angle = 0))
p1

cancer = as.data.frame(cancer_order)
Cancer <- cancer %>% 
  ggplot(aes(x = cancer_types, y = "",,fill=cancer_types))+
  scale_fill_manual(values = cancer_colors)+
  scale_x_discrete(position="bottom")+
  geom_tile() + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8,color = "black",vjust = 0.5,hjust = 1,angle = 90),
    axis.text.y = element_blank(), 
    panel.grid = element_blank(),
    plot.margin = margin(0, 0, 0, 0),
    legend.position = 'none'
  ) +
  labs(x = NULL, y = NULL, fill = "",title = "")
plot=p1 %>%
  insert_bottom(Cancer,height=.03)
plot

pdf('Fs1d.pdf',width = 15,height = 6)
print(plot)
dev.off()
