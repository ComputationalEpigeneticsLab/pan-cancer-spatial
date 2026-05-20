##bdy，hop1代谢差异####
library(ggridges)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

metabolic_df = read.delim('/scFEA/selectscFEA_flux.txt',check.names = FALSE);gc()
bdylist_uni = read.delim('/metabolism/bdy_hop10_unique.txt',check.names = FALSE);gc()
loc = read.delim('/ST/ST_merge_res0.8_metadatainfor_BdyCoreBud.txt',check.names = FALSE);gc()

bdy = metabolic_df[rownames(loc)[which(loc$FinalLocalType == "Boundary")],];gc()
hop1 = metabolic_df[rownames(bdylist_uni)[which(bdylist_uni$step == "step_1")],];gc()
proximal = metabolic_df[rownames(bdylist_uni)[which(bdylist_uni$step%in%paste0('step_',1:3))],];gc()
distal = metabolic_df[rownames(bdylist_uni)[which(bdylist_uni$step%in%paste0('step_',4:10))],];gc()

bdy$type='bdy'
hop1$type='hop1'
proximal$type='proximal'
distal$type='distal'

plotdata = rbind(bdy,hop1,proximal,distal)
plotdata$all <- rowMeans(plotdata[, 1:70], na.rm = TRUE)

select_md=read.csv('/scFEA/ref/Human_M168_energy_anabolism_filtered.csv')
select_md$name = apply(select_md,1,function(x){paste0(x[3]," -> ", x[6])})
identical(colnames(plotdata)[1:70],select_md[,1])

pdf(paste0("/metabolism/bdy&hop1&near&far代谢活性差异(所有通路).pdf"),width = 4.5, height = 4)
for (i in 1:70) {
  pdata = plotdata[,c(i,72)]
  colnames(pdata)[1]='score'
  pdata[which(is.na(pdata$score)),'score']=0
  threshold <- as.numeric(quantile(pdata$score, 0.99, na.rm = TRUE))
  pdata <- pdata[pdata$score <= threshold, ]
  
  bdy_hop1 = wilcox.test(pdata[pdata$type%in%'bdy',]$score,pdata[pdata$type%in%'hop1',]$score)
  near_far = wilcox.test(pdata[pdata$type%in%'proximal',]$score,pdata[pdata$type%in%'distal',]$score)
  pdata$type = factor(pdata$type,levels = c('bdy','hop1','proximal','distal'))
  p=ggplot(pdata, aes(x = score, y=type, fill = type)) +
    geom_density_ridges(alpha = 0.7,color='white',scale = 1.5)+
    theme_bw()+
    labs(x='all', y="",title = '') +
    annotate("text", x = (max(pdata$score)-0.01), y = 3.5,color = "black", size = 3,
             label = ifelse(bdy_hop1$p.value > 2e-16,
                      paste0("p = ", format.pval(bdy_hop1$p.value, digits = 3)),
                      paste0("p ", format.pval(bdy_hop1$p.value, digits = 3)))) +
    annotate("text", x = (max(pdata$score)-0.01), y = 1.5,color = "black", size = 3,
             label = ifelse(near_far$p.value > 2e-16,
                            paste0("p = ", format.pval(near_far$p.value, digits = 3)),
                            paste0("p ", format.pval(near_far$p.value, digits = 3)))) +
    scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))+
    scale_fill_manual(values = c("proximal"="#827CBA","distal"="#BBB7D8","bdy"="#eeb46c","hop1"="#a0d5e8"))+
    theme(plot.title = element_text(size = 10,colour = "black",face="bold",hjust = 0.5),
          axis.title.x = element_text(size = 9,color = "black",face = "bold", vjust = 0.5, hjust = 0.5,angle = 0),
          legend.title = element_text(color="black",size=9,face="bold"),
          legend.text = element_text(color="black",size = 8),
          legend.position = 'none',
          axis.text.x = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0),
          axis.text.y = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0))
  print(p)
}
dev.off()