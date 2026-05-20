library(ggridges)
library(openxlsx)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)

bdylist_uni = readRDS('E:/panCancerST/metabolism/bdy_hop10_unique.rds');gc()
RCTD = readRDS('E:/panCancerST/Deconvolution/RCTD_result.rds')
RCTD_df <- do.call(rbind, RCTD)

metabolic_df = readRDS('E:/panCancerST/scFEA/scFEA_flux_selectmean.rds');gc()
metabolic_df = as.data.frame(metabolic_df)
rownames(metabolic_df) = lapply(rownames(metabolic_df),function(x){gsub('\\.','_',x)}) %>% unlist()
metabolic_df = metabolic_df %>% rownames_to_column()

plots = list()
for (i in 1:10) {
  aa = metabolic_df[metabolic_df$rowname%in%intersect(bdylist_uni[[i]],metabolic_df$rowname),]
  aa$step = paste0('step_',i)
  plots[[i]] = aa
}
plotdata = do.call(rbind,plots)
#plotdata = plotdata[-c(which(is.na(plotdata$metabolic_df))),]
plotdata$type = ifelse(plotdata$step%in%paste0('step_',1:3),'proximal','distal')
#plotdata$step = factor(plotdata$step,levels = paste0('step_',1:10))

TAM_cell = RCTD_df[RCTD_df$celltype%in%c('TAM'),'rowname']
CAF_cell = RCTD_df[RCTD_df$celltype%in%c('CAF'),'rowname']
TAM = plotdata[plotdata$rowname%in%TAM_cell,]
CAF = plotdata[plotdata$rowname%in%CAF_cell,]
TAM$celltype='TAM'
CAF$celltype='CAF'
plotdata$celltype='all'
pdata = rbind(plotdata,TAM,CAF)
pdata$cell_loc = apply(pdata,1,function(x){paste0(x[5],"_",x[4])})

pdata = pdata[-c(which(pdata$metabolic_df>0.015)),]

##近远端
TAM_test = wilcox.test(pdata[pdata$cell_loc%in%'TAM_proximal',]$metabolic_df,pdata[pdata$cell_loc%in%'TAM_distal',]$metabolic_df)
CAF_test = wilcox.test(pdata[pdata$cell_loc%in%'CAF_proximal',]$metabolic_df,pdata[pdata$cell_loc%in%'CAF_distal',]$metabolic_df)
all_test = wilcox.test(pdata[pdata$cell_loc%in%'all_proximal',]$metabolic_df,pdata[pdata$cell_loc%in%'all_distal',]$metabolic_df)

p=ggplot(pdata, aes(x = metabolic_df, y=cell_loc, fill = type)) +
  geom_density_ridges(alpha = 0.7,color='white',scale = 1.5)+
  theme_bw()+
  labs(x="scFEA", y="",title = '') +
  annotate("text", x = 0.01, y = 5.5,label = paste0("p =", format.pval(TAM_test$p.value, digits = 3)),color = "black", size = 3) +
  annotate("text", x = 0.01, y = 3.5,label = paste0("p ", format.pval(CAF_test$p.value, digits = 3)),color = "black", size = 3) +
  annotate("text", x = 0.01, y = 1.5,label = paste0("p ", format.pval(all_test$p.value, digits = 3)),color = "black", size = 3) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))+
  scale_fill_manual(values = c("proximal"="#827CBA","distal"="#BBB7D8"))+
  theme(plot.title = element_text(size = 10,colour = "black",face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 9,color = "black",face = "bold", vjust = 0.5, hjust = 0.5,angle = 0),
        legend.title = element_text(color="black",size=9,face="bold"),
        legend.text = element_text(color="black",size = 8),
        legend.position = 'none',
        axis.text.x = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0))
p
ggsave("all-score_近远端代谢活性差异(部分代谢通路).pdf",p,device = "pdf",
       path = "E:/panCancerST/metabolism/score",height = 3,width = 3.3)