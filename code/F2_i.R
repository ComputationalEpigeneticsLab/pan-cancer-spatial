################################################################################################
####免疫配体表达热图##############################################

all_immLigand<-read.delim(paste0(dir_out,'all_immLigand_CoreBdyBud.txt'),stringsAsFactors = F,check.names = F)
apply(all_immLigand[,5:17],2,function(x){
  length(which(x==0))
})
all_immLigand<-all_immLigand[,-18]

cancer<-unique(all_immLigand$cancer)

all_exp_mean<-c()
all_p_value<-c()
library(tidyverse)

for(i in 1:length(cancer)){#i=1
  cancer_Imm<-all_immLigand[which(all_immLigand$cancer==cancer[i]),]
  exp_mean<-aggregate(cancer_Imm[,5:17],by=list(cancer_Imm$LocalType),mean)
  exp_mean<-exp_mean[which(exp_mean$Group.1=='Core'),]
  rownames(exp_mean)<-cancer[i]
  exp_mean<-exp_mean[,-1]
  all_exp_mean<-rbind(all_exp_mean,exp_mean)
  
  core_site<-which(cancer_Imm$LocalType=='Core')
  BdyDis_site<-which(cancer_Imm$LocalType=='Boundary'|cancer_Imm$LocalType=='Budding')
  
  p_value<-apply(cancer_Imm[,5:17],2,function(x){
    wilcox.test(x[core_site],x[BdyDis_site],alternative = 'greater')[["p.value"]]
  }) %>% t() %>% as.data.frame()
  rownames(p_value)<-cancer[i]
  all_p_value<-rbind(all_p_value,p_value)
}

all_exp_mean<-t(all_exp_mean)
all_p_value<-t(all_p_value)

dis_dot<-ifelse(as.matrix(all_p_value) < 0.05, "*", "")
# dis_dot<-ifelse(as.matrix(all_p_value) < 0.01, "**", "")
# dis_dot<-ifelse(as.matrix(all_p_value) < 0.001, "***", "")
#diag(dis_dot)<-""
dis_dot[is.na(dis_dot)]<-""
range(all_exp_mean,na.rm = T)
heat_data<-scale(all_exp_mean)
heat_data<-apply(all_exp_mean, 1,function(x){
  x<-x-min(x)
  x<-x/max(x)
}) %>% t()
range(heat_data)
bk<-seq(0,1,length.out=100)
color_pheatmap<-c(colorRampPalette(c("white",'#E9C1C6'))(5),
                  colorRampPalette(c("#E9C1C6",'#b5414c'))(95)) ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(heat_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = F,
                      cluster_cols = T,
                      treeheight_row = T,treeheight_col = T,
                      display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "white",###black
                      fontsize = 10,
                      cellwidth=15,
                      cellheight=15,
                      main = "Imm_ligand_CoreVsBdyDis",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('D:/pan_cancer/0_修稿/Fs2/ImmLigand_CoreVsBdyDis.pdf',width = 5.5,height = 4)
print(p)
dev.off()


