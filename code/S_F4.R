####MP挑选通路
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


###数据来源 
dir_MP<-'E:/Mirror/ST_analysis/data/10X Visium/NMF_module/Function-geneset/'
file_MP<-list.files(pattern = 'txt',path = dir_MP,recursive = T)
MP_name<-unlist(lapply(strsplit(file_MP,'_pValue'),function(x)x[1]))
MP_order<-c('MP_11','MP_8','MP_10','MP_2','MP_14','MP_3','MP_12','MP_9','MP_6','MP_7','MP_4','MP_5','MP_13')
#names(file_MP)<-MP_name
#file_MP<-file_MP[c('MP11','MP8','MP10','MP2','MP14','MP3','MP12','MP9','MP6','MP7','MP4','MP5','MP13')]

all_MP_data<-c()
for(i in 1:length(file_MP)){
  #i=12
  MP_data<-read.delim(paste0(dir_MP,file_MP[i]),stringsAsFactors = F,check.names = F)
  # MP_data<-MP_data[!duplicated(MP_data$pathway),]
  # rownames(MP_data)<-MP_data$pathway
  # MP_data<-MP_data[,-1]
  # write.table(MP_data,'E:/Mirror/ST_analysis/data/10X Visium/NMF_module/Function-geneset/PATHWAY_108_df_pValue.txt',
  #             quote = F,sep = '\t')
  rownames(MP_data)<-paste0(MP_name[i],'_',rownames(MP_data))
  all_MP_data<-rbind(all_MP_data,MP_data)
  
}
colnames(all_MP_data)<-paste0('MP_',1:13)

MP_order<-c('MP_11','MP_8','MP_10','MP_2','MP_14','MP_3','MP_12','MP_9',
            'MP_6',
            'MP_7','MP_4','MP_5','MP_13')
MP_order<-c('MP_8','MP_9','MP_2','MP_3','MP_6','MP_7','MP_14','MP_4',
            'MP_5',
            'MP_12','MP_13','MP_10','MP_11')
MP_order<-c('MP_7','MP_8','MP_1','MP_2','MP_5','MP_6','MP_13','MP_3','MP_4','MP_10','MP_12','MP_9','MP_13')
all_MP_data<-all_MP_data[,MP_order]
all_MP_data2<-all_MP_data[c(lapply(c('HALLMARK_ANDROGEN_RESPONSE','GOBP_CELLULAR_RESPONSE_TO_STEROID_HORMONE_STIMULUS',
                                     'GOBP_IRON_ION_TRANSPORT','GOBP_STEROID_HORMONE_MEDIATED_SIGNALING_PATHWAY',
                                     'pEMT'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP14 EMT-III','MP10 Protein maturation','GOBP_PROTEIN_TARGETING_TO_VACUOLE_INVOLVED_IN_AUTOPHAGY',
                                     'GOBP_PEPTIDE_METABOLIC_PROCESS','GOBP_POSITIVE_REGULATION_OF_TAU_PROTEIN_KINASE_ACTIVITY'),
                                   function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('Stress','GOBP_PROTEIN_REFOLDING','GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN',
                                     'GOBP_APOPTOTIC_SIGNALING_PATHWAY','GOBP_APOPTOTIC_PROCESS$'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('HALLMARK_COAGULATION','GOBP_NEGATIVE_REGULATION_OF_COAGULATION',
                                     'GOBP_REGULATION_OF_COAGULATION','GOBP_NEGATIVE_REGULATION_OF_WOUND_HEALING',
                                     'GOBP_ACUTE_PHASE_RESPONSE'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('GOBP_DEFENSE_RESPONSE$','GOBP_ACUTE_INFLAMMATORY_RESPONSE$',
                                     'GOBP_ACUTE_INFLAMMATORY_RESPONSE$','GOBP_PROTEIN_CONTAINING_COMPLEX_REMODELING',
                                     'Interferon'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP36 IG','GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY',
                                     'GOBP_ADAPTIVE_IMMUNE_RESPONSE_BASED_ON_SOMATIC_RECOMBINATION_OF_IMMUNE_RECEPTORS_BUILT_FROM_IMMUNOGLOBULIN_SUPERFAMILY_DOMAINS',
                                     'GOBP_LEUKOCYTE_MEDIATED_IMMUNITY','GOBP_IMMUNE_EFFECTOR_PROCESS'),
                                   function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP12 EMT-I','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                                     'MES1','^MP1$','Mesenchymal'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('HALLMARK_MYC_TARGETS_V1','Senescent cells','Senescence-Associated Secretory Phenotype',
                                     'Genome Instability and Mutation','GOBP_PROTEIN_FOLDING$'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP17 Interferon','GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_PEPTIDE_ANTIGEN',
                                     'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN',
                                     'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN$',
                                     'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II'),
                                   function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP18 Interferon','GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN',
                                     'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN',
                                     'GOBP_BIOLOGICAL_PROCESS_INVOLVED_IN_INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS',
                                     'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('^MP6$','MP19 Epithelial','Squamous','^MP8$','Epi dif. 1'),function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('MP25 Astrocytes','^AC$','GOBP_OXIDATIVE_PHOSPHORYLATION',
                                     'GOBP_AEROBIC_RESPIRATION','GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT'),
                                   function(x) grep(x,rownames(all_MP_data))) %>% unlist(),
                            lapply(c('GOBP_CELLULAR_RESPIRATION','GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN',
                                     'GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS',
                                     'GOBP_PROTON_TRANSMEMBRANE_TRANSPORT','GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY'),
                                   function(x) grep(x,rownames(all_MP_data))) %>% unlist()),]


select_pathway<-c()
for(j in 1:13){
  all_MP_data<-all_MP_data[order(all_MP_data[,j]),]
  bb<-all_MP_data[setdiff(rownames(all_MP_data),select_pathway),]
  #top_way<-rownames(all_MP_data)[1:5]
  
  select_pathway<-c(select_pathway,rownames(bb)[1:5])
}
select_pathway<-unique(select_pathway)


plot_data<-all_MP_data[select_pathway,]
plot_data<-all_MP_data2
plot_data<-(-log10(plot_data))

range(plot_data)
bk<-seq(0,10,length.out=100)
color_pheatmap<-c(#colorRampPalette(c("#0669AD",'#89BDD9'))(20),
  #colorRampPalette(c("#89BDD9",'white'))(30),
  colorRampPalette(c("white",'#E9C1C6'))(20),
  colorRampPalette(c("#E9C1C6",'#BF404D'))(80)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(plot_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = F,
                      cluster_cols = F,
                      #treeheight_row = T,treeheight_col = T,
                      #display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "black",
                      #fontsize = 10,
                      #cellwidth=15,
                      #cellheight=15,
                      main = "MP_pathway",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/re/MP_score/MP_GOBP_enrich.pdf',width = 18,height = 8)
print(p)
dev.off()




col_fun <- circlize::colorRamp2(
  c(0,3,10), 
  c( "white","#E9C1C6","#BF404D")
)
p2<-Heatmap(matrix = as.matrix(plot_data),
            name = '-logP',
            col = col_fun,row_names_gp = gpar(fontsize = 8),
            border='grey',
            cluster_rows =F,
            column_title = 'MP_pathway',
            cluster_columns = F
)
print(p2)

pdf('E:/Mirror/ST_analysis/pic/re/MP_score/MP_GOBP_enrich2.pdf',width = 10,height = 8)
print(p2)
dev.off()

?Heatmap




############################################################################################################################
############################################################################################################################
############################################################################################################################
############################################################################################################################
####最终使用GOBP的画
file_MP[12]
all_MP_data<-read.delim(paste0(dir_MP,file_MP[12]),stringsAsFactors = F,check.names = F)
colnames(all_MP_data)<-paste0('MP_',1:13)

MP_order<-c('MP_7','MP_8','MP_1','MP_2','MP_5','MP_6','MP_13','MP_3','MP_4','MP_11','MP_12','MP_9','MP_10')
all_MP_data<-all_MP_data[,MP_order]


select_pathway<-c()
for(j in 1:13){
  all_MP_data<-all_MP_data[order(all_MP_data[,j]),]
  bb<-all_MP_data[setdiff(rownames(all_MP_data),select_pathway),]
  #top_way<-rownames(all_MP_data)[1:5]
  
  select_pathway<-c(select_pathway,rownames(bb)[1:5])
}
select_pathway<-unique(select_pathway)


plot_data<-all_MP_data[select_pathway,]
plot_data<-(-log10(plot_data))

range(plot_data)
bk<-seq(0,10,length.out=100)
col_fun <- colorRamp2(
  c(0, 5, 15, 20, 100), 
  c("white","#FBF4C3",'#FBF4C3','#CC4476','#251F49')
)
color_pheatmap<-c(#colorRampPalette(c("#0669AD",'#89BDD9'))(20),
  #colorRampPalette(c("#89BDD9",'white'))(30),
  colorRampPalette(c("white",'#E9C1C6'))(20),
  colorRampPalette(c("#E9C1C6",'#BF404D'))(80)) ###"#CC281B"
#color_pheatmap<-colorRampPalette(c("#4575B4","white","#FF0033"))(100)  ###"#CC281B"
p<-pheatmap::pheatmap(as.matrix(plot_data), 
                      scale = "none",
                      color=color_pheatmap,
                      border_color = 'white',
                      # border='white',
                      cluster_rows = F,
                      cluster_cols = F,
                      #treeheight_row = T,treeheight_col = T,
                      #display_numbers = dis_dot,
                      na_col = "grey90",
                      fontsize_number=15,
                      number_color = "black",
                      #fontsize = 10,
                      #cellwidth=15,
                      #cellheight=15,
                      main = "MP_pathway",
                      breaks = bk,
                      name = 'scale_exp'
)
print(p)
pdf('E:/Mirror/ST_analysis/pic/re/MP_score/MP_GOBP_enrich_final.pdf',width = 18,height = 8)
print(p)
dev.off()


col_fun <- colorRamp2(
  c(0, 1, 2, 8, 10), 
  c("white","#FBF4C3",'#FBF4C3','#CC4476','#251F49')
)
p2<-Heatmap(matrix = as.matrix(plot_data),
            name = '-logP',
            col = col_fun,
            row_names_gp = gpar(fontsize = 8),
            border='grey',
            cluster_rows =F,
            column_title = 'MP_pathway',
            cluster_columns = F
)
print(p2)
pdf('E:/Mirror/ST_analysis/pic/re/MP_score/MP_GOBP_enrich_final.pdf',width = 8,height = 8)
print(p2)
dev.off()


















