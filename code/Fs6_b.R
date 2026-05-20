dir_sc<-'/data/SC_data/re/'
file_sc<-list.files(pattern = '_meta.txt',path = dir_sc,recursive = T)
file_sc

for(i in 1:length(file_sc)){
  print(file_sc[i])
  sc_meta<-read.delim(paste0(dir_sc,file_sc[i]),stringsAsFactors = F,check.names = F)
  # print(unique(sc_meta$subCAFTAM))
  print(paste0('B: ','B lymphocytes'%in%unique(sc_meta$subCAFTAM)))
  print(paste0('T: ','T lymphocytes'%in%unique(sc_meta$subCAFTAM)))
}


dir_bdy<-'/data/10X_Visium/new_st_hop10/'
file_near<-list.files(pattern = '_nearSpotStep1to10.txt',path = dir_bdy,recursive = T)
dataset_slice<-unlist(lapply(strsplit(file_near,'_nearSpot'),function(x) x[1]))
length(file_near)
file_near[1:10]
dataset_slice[1:10]

dir_RCTD<-'/data/10X_Visium/RCTD3//'
file_RCTD<-list.files(pattern = 'Deconvolution.txt',path = dir_RCTD,recursive = T)
file_RCTD[1:10]
table(unlist(lapply(strsplit(file_RCTD,'_Deco'),function(x) x[1]))==dataset_slice)
cancer<-unlist(lapply(strsplit(file_near,'/'),function(x)x[1]))
# cancer<-substr(cancer,1,nchar(cancer)-2) %>% unique()
cancer[1:10]
cancer<-unique(cancer)

all_FC<-c()
all_P<-c()
for(i in 1:length(cancer)){#i=11
  file_near_cancer<-file_near[grep(cancer[i],file_near)]
  file_RCTD_cancer<-file_RCTD[grep(cancer[i],file_RCTD)]
  dataset_slice_cancer<-unlist(lapply(strsplit(file_near_cancer,'_nearSpot'),function(x) x[1]))
  
  CAF_TAM<-c()
  Imm_other<-c()
  for(j in 1:length(file_near_cancer)){#j=1
    st_near<-read.delim(paste0(dir_bdy,file_near_cancer[j]),stringsAsFactors = F,check.names = F)
    st_RCTD<-read.delim(paste0(dir_RCTD,file_RCTD_cancer[j]),stringsAsFactors = F,check.names = F)
    
    spot_near<-unlist(strsplit(st_near[st_near$FinalLocalType%in%c('Boundary','Budding'),'step_1'],',')) %>% unique()
    spot_near<-st_near[spot_near,]
    spot_near<-spot_near$cell_name[spot_near$FinalLocalType%in%c('Immune','Normal')]
    
    bdy_dis_RCTD<-st_RCTD[spot_near,]
    CAF_TAM<-c(CAF_TAM,bdy_dis_RCTD$CAF+bdy_dis_RCTD$TAM)
    if(length(intersect(c('B lymphocytes','T lymphocytes'),colnames(bdy_dis_RCTD)))==2){
      Imm_other<-c(Imm_other,apply(bdy_dis_RCTD[,intersect(c('B lymphocytes','T lymphocytes'),colnames(bdy_dis_RCTD))],1,sum))
    }else{
      Imm_other<-c(Imm_other,bdy_dis_RCTD[,intersect(c('B lymphocytes','T lymphocytes'),colnames(bdy_dis_RCTD))])
    }
    
  }
  
  FC<-mean(CAF_TAM)/mean(Imm_other)
  P_value<-wilcox.test(CAF_TAM,Imm_other,alternative ='greater')$p.value
  
  all_FC<-c(all_FC,FC)
  all_P<-c(all_P,P_value)
}


plot_d<-data.frame(logFC=log2(all_FC),
                   logP=-log10(all_P),
                   P_value=all_P,
                   cancer=cancer,
                   x='logFC')
write.table(plot_d,'/data/10X_Visium/RCTD3/step1_compare.txt',
            quote = F,sep = '\t',row.names = F)
plot_d<-read.delim('/Fs6/FigS6ab/step1_compare.txt',
                   stringsAsFactors = F,check.names = F)

plot_d$logP[which(plot_d$logP>10)]<-10
plot_d$color<-'white'
plot_d$color[which(plot_d$P_value<0.05)]<-'#336633'####p值显著的画圈 也可以不画

p_dot<-ggplot(plot_d, aes(x=logFC, y=cancer,size=logFC,fill=logP)) +
  geom_point(shape = 21,  color = plot_d$color) + # 使用shape = 21画圈
  geom_vline(xintercept=c(-log2(1.5),log2(1.5)),lty=2,col="black",lwd=0.6) +
  scale_fill_gradientn(colours = c(colorRampPalette(c("#899DA4","#FBEDD2"))(10),
                                   colorRampPalette(c("#FBEDD2","#C13710"))(90)) )+ #设置填充颜色
  scale_size_continuous(range = c(3, 10))+
  #geom_text(aes(label = num), vjust = 0.5,size = 5)+
  theme_minimal()
print(p_dot)
pdf('/Fs6/cancer_step1_celltypeCompare_T_B.pdf',height = 11,width = 6)
print(p_dot)
dev.off()





