###CytoSpace结果计算 评估细胞类型比例  结果统计
library(ggplot2)

dir_out<-'E:/Mirror/ST_analysis/data/10X Visium/new_st_Cytospace/'
file_cyto<-list.files(pattern = '_cellNum.txt',path = dir_out,recursive = T)
file_cyto[1:10]
length(file_cyto)
dataSlice<-unlist(lapply(strsplit(file_cyto,'_cellNum'),function(x)x[1]))

dir_bdy<-'E:/Mirror/ST_analysis/data/10X Visium/new_st_copykat/'
file_bdy<-paste0(dataSlice,'_BdyCoreBud.txt')

all_cell_num<-c()
for(i in 1:length(dataSlice)){
  # i=1
  st_cyto<-read.delim(paste0(dir_out,file_cyto[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-read.delim(paste0(dir_bdy,file_bdy[i]),stringsAsFactors = F,check.names = F)
  st_bdy<-st_bdy[st_bdy$FinalLocalType%in%c('Immune','Normal'),]
  
  st_cyto<-st_cyto[intersect(rownames(st_cyto),rownames(st_bdy)),grep('CAF',colnames(st_cyto))]
  cell_num<-apply(st_cyto,2,sum)
  cell_num<-data.frame(cellType=names(cell_num),
                       cell_num=cell_num,
                       slice=dataSlice[i])
  
  all_cell_num<-rbind(all_cell_num,cell_num)
}

dir_pic<-'D:/pan_cancer/0_修稿/F4/Fig4f/'
all_cell_num <- read.table(file = paste0(dir_pic,"all_step_CAF.txt"), header = T, sep = "\t")
all_cell_num$cancer<-unlist(lapply(strsplit(all_cell_num$slice,'/'),function(x)x[1]))
unique(all_cell_num$cellType)

all_cell_num <- all_cell_num[all_cell_num$cellType != "",]
all_cell_num$cellType<-factor(all_cell_num$cellType,levels = c('apCAF','iCAF',"ifnCAF","rCAF","mCAF","vCAF","dCAF","tCAF"))

p_compare<-ggplot(all_cell_num,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("apCAF"="#c26701",'iCAF'='#e2a503',"ifnCAF"="#ffc000","rCAF"="#bbc100",
                               "mCAF"="#007182","vCAF"="#00827c","dCAF"="#4a806c","tCAF"="#737f6b"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subCAF_allStroma')
print(p_compare)

p_compare <- ggplot(all_cell_num, aes(x = cancer, y = Freq, fill = cellType)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ step, ncol = 1) +  # 按步骤分面
  scale_fill_manual(values = c("apCAF" = "#c26701",
                               'iCAF' = '#e2a503',
                               "ifnCAF" = "#ffc000",
                               "rCAF" = "#bbc100",
                               "mCAF" = "#007182",
                               "vCAF" = "#00827c",
                               "dCAF" = "#4a806c",
                               "tCAF" = "#737f6b")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 15),
        strip.background = element_rect(fill = "lightgrey")) +
  xlab("Cancer Type") + 
  ylab("Proportion") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle('CAF Subtypes Across Cancer Types and Steps')



pdf(paste0(dir_pic,'subCAF_bar_step1_5.pdf'),height =15,width = 6.5)
print(p_compare)
dev.off()

