###TAM

all_sliceTAM<-read.delim('/Fs6/all_step_TAM.txt',stringsAsFactors = F,check.names = F)



all_sliceTAM$cancer<-unlist(lapply(strsplit(all_sliceTAM$slice,'/'),function(x)x[1]))
all_sliceTAM <- all_sliceTAM[!is.na(all_sliceTAM$cellType) & all_sliceTAM$cellType != "", ]


all_sliceTAM$cellType[all_sliceTAM$cellType=="FCN1_TAM"] <- 'SPP1_TAM'
all_sliceTAM$cellType<-factor(all_sliceTAM$cellType,levels = c("SPP1_TAM",'C1Q_TAM',"Microglial_TAM","blood_TAM"))
p_compare<-ggplot(all_sliceTAM,aes(x=cancer,y=Freq,fill=cellType)) +
  geom_bar(stat = "identity",position="fill") + ###,color="white"  边框
  #coord_flip()+
  scale_fill_manual(values = c("C1Q_TAM"="#71a3c5",'FCN1_TAM'='#cc889f',"SPP1_TAM"="#b83231",
                               "Microglial_TAM"="#80c598","blood_TAM"="#1d804e"))+
  #geom_text(size = 4, position = position_stack(vjust = 0.5),colour = "black")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1),
        panel.background = element_blank(),
        axis.line = element_line(),
        axis.text = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 15))+
  xlab("cancer")+ylab("")+
  #guides(fill = "none")+
  #ylim(0, 1)+
  ggtitle('subTAM')
print(p_compare)



