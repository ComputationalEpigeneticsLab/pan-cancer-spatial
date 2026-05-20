slice_re <- read.table(file = '/Fs3/Texhausted_Bavidity_CoreRatio.txt',
                       stringsAsFactors = F,check.names = F)

dir_picc<-'/Fs3/'

plot_data<-data.frame(core_ratio=rep(slice_re$core_ratio,2),
                      T_exhausted=c(slice_re$T_exhausted_Zhangcl,slice_re$T_exhausted_Xug),
                      exhausted_type=rep(c('T_exhausted_Zhangcl','T_exhausted_Xug'),each=nrow(slice_re)))
plot_data<-plot_data[which(plot_data$T_exhausted!=0),]
plot_data<-plot_data[which(plot_data$exhausted_type=='T_exhausted_Xug'),]

pdf(paste0(dir_picc,'T_exhausted_core2.pdf'),width = 5.5,height = 3.5)
b <- ggplot(plot_data, aes(x = core_ratio, y = T_exhausted))
pp<-b + geom_point(aes(color = exhausted_type),size=0.01)+
  geom_smooth(aes(color = exhausted_type, fill = exhausted_type), method = "lm",se=T) +
  #geom_rug(aes(color =exhausted_type)) +
  scale_color_manual(values = c("#53b654", "#16afe2"))+
  scale_fill_manual(values = c("#53b654", "#16afe2"))+
  ggpubr::stat_cor(aes(color = exhausted_type), label.x = 0.1)+
  theme_classic()+
  ggtitle('T_exhausted')
print(pp)
dev.off()
?stat_cor

library(ggpubr)
all_TLS_data<-read.delim('/Fs3/all_slice_TLS.txt',stringsAsFactors = F,check.names = F)

table(all_TLS_data$cancer)
table(all_TLS_data$FinalLocalType)

plot_TLS<-all_TLS_data[all_TLS_data$FinalLocalType%in%'Core',]
TLS_type<-c('TLS_2021','TLS_31942071','TLS_CancerSRT')

dir_picc<-'/Fs3/'
pdf(paste0(dir_picc,'panCancer_TLS_cor2.pdf'),width = 5.5,height = 3.5)
#png(paste0(dir_picc,'panCancer_TLS_cor2.png'),width = 300,height = 300)
i='TLS_31942071'
p_data<-plot_TLS[,c(i,'CytoTRACE')]
colnames(p_data)<-c('TLS','stem')

b <- ggplot(p_data, aes(x = stem, y = TLS))
pp<-b + geom_point(color = "#f6c046",size=0.01,alpha=0.5)+
  geom_smooth(method = "lm", color = "black", fill = "gray60",size=0.7,se=T)+
  #geom_rug(aes(color =TLS_type)) +
  # scale_color_manual(values = c("#f6c046"))+
  # scale_fill_manual(values = c("#f6c046"))+
  stat_cor(method = "pearson", label.x = 0, label.y = 0.5,size=6)+
  theme_classic()+
  theme(axis.text = element_text(color = 'black'))+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=15),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 15) )
  #+ggtitle(paste0('pan_cancer_',i))
print(pp)

dev.off()