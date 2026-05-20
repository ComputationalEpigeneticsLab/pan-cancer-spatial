
all_AUC <- read.table(file='/SpaCET_AUC.txt',sep = '\t',header = T)

####AUC大于0.7的算个比例
th_index<-0.7

####环形条状图
plot_data<-data.frame(AUC=c('>0.7','<=0.7'),
                      num=c(length(which(all_AUC$AUC>0.7)),length(which(all_AUC$AUC<=0.7))))###c(205,24)
plot_data$value<-plot_data$num/sum(plot_data$num)
plot_data$AUC<-factor(plot_data$AUC,levels = c('>0.7','<=0.7'))

p5 <- ggplot(plot_data,aes(x = 1, y = value, fill = AUC)) +
  geom_col(colour = "white")+ 
  coord_polar(theta = "y", start = 1.65) +
  geom_text(aes(label = paste0(round(value * 100, 2), "%")),
            position = position_fill(vjust = 0.5)) +
  scale_fill_manual(values=c(">0.7"="#f4b184","<=0.7"="#8fabdd"))+
  xlim(c(-0.2, 2)) +
  theme(
    panel.background = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
print(p5)
pdf(paste0("/Fs1/",'AUC_pie_0.7.pdf'),width = 7,height = 7)
print(p5)
dev.off()




