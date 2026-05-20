####各区域相关marker的柱状散点图绘制
###bdy
select_gene<-'CXCL17'
plot_data<-data.frame(gene=st_rds@assays[["Spatial"]]@counts[select_gene,],
                      localType=st_bdy$FinalLocalType)
spot_type<-c('Core',"Boundary","Dispersion",'Immune', 'Normal')
plot_data<-plot_data[plot_data$localType%in%spot_type,]
plot_data<-lapply(c('Core',"Boundary","Dispersion",'Immune', 'Normal'),function(x){#x='Core'
  xx<-plot_data$gene[which(plot_data$localType==x)]
  xx<-xx[which(xx<min(boxplot.stats(xx)$out))]
  xx_data<-data.frame(gene=xx,localType=x)
})
plot_data<-do.call(rbind,plot_data)
plot_data<-mutate(plot_data,Legend = factor(plot_data$localType, levels = spot_type))

median(plot_data$gene[which(plot_data$Legend=='Dispersion')])
median(plot_data$gene[which(plot_data$Legend=='Boundary')])
e <- ggplot(plot_data, aes(x = Legend, y = gene,fill=Legend))+ 
  geom_jitter(aes(color = Legend),size = 1,alpha=0.6,show.legend = F,
              position=position_jitterdodge(jitter.width = 0.8, 
                                            jitter.height = 0, 
                                            dodge.width = 1)) + # 不重叠的散点图
  stat_summary(fun.data = "median_q1q3", geom = "errorbar", width = 0.3, size = 0.5,position = position_dodge(0.8),show.legend = F) + # 误差棒，中位数，25%和75%分位数
  stat_summary(aes(fill = Legend), fun.y = median, geom = "crossbar", width = 0.6, size = 0.3,position = position_dodge(0.8),show.legend = F) + # 中位数水平线
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'))+
  scale_fill_manual(values = c("Boundary"="#FF9900","Core"="#990033","Dispersion"="#CD5A5A","Immune"="#5477AF","Normal"="#669933"))+
  scale_color_manual(values = c("Boundary"="#FF9900","Core"="#990033","Dispersion"="#CD5A5A","Immune"="#5477AF","Normal"="#669933"))+
  ggtitle(paste0('Bdy_',select_gene))+
  theme(plot.title = element_text(hjust = 0.4))+
  theme(plot.title = element_text(size = 12))+
  stat_compare_means(label.x = 1)+
  #theme_bw()+
  theme(axis.title.x = element_text(size=12),axis.text.x = element_text(size=10),
        axis.title.y =element_text(size=12),axis.text.y = element_text(size = 10) )
print(e)
?geom_jitter