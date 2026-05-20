wc_test = wilcox.test(MP_eco[MP_eco$eco_new%in%'eco3-4',]$metabolic, MP_eco[MP_eco$eco_new%in%'other',]$metabolic)
MP_eco$eco_new = factor(MP_eco$eco_new,levels = c('other','eco3-4'))
pdata = MP_eco[MP_eco$metabolic<0.018,]
p=ggplot(pdata, aes(x = metabolic, y=eco_new, fill = eco_new)) +
  #geom_density(aes(color = eco_new))+
  geom_density_ridges(alpha = 0.7,color='white',scale = 1.5)+
  #theme_classic(base_line_size = 1) +
  theme_bw()+
  labs(x="scFEA", y="density",title = 'Core ecosystem') +
  annotate("text", x = 0.014, y = 3,
           label = paste0("p ", format.pval(wc_test$p.value, digits = 3)),
           color = "black", size = 3) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))+
  scale_fill_manual(values = c('other'='#7AADD2','eco3-4'='#D390A7'))+
  theme(plot.title = element_text(size = 10,colour = "black",face="bold",hjust = 0.5),
        axis.title.x = element_text(size = 9,color = "black",face = "bold", vjust = 0.5, hjust = 0.5,angle = 0),
        legend.title = element_text(color="black",size=9,face="bold"),
        legend.text = element_text(color="black",size = 8),
        legend.position = 'none',
        axis.text.x = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0),
        axis.text.y = element_text(size = 8,color = "black",vjust = 0.5,hjust = 0.5,angle = 0))
p
ggsave("scFEA_Core中生态型代谢活性差异(部分代谢通路-分组)-1.pdf",p,device = "pdf",
       path = "E:/panCancerST/metabolism/scFEA/山脊图",height = 2,width = 3)