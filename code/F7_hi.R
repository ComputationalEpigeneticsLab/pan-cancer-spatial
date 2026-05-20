##public代谢可视化箱式图
immune_metabolic_adaptation <- read.table(file = "/metabolism/hop1-3_scFEA_immune_metabolic_adaptation.txt",header = T,sep = ",")
scFEA_tumor_anabolism <- read.table(file = "/metabolism/hop1-3_scFEA_tumor_anabolism.txt",header = T,sep = ",")
immune_metabolic_adaptation <- immune_metabolic_adaptation[,1:2]
scFEA_tumor_anabolism <- scFEA_tumor_anabolism[,1:2]

##immune_metabolic_adaptation
r_values <- immune_metabolic_adaptation$r_immune_metabolic_adaptation # 填入所有r值
nr_values <- immune_metabolic_adaptation$nr_immune_metabolic_adaptation  # 填入所有nr值

r_values <- as.numeric(unlist(strsplit(r_values, ",")))
nr_values <- as.numeric(unlist(strsplit(nr_values, ",")))

# 创建数据框
data <- data.frame(
  value = c(r_values, nr_values),
  group = c(rep("r", length(r_values)), 
            rep("nr", length(nr_values)))
)


p_box1 <- ggplot(data, aes(x = group, y = value, fill = group)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = group), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 0.7, outlier.shape = NA) +
  
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  scale_color_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  ggtitle("immune_metabolic_adaptation") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",
                     label.x = 1.5,
                     label.y = 0.11) +
  coord_cartesian(ylim = c(0, 0.12)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p_box1)


pdf('/metabolism/immune_metabolic_adaptation_public.pdf',width = 4,height = 4)
print(p_box1)
dev.off()


##scFEA_tumor_anabolism
r_values <- scFEA_tumor_anabolism$r_tumor_anabolism # 填入所有r值
nr_values <- scFEA_tumor_anabolism$nr_tumor_anabolism  # 填入所有nr值

r_values <- as.numeric(unlist(strsplit(r_values, ",")))
nr_values <- as.numeric(unlist(strsplit(nr_values, ",")))

# 创建数据框
data <- data.frame(
  value = c(r_values, nr_values),
  group = c(rep("r", length(r_values)), 
            rep("nr", length(nr_values)))
)

p_box2 <- ggplot(data, aes(x = group, y = value, fill = group)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = group), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 0.7, outlier.shape = NA) +
  
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  scale_color_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  ggtitle("scFEA_tumor_anabolism") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",
                     label.x = 1.5,
                     label.y = 0.11) +
  coord_cartesian(ylim = c(0, 0.12)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p_box2)


pdf('/metabolism/scFEA_tumor_anabolism_public.pdf',width = 4,height = 4)
print(p_box2)
dev.off()









##ESCC代谢可视化箱式图
immune_metabolic_adaptation <- read.table(file = "/metabolism/hop1-3_scFEA_immune_metabolic_adaptation_escc.txt",header = T,sep = ",")
scFEA_tumor_anabolism <- read.table(file = "/metabolism/hop1-3_scFEA_tumor_anabolism_escc.txt",header = T,sep = ",")

##immune_metabolic_adaptation
r_values <- immune_metabolic_adaptation$r_immune_metabolic_adaptation # 填入所有r值
nr_values <- immune_metabolic_adaptation$nr_immune_metabolic_adaptation  # 填入所有nr值

r_values <- as.numeric(unlist(strsplit(r_values, ",")))
nr_values <- as.numeric(unlist(strsplit(nr_values, ",")))

# 创建数据框
data <- data.frame(
  value = c(r_values, nr_values),
  group = c(rep("r", length(r_values)), 
            rep("nr", length(nr_values)))
)


p_box3 <- ggplot(data, aes(x = group, y = value, fill = group)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = group), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 0.7, outlier.shape = NA) +
  
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  scale_color_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  ggtitle("immune_metabolic_adaptation_escc") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",
                     label.x = 1.5,
                     label.y = 0.06) +
  coord_cartesian(ylim = c(0.015, 0.07)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p_box3)


pdf('/metabolism/immune_metabolic_adaptation_escc.pdf',width = 4,height = 4)
print(p_box3)
dev.off()


##scFEA_tumor_anabolism
r_values <- scFEA_tumor_anabolism$r_tumor_anabolism # 填入所有r值
nr_values <- scFEA_tumor_anabolism$nr_tumor_anabolism  # 填入所有nr值

r_values <- as.numeric(unlist(strsplit(r_values, ",")))
nr_values <- as.numeric(unlist(strsplit(nr_values, ",")))

# 创建数据框
data <- data.frame(
  value = c(r_values, nr_values),
  group = c(rep("r", length(r_values)), 
            rep("nr", length(nr_values)))
)

p_box4 <- ggplot(data, aes(x = group, y = value, fill = group)) + 
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(aes(fill = group), color = 'black', width = 0.8,
               position = position_dodge(0.9), alpha = 0.7, outlier.shape = NA) +
  
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = 'black')) +
  scale_fill_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  scale_color_manual(values = c("nr" = "#15629e", "r" = "#ddb424")) +
  ggtitle("scFEA_tumor_anabolism_escc") +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox.test",
                     label.x = 1.5,
                     label.y = 0.06) +
  coord_cartesian(ylim = c(0.015, 0.07)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10))

print(p_box4)


pdf('/metabolism/scFEA_tumor_anabolism_escc.pdf',width = 4,height = 4)
print(p_box4)
dev.off()


