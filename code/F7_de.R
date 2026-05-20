hop_infor_imm <- read.table("\\new\\hop_infor_imm.txt",sep = "\t",header = T)



hop_infor_imm$cancer_type <- sapply(strsplit(hop_infor_imm$X, "_"), function(x) x[1])
hop_infor_imm$gse_id <- sapply(strsplit(hop_infor_imm$X, "_"), function(x) x[2])
hop_infor_imm$gsm_id <- sapply(strsplit(hop_infor_imm$X, "_"), function(x) x[3])

hop_infor_imm$slice <- paste0(hop_infor_imm$cancer_type,"_",hop_infor_imm$gse_id,"_",hop_infor_imm$gsm_id)

g_r_noTIB=c("GBM_GSE235672_GSM7507327","GBM_GSE235672_GSM7507330",
            "LIHC_GSE238264_GSM7661255","LIHC_GSE238264_GSM7661256","LIHC_GSE238264_GSM7661257","LIHC_GSE238264_GSM7661258",'LIHC_lihc03_slice7' )

g_nr_TIB=c("GBM_GSE235672_GSM7507312","GBM_GSE235672_GSM7507328","LIHC_GSE238264_GSM7661260",'LIHC_lihc03_slice5')

g_nr_noTIB=c("GBM_GSE235672_GSM7507323","GBM_GSE235672_GSM7507329","LIHC_GSE238264_GSM7661259","LIHC_GSE238264_GSM7661261")

g_nr <-c(g_nr_TIB,g_nr_noTIB)

hop_infor_imm_r <- hop_infor_imm[hop_infor_imm$slice %in% g_r_noTIB,]
hop_infor_imm_nr <- hop_infor_imm[hop_infor_imm$slice %in% g_nr,]

rctd_data <- read.table("\\new\\deconvolution_rctd_infor_imm.txt",sep = "\t",header = T)


library(tidyr)  # 提供 separate_rows() 函数
library(dplyr)  # 提供管道操作符 %>%

threehop_infor_imm_r <- hop_infor_imm_r[,c(6,7,8,21)]
step1_detail <- threehop_infor_imm_r %>%
  select(slice, step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) %>%
  mutate(step = "step_1",
         full_cell_name = paste0(slice, "_", step_1))

step2_detail <- threehop_infor_imm_r %>%
  select(slice, step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) %>%
  mutate(step = "step_2",
         full_cell_name = paste0(slice, "_", step_2))

step3_detail <- threehop_infor_imm_r %>%
  select(slice, step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) %>%
  mutate(step = "step_3",
         full_cell_name = paste0(slice, "_", step_3))

# 合并所有详细数据
colnames(step1_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step2_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step3_detail) <- c("slice" ,"cell","step","full_cell_name")
all_steps_detail_r <- rbind(step1_detail, step2_detail, step3_detail)
all_steps_detail_r <- all_steps_detail_r[!duplicated(all_steps_detail_r[, 4]), ]
threehop_imm_r <- merge(all_steps_detail_r,rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)


threehop_infor_imm_nr <- hop_infor_imm_nr[,c(6,7,8,21)]
step1_detail <- threehop_infor_imm_nr %>%
  select(slice, step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) %>%
  mutate(step = "step_1",
         full_cell_name = paste0(slice, "_", step_1))

step2_detail <- threehop_infor_imm_nr %>%
  select(slice, step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) %>%
  mutate(step = "step_2",
         full_cell_name = paste0(slice, "_", step_2))

step3_detail <- threehop_infor_imm_nr %>%
  select(slice, step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) %>%
  mutate(step = "step_3",
         full_cell_name = paste0(slice, "_", step_3))

# 合并所有详细数据
colnames(step1_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step2_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step3_detail) <- c("slice" ,"cell","step","full_cell_name")
all_steps_detail_nr <- rbind(step1_detail, step2_detail, step3_detail)
all_steps_detail_nr <- all_steps_detail_nr[!duplicated(all_steps_detail_nr[, 4]), ]
threehop_imm_nr <- merge(all_steps_detail_nr,rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)



##三组
hop_infor_imm_nr_TIB <- hop_infor_imm[hop_infor_imm$slice %in% g_nr_TIB,]
hop_infor_imm_nr_noTIB <- hop_infor_imm[hop_infor_imm$slice %in% g_nr_noTIB,]


threehop_infor_imm_nr_TIB <- hop_infor_imm_nr_TIB[,c(6,7,8,21)]
step1_detail <- threehop_infor_imm_nr_TIB %>%
  select(slice, step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) %>%
  mutate(step = "step_1",
         full_cell_name = paste0(slice, "_", step_1))

step2_detail <- threehop_infor_imm_nr_TIB %>%
  select(slice, step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) %>%
  mutate(step = "step_2",
         full_cell_name = paste0(slice, "_", step_2))

step3_detail <- threehop_infor_imm_nr_TIB %>%
  select(slice, step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) %>%
  mutate(step = "step_3",
         full_cell_name = paste0(slice, "_", step_3))

# 合并所有详细数据
colnames(step1_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step2_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step3_detail) <- c("slice" ,"cell","step","full_cell_name")
all_steps_detail_nr_TIB <- rbind(step1_detail, step2_detail, step3_detail)
all_steps_detail_nr_TIB <- all_steps_detail_nr_TIB[!duplicated(all_steps_detail_nr_TIB[, 4]), ]
threehop_imm_nr_TIB <- merge(all_steps_detail_nr_TIB,rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)



threehop_infor_imm_nr_noTIB <- hop_infor_imm_nr_noTIB[,c(6,7,8,21)]
step1_detail <- threehop_infor_imm_nr_noTIB %>%
  select(slice, step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) %>%
  mutate(step = "step_1",
         full_cell_name = paste0(slice, "_", step_1))

step2_detail <- threehop_infor_imm_nr_noTIB %>%
  select(slice, step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) %>%
  mutate(step = "step_2",
         full_cell_name = paste0(slice, "_", step_2))

step3_detail <- threehop_infor_imm_nr_noTIB %>%
  select(slice, step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) %>%
  mutate(step = "step_3",
         full_cell_name = paste0(slice, "_", step_3))

# 合并所有详细数据
colnames(step1_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step2_detail) <- c("slice" ,"cell","step","full_cell_name")
colnames(step3_detail) <- c("slice" ,"cell","step","full_cell_name")
all_steps_detail_nr_noTIB <- rbind(step1_detail, step2_detail, step3_detail)
all_steps_detail_nr_noTIB <- all_steps_detail_nr_noTIB[!duplicated(all_steps_detail_nr_noTIB[, 4]), ]
threehop_imm_nr_noTIB <- merge(all_steps_detail_nr_noTIB,rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)

########两组nr_tib vs r
###B
b_r <- data.frame(
  group = "threehop_imm_r",
  count = threehop_imm_r$B.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
b_nr <- data.frame(
  group = "threehop_imm_nr_TIB",
  count = threehop_imm_nr_TIB$B.lymphocytes
)

# 合并两个数据框
b_threestep <- rbind(b_r, b_nr)

p1 <- ggplot(b_threestep, aes(x = group, y = count, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr_TIB" = "#15629e")) +
  labs(title = "Boxplot of B.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.15) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.15) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p1)
pdf("\\new\\Boxplot_Blymphocytes_public_nr_TIBvsr1.pdf",width = 5,height = 7)
print(p1)
dev.off()

###T
t_r <- data.frame(
  group = "threehop_imm_r",
  count = threehop_imm_r$T.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
t_nr <- data.frame(
  group = "threehop_imm_nr_TIB",
  count = threehop_imm_nr_TIB$T.lymphocytes
)

# 合并两个数据框
t_threestep <- rbind(t_r, t_nr)

p2 <- ggplot(t_threestep, aes(x = group, y = count, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr_TIB" = "#15629e")) +
  labs(title = "Boxplot of T.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.1) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p2)

pdf("\\new\\Boxplot_Tlymphocytes_public_nr_TIBvsr1.pdf",width = 5,height = 7)
print(p2)
dev.off()

########两组nr vs r
###B
b_r <- data.frame(
  group = "threehop_imm_r",
  count = threehop_imm_r$B.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
b_nr <- data.frame(
  group = "threehop_imm_nr",
  count = threehop_imm_nr$B.lymphocytes
)

# 合并两个数据框
b_threestep <- rbind(b_r, b_nr)

p1 <- ggplot(b_threestep, aes(x = group, y = count, fill = group)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr" = "#15629e")) +
  labs(title = "Boxplot of B.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.15) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.15) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


print(p1)
pdf("\\new\\Boxplot_Blymphocytes_public.pdf",width = 5,height = 7)
print(p1)
dev.off()

###T
t_r <- data.frame(
  group = "threehop_imm_r",
  count = threehop_imm_r$T.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
t_nr <- data.frame(
  group = "threehop_imm_nr",
  count = threehop_imm_nr$T.lymphocytes
)

# 合并两个数据框
t_threestep <- rbind(t_r, t_nr)

p2 <- ggplot(t_threestep, aes(x = group, y = count, fill = group)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr" = "#15629e")) +
  labs(title = "Boxplot of T.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.1) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.1) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p2)

pdf("\\new\\Boxplot_Tlymphocytes_public.pdf",width = 5,height = 7)
print(p2)
dev.off()




###ESCC
XYZ_rctd_data <- read.table("\\new\\ESCC\\XYZ_slice1_Deconvolution_1_2_9_10_12_public.txt",sep = "\t",header = T)
BGM_rctd_data <- read.table("\\new\\ESCC\\BGM_slice1_Deconvolution_1_2_9_10_12_public.txt",sep = "\t",header = T)

XYZ_rctd_data$X <- rownames(XYZ_rctd_data)
BGM_rctd_data$X <- rownames(BGM_rctd_data)
ESCC_rctd_data <- rbind(XYZ_rctd_data,BGM_rctd_data)

ESCC_R_hop_infor_imm <- read.table("\\new\\ESCC\\XYZ_slice1_nearSpotStep1to10.txt",sep = "\t",header = T)
ESCC_NR_hop_infor_imm <- read.table("\\new\\ESCC\\BGM_slice1_nearSpotStep1to10.txt",sep = "\t",header = T)

ESCC_R_hop_infor_imm <- ESCC_R_hop_infor_imm[,c(1,7,8,9,10,11)]
ESCC_NR_hop_infor_imm <- ESCC_NR_hop_infor_imm[,c(1,7,8,9,10,11)]


step1_detail <- ESCC_R_hop_infor_imm %>%
  select(step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) 

step2_detail <- ESCC_R_hop_infor_imm %>%
  select(step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) 

step3_detail <- ESCC_R_hop_infor_imm %>%
  select(step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) 
colnames(step1_detail) <- c("full_cell_name")
colnames(step2_detail) <- c("full_cell_name")
colnames(step3_detail) <- c("full_cell_name")

ESCC_threehop_imm_r <- rbind(step1_detail,step2_detail,step3_detail)
ESCC_threehop_imm_r <- unique(ESCC_threehop_imm_r)
ESCC_threehop_imm_r <- merge(ESCC_threehop_imm_r,XYZ_rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)



step1_detail <- ESCC_NR_hop_infor_imm %>%
  select(step_1) %>%
  separate_rows(step_1, sep = ",") %>%
  mutate(step_1 = trimws(step_1)) 

step2_detail <- ESCC_NR_hop_infor_imm %>%
  select(step_2) %>%
  separate_rows(step_2, sep = ",") %>%
  mutate(step_2 = trimws(step_2)) 

step3_detail <- ESCC_NR_hop_infor_imm %>%
  select(step_3) %>%
  separate_rows(step_3, sep = ",") %>%
  mutate(step_3 = trimws(step_3)) 
colnames(step1_detail) <- c("full_cell_name")
colnames(step2_detail) <- c("full_cell_name")
colnames(step3_detail) <- c("full_cell_name")

ESCC_threehop_imm_nr <- rbind(step1_detail,step2_detail,step3_detail)
ESCC_threehop_imm_nr <- unique(ESCC_threehop_imm_nr)
ESCC_threehop_imm_nr <- merge(ESCC_threehop_imm_nr,BGM_rctd_data,by.x = "full_cell_name",by.y = "X",all.x = T)

########两组
###B
b_r <- data.frame(
  group = "threehop_imm_r",
  count = ESCC_threehop_imm_r$B.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
b_nr <- data.frame(
  group = "threehop_imm_nr",
  count = ESCC_threehop_imm_nr$B.lymphocytes
)

# 合并两个数据框
b_threestep <- rbind(b_r, b_nr)

p5 <- ggplot(b_threestep, aes(x = group, y = count, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr" = "#15629e")) +
  labs(title = "Boxplot of B.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.0005) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.0005) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p5)
pdf("\\new\\Boxplot_Blymphocytes_ESCC.pdf",width = 5,height = 7)
print(p5)
dev.off()

###T
t_r <- data.frame(
  group = "threehop_imm_r",
  count = ESCC_threehop_imm_r$T.lymphocytes
)

# 创建三个hop_imm_nr的B细胞数据
t_nr <- data.frame(
  group = "threehop_imm_nr",
  count = ESCC_threehop_imm_nr$T.lymphocytes
)

# 合并两个数据框
t_threestep <- rbind(t_r, t_nr)

p6 <- ggplot(t_threestep, aes(x = group, y = count, fill = group)) +
  geom_boxplot(alpha = 0.7, width = 0.6, 
               outlier.shape = NA,  # 去掉离散点
               na.rm = TRUE) +      # 忽略缺失值
  scale_fill_manual(values = c("threehop_imm_r" = "#ddb424", 
                               "threehop_imm_nr" = "#15629e")) +
  labs(title = "Boxplot of T.lymphocytes",
       x = "Group",
       y = "Count",
       fill = "Group") +
  # 添加p值
  stat_compare_means(method = "wilcox.test",  # 使用Wilcoxon检验（非参数）
                     label = "p.format",      # 显示p值格式
                     label.x = 1.5,           # p值标签的x位置
                     label.y = 0.15) +        # p值标签的y位置
  # 设置y轴范围
  ylim(0, 0.15) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p6)

pdf("\\new\\Boxplot_Tlymphocytes_ESCC.pdf",width = 5,height = 7)
print(p6)
dev.off()

