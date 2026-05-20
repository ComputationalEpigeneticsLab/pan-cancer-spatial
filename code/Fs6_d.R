###################################################
####按每个step画
pdf('D:/pan_cancer/0_修稿/Fs6/subTAM_step_SPP1&FN1.pdf', width = 6, height = 10)

# 确保步骤顺序正确
steps <- unique(all_step_TAM$step)
steps <- steps[order(steps)]  # 按step1, step2...排序

for(i in 1:length(steps)){
  # 获取当前step的数据
  current_step <- steps[i]
  stat_data <- all_step_TAM[all_step_TAM$step == current_step, ]
  
  # 转换数据格式
  # 使用acast时指定fun.aggregate = sum
  p_data <- as.data.frame(
    reshape2::acast(
      stat_data[, c('cellType', 'slice', 'ratio')], 
      slice ~ cellType,
      value.var = "ratio",
      fun.aggregate = sum,  # 改为sum
      na.rm = TRUE
    )
  )
  
  # 确保列顺序正确
  required_cols <- c('SPP1_TAM', 'C1Q_TAM', 'Microglial_TAM', 'blood_TAM')
  missing_cols <- setdiff(required_cols, colnames(p_data))
  
  # 如果缺少某些列，添加全为0的列
  for(col in missing_cols) {
    p_data[[col]] <- 0
  }
  
  # 按所需顺序排列列
  p_data <- p_data[, required_cols]
  
  # 提取癌症类型 - 更健壮的提取方法
  p_data$cancer <- sapply(strsplit(rownames(p_data), "/"), function(x) {
    # 提取第一部分（如BRCA）
    cancer_type <- x[1]
    
    # 清理可能的数字后缀（如BRCA01 -> BRCA）
    cancer_type <- gsub("[0-9]+$", "", cancer_type)  # 移除末尾数字
    cancer_type <- toupper(cancer_type)  # 转为大写
    
    # 特殊处理一些情况
    if(cancer_type == "BRCA/BRCA") cancer_type <- "BRCA"
    if(cancer_type == "BRCA") cancer_type <- "BRCA"
    
    return(cancer_type)
  })
  
  # 标准化癌症类型（确保在颜色映射中存在）
  p_data$cancer <- factor(p_data$cancer, levels = cancer_types)
  all_colors <- cancer_colors
  
  # 创建行注释
  la <- rowAnnotation(
    df = data.frame(cancer = p_data$cancer),
    col = list(cancer = all_colors)
    
  )
  
  # 准备热图数据
  heat_data <- as.matrix(p_data[, 1:4])
  
  # 设置颜色标度
  col_fun <- colorRamp2(
    seq(0, 0.4, length.out = 10), 
    rev(c('#ab2c74','#B84687','#C5619B','#D98AB8','#e7a5cc','#ECCBD6','#EEDADA','#f0eadf','#bbdf8e','#82AF5D'))
  )
  
  # 创建热图
  p1 <- Heatmap(
    heat_data,
    col = col_fun,
    left_annotation = la,
    show_row_names = FALSE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = "ratio",
    column_title = paste(current_step, "- TAM Subtype Ratios"),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm"),
    row_split = p_data$cancer,
    row_title = "Samples",
    heatmap_legend_param = list(
      title = "Ratio",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )
  
  print(p1)
}

dev.off()


