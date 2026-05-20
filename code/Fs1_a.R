cancer_types <- c("GBM", "HNSCC", "HN-AS", "DIPG", "PN", "PCNSL", "DLBCL", "LGACC", "OSCC", 
                  "NPC", "TC", "DSRCT", "PTCL", "BRCA", "IPMN", "LUSC", "LUAD", "LNC", 
                  "LAM", "LIHC", "HB", "GC", "GIST", "CRC", "RCC", "PDAC", "EC", "OVCA", 
                  "HGSC", "CESC", "PRAD", "MIBC", "OS", "CSCC", "SKCM")
colors <- c("#915824", "#AD8031", "#C9A063", "#C5A97F", "#CBBE82", "#BAB64D", "#C7C458", 
            "#907F28", "#D8E7A3", "#E4E653", "#D79E63", "#EDC824", "#DF7561", "#91302C", 
            "#B66C6A", "#C693C1", "#C7A1C9", "#B8A3B9", "#B8AECA", "#55A7AE", "#89CDC9", 
            "#649DD4", "#C1E4F8", "#D79E63", "#81B464", "#506248", "#BF5099", "#D89AB1", 
            "#F5B8C9", "#F5B2B1", "#488B45", "#EA9D19", "#9F8A83", "#C4CDB5", "#54626A")

#颜色
cancer_data <- data.frame(
  type = factor(cancer_types, levels = cancer_types), # 保持顺序
  color = colors
)

data <- read.table(file = "/ST_merge_res0.8_metadatainfor_BdyCoreBud.txt",sep = "\t",header = T)
aa<-data[!duplicated(data[,c(1)]),]
cancer_slice <- as.data.frame(table(aa$cancer))

df <- merge(cancer_slice,cancer_data,by.x = "Var1",by.y = "type")

#癌症类型顺序
cancer_order <- c("GBM", "HNSCC", "HN-AS", "DIPG", "PN", "PCNSL", "DLBCL", "LGACC", "OSCC", 
                  "NPC", "TC", "DSRCT", "PTCL", "BRCA", "IPMN", "LUSC", "LUAD", "LNC", 
                  "LAM", "LIHC", "HB", "GC", "GIST", "CRC", "RCC", "PDAC", "EC", "OVCA", 
                  "HGSC", "CESC", "PRAD", "MIBC", "OS", "CSCC", "SKCM")

df$Var1 <- factor(df$Var1, levels = cancer_order)
df <- df[order(df$Var1), ]

par(mar = c(12, 5, 4, 2) + 0.1)  

bar_positions <- gap.barplot1(df, 
                              y.cols = 2,
                              brk.type = "zigzag",
                              col = as.character(df$color),
                              brk.srt = 400,
                              brk.size = 0.5,
                              brk.lwd = 2,
                              max.fold = 5,
                              ratio = 0.5,
                              gap.width = 0.5,
                              cex.error = 1)


axis(1, at = bar_positions, labels = as.character(df$Var1), 
     las = 2, cex.axis = 0.7, tick = FALSE, line = -0.5)

title(main = "Number of slices", 
      xlab = "",  
      ylab = "Frequency / Count",
      cex.main = 1.2)

#可以添加参考线
abline(h = 0, col = "gray50")



