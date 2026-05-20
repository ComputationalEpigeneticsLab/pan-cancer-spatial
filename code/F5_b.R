###cellchat分开统计
library(CellChat)
library(Seurat)
#library(SeuratData)
library(tidyverse)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(reshape2)

bdy_step1CAF_divi<-read.delim('/F5/bdy_step1CAF_diviLR.txt',stringsAsFactors=F,check.names=F)
bdy_step1TAM_divi<-read.delim('/F5/bdy_step1TAM_diviLR.txt',stringsAsFactors=F,check.names=F)

matrix_CAF_divi<-read.delim('/F5/matrix_CAF_diviLR.txt',stringsAsFactors = F,check.names = F)
matrix_TAM_divi<-read.delim('/F5/matrix_TAM_diviLR.txt',stringsAsFactors = F,check.names = F)
#sum(matrix_CAF_divi$LR_num)

bdy_step1CAF_divi$L_R<-'source'
bdy_step1CAF_divi$L_R[which(bdy_step1CAF_divi$target=='CAF')]<-'target'
bdy_CAFpath<-as.data.frame.array(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$L_R))
bdy_CAFpath$sum<-bdy_CAFpath$source+bdy_CAFpath$target

bdy_step1TAM_divi$L_R<-'source'
bdy_step1TAM_divi$L_R[which(bdy_step1TAM_divi$target=='TAM')]<-'target'
bdy_TAMpath<-as.data.frame.array(table(bdy_step1TAM_divi$pathway_name,bdy_step1TAM_divi$L_R))
bdy_TAMpath$sum<-bdy_TAMpath$source+bdy_TAMpath$target

matrix_CAF_divi$L_R<-'source'
matrix_CAF_divi$L_R[which(matrix_CAF_divi$target=='CAF')]<-'target'
matrix_CAFpath<-as.data.frame.array(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$L_R))
matrix_CAFpath$sum<-matrix_CAFpath$source+matrix_CAFpath$target

matrix_TAM_divi$L_R<-'source'
matrix_TAM_divi$L_R[which(matrix_TAM_divi$target=='TAM')]<-'target'
matrix_TAMpath<-as.data.frame.array(table(matrix_TAM_divi$pathway_name,matrix_TAM_divi$L_R))
matrix_TAMpath$sum<-matrix_TAMpath$source+matrix_TAMpath$target


dir_pic<-'/F5/Fig5b'
####matrix_TAM_divi#######################################################################################
bar_data<-as.data.frame(table(matrix_TAM_divi$pathway_name,matrix_TAM_divi$L_R))
bar_data<-bar_data[which(bar_data$Freq>0),]
bar_data<-bar_data[bar_data$Var1%in%rownames(matrix_TAMpath)[order(matrix_TAMpath$sum,decreasing = T)][1:10],]
bar_data$Var1<-factor(bar_data$Var1,levels = rownames(matrix_TAMpath)[order(matrix_TAMpath$sum,decreasing = T)][1:10])
bar_data$logFreq<-log10(bar_data$Freq)
bar_data<-bar_data[order(bar_data$Var1),]
bar_data<-bar_data[order(bar_data$Var2),]
sum_fre<-bar_data$Freq[bar_data$Var2%in%'source']+bar_data$Freq[bar_data$Var2%in%'target']
bar_data$lab<-bar_data$Freq/sum_fre
bar_data$lab<-paste0(round(bar_data$lab * 100, 1), "%")
bar_data$color[which(bar_data$Var2== "source")]<-'#fcbe6b'
bar_data$color[which(bar_data$Var2== "target")]<-'#1578b4'

p_dot<-ggplot(bar_data, aes(x=Var1, y=Var2,size=logFreq, fill = color)) +
  geom_point(
    #shape = 21,  
    color = bar_data$color) + # 使用shape = 21画圈
  scale_color_manual(values = c("source" = '#fcbe6b',
                                "target" = '#1578b4') )+ #设置填充颜色
  #scale_size_continuous(range = c(2, 8))+
  theme_minimal()
print(p_dot)
pdf(paste0(dir_pic,'/dot_matrix_TAM_divi1.pdf'),width = 6,height = 2)##饼图 width = 6,height = 5
print(p_dot)
dev.off()

matrix_TAM <- bar_data
matrix_TAM$resource <- "matrix_TAM"

####matrix_CAF_divi#######################################################################################
bar_data<-as.data.frame(table(matrix_CAF_divi$pathway_name,matrix_CAF_divi$L_R))
bar_data<-bar_data[which(bar_data$Freq>0),]
bar_data<-bar_data[bar_data$Var1%in%rownames(matrix_CAFpath)[order(matrix_CAFpath$sum,decreasing = T)][1:10],]
bar_data$Var1<-factor(bar_data$Var1,levels = rownames(matrix_CAFpath)[order(matrix_CAFpath$sum,decreasing = T)][1:10])
bar_data$logFreq<-log10(bar_data$Freq)
bar_data<-bar_data[order(bar_data$Var1),]
bar_data<-bar_data[order(bar_data$Var2),]
sum_fre<-bar_data$Freq[bar_data$Var2%in%'source']+bar_data$Freq[bar_data$Var2%in%'target']
bar_data$lab<-bar_data$Freq/sum_fre
bar_data$lab<-paste0(round(bar_data$lab * 100, 1), "%")
bar_data$color[which(bar_data$Var2== "source")]<-'#fcbe6b'
bar_data$color[which(bar_data$Var2== "target")]<-'#1578b4'

p_dot<-ggplot(bar_data, aes(x=Var1, y=Var2,size=logFreq, fill = color)) +
  geom_point(
    #shape = 21,  
    color = bar_data$color) + # 使用shape = 21画圈
  scale_color_manual(values = c("source" = '#fcbe6b',
                                "target" = '#1578b4') )+ #设置填充颜色
  #scale_size_continuous(range = c(2, 8))+
  theme_minimal()
print(p_dot)
pdf(paste0(dir_pic,'/dot_matrix_CAF_divi1.pdf'),width = 6,height = 2)##饼图 width = 6,height = 5
print(p_dot)
dev.off()

matrix_CAF <- bar_data
matrix_CAF$resource <- "matrix_CAF"

####bdy_step1CAF_divi#######################################################################################
bar_data<-as.data.frame(table(bdy_step1CAF_divi$pathway_name,bdy_step1CAF_divi$L_R))
bar_data<-bar_data[which(bar_data$Freq>0),]
bar_data<-bar_data[bar_data$Var1%in%rownames(bdy_CAFpath)[order(bdy_CAFpath$sum,decreasing = T)][1:10],]
bar_data$Var1<-factor(bar_data$Var1,levels = rownames(bdy_CAFpath)[order(bdy_CAFpath$sum,decreasing = T)][1:10])
bar_data$logFreq<-log10(bar_data$Freq)
bar_data<-bar_data[order(bar_data$Var1),]
bar_data<-bar_data[order(bar_data$Var2),]
sum_fre<-bar_data$Freq[bar_data$Var2%in%'source']+bar_data$Freq[bar_data$Var2%in%'target']
bar_data$lab<-bar_data$Freq/sum_fre
bar_data$lab<-paste0(round(bar_data$lab * 100, 1), "%")
bar_data$color[which(bar_data$Var2== "source")]<-'#fcbe6b'
bar_data$color[which(bar_data$Var2== "target")]<-'#1578b4'

p_dot<-ggplot(bar_data, aes(x=Var1, y=Var2,size=logFreq, fill = color)) +
  geom_point(
    #shape = 21,  
    color = bar_data$color) + # 使用shape = 21画圈
  scale_color_manual(values = c("source" = '#fcbe6b',
                                "target" = '#1578b4') )+ #设置填充颜色
  #scale_size_continuous(range = c(2, 8))+
  theme_minimal()
print(p_dot)
pdf(paste0(dir_pic,'/dot_bdy_step1CAF_divi1.pdf'),width = 6,height = 2)##饼图 width = 6,height = 5
print(p_dot)
dev.off()

bdy_step1CAF <- bar_data
bdy_step1CAF$resource <- "bdy_step1CAF"

####bdy_step1TAM_divi#######################################################################################
bar_data<-as.data.frame(table(bdy_step1TAM_divi$pathway_name,bdy_step1TAM_divi$L_R))
bar_data<-bar_data[which(bar_data$Freq>0),]
bar_data<-bar_data[bar_data$Var1%in%rownames(bdy_TAMpath)[order(bdy_TAMpath$sum,decreasing = T)][1:10],]
bar_data$Var1<-factor(bar_data$Var1,levels = rownames(bdy_TAMpath)[order(bdy_TAMpath$sum,decreasing = T)][1:10])
bar_data$logFreq<-log10(bar_data$Freq)
bar_data<-bar_data[order(bar_data$Var1),]
bar_data<-bar_data[order(bar_data$Var2),]
sum_fre<-bar_data$Freq[bar_data$Var2%in%'source']+bar_data$Freq[bar_data$Var2%in%'target']
bar_data$lab<-bar_data$Freq/sum_fre
bar_data$lab<-paste0(round(bar_data$lab * 100, 1), "%")
bar_data$color[which(bar_data$Var2== "source")]<-'#fcbe6b'
bar_data$color[which(bar_data$Var2== "target")]<-'#1578b4'

p_dot<-ggplot(bar_data, aes(x=Var1, y=Var2,size=logFreq, fill = color)) +
  geom_point(
    #shape = 21,  
    color = bar_data$color) + # 使用shape = 21画圈
  scale_color_manual(values = c("source" = '#fcbe6b',
                                "target" = '#1578b4') )+ #设置填充颜色
  #scale_size_continuous(range = c(2, 8))+
  theme_minimal()
print(p_dot)
pdf(paste0(dir_pic,'/dot_bdy_step1TAM_divi1.pdf'),width = 6,height = 2)##饼图 width = 6,height = 5
print(p_dot)
dev.off()

bdy_step1TAM <- bar_data
bdy_step1TAM$resource <- "bdy_step1TAM"

#cbind
bar_data <- rbind(matrix_TAM,matrix_CAF)
bar_data <- rbind(bar_data,bdy_step1TAM)
bar_data <- rbind(bar_data,bdy_step1CAF)

bar_data$name <- paste(bar_data$resource, bar_data$Var2, sep="_")

p_dot<-ggplot(bar_data, aes(x=Var1, y=name,size=logFreq, fill = color)) +
  geom_point(
    #shape = 21,  
    color = bar_data$color) + # 使用shape = 21画圈
  scale_color_manual(values = c("source" = '#fcbe6b',
                                "target" = '#1578b4') )+ #设置填充颜色
  #scale_size_continuous(range = c(2, 8))+
  theme_minimal()
print(p_dot)
pdf(paste0(dir_pic,'/dot.pdf'),width = 7.5,height = 2.5)##饼图 width = 6,height = 5
print(p_dot)
dev.off()
