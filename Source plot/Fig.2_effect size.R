rm(list = ls())
library(ggplot2)
library(reshape2)
library(devtools)
library(dplyr)
library(grid)
library(cowplot)
library(ComplexHeatmap)
library(igraph)
library(mlbench)   
library(psych)
library(RColorBrewer)
library(circlize)
library(dendextend)
library(openxlsx)
library(readr)

setwd("~/.")
# for females
dat<-read.csv('./results_female/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
name<-colnames(m2)

data <- read_tsv("ICD10_coding.tsv")
rownames(data)<-data$coding
name2<-data[name,]$meaning

mat1<-m2
colnames(mat1)<-name2
rownames(mat1)<-name2

col_fun1 = colorRamp2(c(-1, -0.05,-0.01,0,0.01,0.05,1), c("navy","blue","lightskyblue3",'white', "orange","red", "firebrick"))
pdf("results/figures/fig2_effect size_a.pdf", width  =7.5, height = 9)
par(mfrow = c(1,1), mar = c(0, 0, 0, 0) , oma = c(0, 0, 0, 0))
circos.par(gap.after = c(15),track.margin=c(0,0))
circos.heatmap(mat1, col = col_fun1,na.col = "black",dend.side = "inside",rownames.side = "outside",rownames.cex = 0.25,
               track.height = 0.45,cluster = FALSE, cell.lwd = 0.05)

circos.text(-3,165, labels = name[1], cex = 0.25, facing = "inside")
circos.text(-3,-5, labels = name[length(name)], cex = 0.25, facing = "inside")
arrows(x0=0.45,y0=0.05,x1=0.05,y1=0.005,length = 0.03)


circos.clear()
text(labels='a',x=-0.72,y=0.825,cex = 0.75,font=2)

lgd = Legend(at<-c(-1,-0.1,-0.01,0,0.01,0.1,1),title = "", labels_gp = gpar(fontsize = 3),break_dist=20,tick_length= unit(0, "mm"),
             grid_width = unit(2, "mm"), col_fun = col_fun1,legend_height = unit(20, "mm"))
pushViewport(viewport(x = unit(0.35, "npc"), y = unit(0.65, "npc"), just = c("left", "top")))
grid.draw(lgd)
popViewport()
dev.off()


## for males
rm(list = ls())
dat<-read.csv('./results_male/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
label <- !stringr::str_detect(colnames(m2), '[P]|[R-T]|[V-Z]\\d{2}')
m2 <- m2[, label]
m2 <- m2[label, ]
name<-colnames(m2)

data <- read_tsv("ICD10_coding.tsv")
rownames(data)<-data$coding
name2<-data[name,]$meaning

mat1<-m2
colnames(mat1)<-name2
rownames(mat1)<-name2

col_fun1 = colorRamp2(c(-1, -0.05,-0.01,0,0.01,0.05,1), c("navy","blue","lightskyblue3",'white', "orange","red", "firebrick"))
pdf("results/figures/fig2_effect size_b.pdf", width  =7.5, height = 9)
par(mfrow = c(1,1), mar = c(0, 0, 0, 0) , oma = c(0, 0, 0, 0))
circos.par(gap.after = c(15),track.margin=c(0,0))
circos.heatmap(mat1, col = col_fun1,na.col = "black",dend.side = "inside",rownames.side = "outside",rownames.cex = 0.25,
               track.height = 0.45,cluster = FALSE, cell.lwd = 0.05)
circos.text(-3,165, labels = name[1], cex = 0.25, facing = "inside")
circos.text(-3,-5, labels = name[length(name)], cex = 0.25, facing = "inside")
arrows(x0=0.45,y0=0.05,x1=0.05,y1=0.005,length = 0.03)
circos.clear()
text(labels='b',x=-0.72,y=0.825,cex = 0.75,font=2)
dev.off()
