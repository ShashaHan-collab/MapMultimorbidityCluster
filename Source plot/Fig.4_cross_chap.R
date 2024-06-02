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
library(plotly)

setwd("~/.")
###############################################
##for female
dat<-read.csv('./results_female/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
Results<-m2

## find the pairs of bi-directional effects with --> + and  <-- +
EstResults <- Results
EstBiPos <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 < 0.01] <- NA
  est2[est2 < 0.01] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPos <- rbind(EstBiPos, EstBi.sub)
  #}
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

## find the pairs of bi-directional effects with --> - and  <-- -
EstResults <- Results
EstBiNeg <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 > -0.01] <- NA
  est2[est2 > -0.01] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiNeg <- rbind(EstBiNeg, EstBi.sub)
  #}
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}


## find the pairs of bi-directional effects with --> + and  <-- - or  --> - and  <-- +
## and reverse  '--> - and  <-- +' to '--> + and  <-- - or  --> -'
EstResults <- Results
EstBiPosNeg <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  ## --> + and  <-- - 
  
  est1[est1 < 0.01] <- NA
  est2[est2 > 0] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosNeg <- rbind(EstBiPosNeg, EstBi.sub)
  #}
  
  ## --> -  and  <-- +
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 > 0] <- NA
  est2[est2 < 0.01] <- NA
  
  ## reverse the order
  EstBi.sub <- cbind(est2,est1)%>% 
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  ## reverse the order
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(rownames(EstBi.sub),'\u2014',code_i))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosNeg <- rbind(EstBiPosNeg, EstBi.sub)
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

## find the pairs of bi-directional effects with --> + and  <-- 0 or  --> 0 and  <-- +
## and reverse  '--> 0 and  <-- +' to '--> + and  <-- 0'
EstResults <- Results
EstBiPosZero <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  ## --> + and  <-- - or zeros
  
  est1[est1 < 0.01] <- NA
  est2[est2 !=0] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosZero <- rbind(EstBiPosZero, EstBi.sub)
  #}
  
  ## --> - or zeros and  <-- +
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 !=0] <- NA
  est2[est2 < 0.01] <- NA
  
  ## reverse the order
  EstBi.sub <- cbind(est2,est1)%>% 
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  ## reverse the order
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(rownames(EstBi.sub),'\u2014',code_i))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosZero <- rbind(EstBiPosZero, EstBi.sub)
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}


###############################
## check the cross-chapter pairs
EstBiPosCross <- EstBiPos %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))


EstBiPosNegCross <- EstBiPosNeg %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

##bi-directional progress
BIPO<-EstBiPosCross%>%
  group_by(chapter)%>%
  summarize(count = n())%>%
  arrange(-count)
posplit<-strsplit(BIPO$chapter,'')
biposplit<-data.frame(matrix(,nrow=nrow(BIPO),ncol=3))
colnames(biposplit)<-c('source','target','weight')
for (i in 1:nrow(BIPO)) {
  te<-posplit[[i]]
  if(length(te)==2){
    biposplit$source[i]=te[1]
    biposplit$target[i]=te[2]}
  else{
    numposi<-which(te==9)
    if(numposi==2){
      biposplit$source[i]=paste0(te[1],te[2])
      biposplit$target[i]=te[3]
    }else{
      biposplit$source[i]=te[1]
      biposplit$target[i]=paste0(te[2],te[3])
    }
  }
}
biposplit$weight<-BIPO$count
biposplit<-rbind(biposplit,c('O','O',0))
biposplit$weight<-as.numeric(biposplit$weight)
chapter<-c('A00-B99',
           'C00-D48',
           'D50-D89',
           'E00-E90',
           'F00-F99',
           'G00-G99',
           'H00-H59',
           'H60-H95',
           'I00-I99',
           'J00-J99',
           'K00-K93',
           'L00-L99',
           'M00-M99',
           'N00-N99',
           'O00–O99')
name<-data.frame(c(biposplit$source,biposplit$target))
ordername<-unique(sort(name[,1]))
matchtable<-data.frame(ordername,chapter)
rownames(matchtable)<-ordername
for (i in 1:nrow(biposplit)) {
  temp1<-biposplit[i,1]
  temp2<-biposplit[i,2]
  biposplit[i,1]<-paste0('',matchtable[temp1,2])
  biposplit[i,2]<-paste0('',matchtable[temp2,2])
}
edges3<-biposplit
nodes3<-data.frame(unique(c(edges3$source,edges3$target)))
colnames(nodes3)<-"label"


## uni-directional
BIPO<-EstBiPosNegCross%>%
  group_by(chapter)%>%
  summarize(count = n())%>%
  arrange(-count)
posplit<-strsplit(BIPO$chapter,'')
biposplit<-data.frame(matrix(,nrow=nrow(BIPO),ncol=3))
colnames(biposplit)<-c('source','target','weight')
for (i in 1:nrow(BIPO)) {
  te<-posplit[[i]]
  if(length(te)==2){
    biposplit$source[i]=te[1]
    biposplit$target[i]=te[2]}
  else{
    numposi<-which(te==9)
    if(numposi==2){
      biposplit$source[i]=paste0(te[1],te[2])
      biposplit$target[i]=te[3]
    }else{
      biposplit$source[i]=te[1]
      biposplit$target[i]=paste0(te[2],te[3])
    }         
  }  
}
biposplit$weight<-BIPO$count
chapter<-c('A00-B99',
           'C00-D48',
           'D50-D89',
           'E00-E90',
           'F00-F99',
           'G00-G99',
           'H00-H59',
           'H60-H95',
           'I00-I99',
           'J00-J99',
           'K00-K93',
           'L00-L99',
           'M00-M99',
           'N00-N99',
           'O00–O99')

name<-data.frame(c(biposplit$source,biposplit$target))
ordername<-unique(sort(name[,1]))
matchtable<-data.frame(ordername,chapter)
rownames(matchtable)<-ordername
for (i in 1:nrow(biposplit)) {
  temp1<-biposplit[i,1]
  temp2<-biposplit[i,2]
  biposplit[i,1]<-paste0('',matchtable[temp1,2])
  biposplit[i,2]<-paste0('',matchtable[temp2,2]) 
}
edges4<-biposplit
nodes4<-data.frame(unique(c(edges4$source,edges4$target)))
colnames(nodes4)<-"label"

###############################################
### for males
Results <-read.csv('./results_male/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
Results<-m2

## find the pairs of bi-directional effects with --> + and  <-- +
EstResults <- Results
EstBiPos <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 < 0.01] <- NA
  est2[est2 < 0.01] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPos <- rbind(EstBiPos, EstBi.sub)
  #}
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

## find the pairs of bi-directional effects with --> - and  <-- -
EstResults <- Results
EstBiNeg <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 > -0.01] <- NA
  est2[est2 > -0.01] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiNeg <- rbind(EstBiNeg, EstBi.sub)
  #}
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

## find the pairs of bi-directional effects with --> + and  <-- - or  --> - and  <-- +
## and reverse  '--> - and  <-- +' to '--> + and  <-- - or  --> -'
EstResults <- Results
EstBiPosNeg <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  ## --> + and  <-- - 
  
  est1[est1 < 0.01] <- NA
  est2[est2 > 0] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosNeg <- rbind(EstBiPosNeg, EstBi.sub)
  #}
  
  ## --> -  and  <-- +
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 > 0] <- NA
  est2[est2 < 0.01] <- NA
  
  ## reverse the order
  EstBi.sub <- cbind(est2,est1)%>% 
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  ## reverse the order
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(rownames(EstBi.sub),'\u2014',code_i))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosNeg <- rbind(EstBiPosNeg, EstBi.sub)
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

## find the pairs of bi-directional effects with --> + and  <-- 0 or  --> 0 and  <-- +
## and reverse  '--> 0 and  <-- +' to '--> + and  <-- 0'
EstResults <- Results
EstBiPosZero <- NULL

while (sum(!is.na(EstResults), na.rm = TRUE)){
  i = 1
  code_i <- colnames(EstResults)[i]
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  ## --> + and  <-- - or zeros
  
  est1[est1 < 0.01] <- NA
  est2[est2 !=0] <- NA
  
  EstBi.sub <- cbind(est1,est2)%>%
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(code_i,'\u2014',rownames(EstBi.sub)))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosZero <- rbind(EstBiPosZero, EstBi.sub)
  #}
  
  ## --> - or zeros and  <-- +
  est1 <- EstResults[i, ] %>% unlist()
  est2 <- EstResults[, i]
  
  est1[est1 !=0] <- NA
  est2[est2 < 0.01] <- NA
  
  ## reverse the order
  EstBi.sub <- cbind(est2,est1)%>% 
    data.frame(.)%>%na.omit(.)
  names(EstBi.sub) <- c('forward','backward')
  
  ## reverse the order
  EstBi.sub <- EstBi.sub %>%
    mutate(pair = paste0(rownames(EstBi.sub),'\u2014',code_i))
  
  #if (dim(EstBi.sub)[1] >0){
  EstBiPosZero <- rbind(EstBiPosZero, EstBi.sub)
  
  EstResults <- EstResults[, -i] %>% data.frame()
  EstResults <- EstResults[-i, ]
}

###############################
## check the cross-chapter pairs
EstBiPosCross <- EstBiPos %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

EstBiPosNegCross <- EstBiPosNeg %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

##bi-directional
BIPO<-EstBiPosCross%>%
  group_by(chapter)%>%
  summarize(count = n())%>%
  arrange(-count)
posplit<-strsplit(BIPO$chapter,'')
biposplit<-data.frame(matrix(,nrow=nrow(BIPO),ncol=3))
colnames(biposplit)<-c('source','target','weight')
for (i in 1:nrow(BIPO)) {
  te<-posplit[[i]]
  if(length(te)==2){
    biposplit$source[i]=te[1]
    biposplit$target[i]=te[2]}
  else{
    numposi<-which(te==9)
    if(numposi==2){
      biposplit$source[i]=paste0(te[1],te[2])
      biposplit$target[i]=te[3]
    }else{
      biposplit$source[i]=te[1]
      biposplit$target[i]=paste0(te[2],te[3])
    }   
  } 
}
biposplit$weight<-BIPO$count
#biposplit<-rbind(biposplit,c('O','O',0))
biposplit$weight<-as.numeric(biposplit$weight)
chapter<-c('A00-B99',
           'C00-D48',
           'D50-D89',
           'E00-E90',
           'F00-F99',
           'G00-G99',
           'H00-H59',
           'H60-H95',
           'I00-I99',
           'J00-J99',
           'K00-K93',
           'L00-L99',
           'M00-M99',
           'N00-N99',
           'U00–U49')

name<-data.frame(c(biposplit$source,biposplit$target))
ordername<-unique(sort(name[,1]))
matchtable<-data.frame(ordername,chapter)
rownames(matchtable)<-ordername
for (i in 1:nrow(biposplit)) {
  temp1<-biposplit[i,1]
  temp2<-biposplit[i,2]
  biposplit[i,1]<-paste0('',matchtable[temp1,2])
  biposplit[i,2]<-paste0('',matchtable[temp2,2])
}
edges1<-biposplit
nodes1<-data.frame(unique(c(edges1$source,edges1$target)))
colnames(nodes1)<-"label"

## uni-directional
BIPO<-EstBiPosNegCross%>%
  group_by(chapter)%>%
  summarize(count = n())%>%
  arrange(-count)
posplit<-strsplit(BIPO$chapter,'')
biposplit<-data.frame(matrix(,nrow=nrow(BIPO),ncol=3))
colnames(biposplit)<-c('source','target','weight')
for (i in 1:nrow(BIPO)) {
  te<-posplit[[i]]
  if(length(te)==2){
    biposplit$source[i]=te[1]
    biposplit$target[i]=te[2]}
  else{
    numposi<-which(te==9)
    if(numposi==2){
      biposplit$source[i]=paste0(te[1],te[2])
      biposplit$target[i]=te[3]
    }else{
      biposplit$source[i]=te[1]
      biposplit$target[i]=paste0(te[2],te[3])
    }   
  } 
}
biposplit$weight<-BIPO$count
chapter<-c('A00-B99',
           'C00-D48',
           'D50-D89',
           'E00-E90',
           'F00-F99',
           'G00-G99',
           'H00-H59',
           'H60-H95',
           'I00-I99',
           'J00-J99',
           'K00-K93',
           'L00-L99',
           'M00-M99',
           'N00-N99',
           'U00–U49')

name<-data.frame(c(biposplit$source,biposplit$target))
ordername<-unique(sort(name[,1]))
matchtable<-data.frame(ordername,chapter)
rownames(matchtable)<-ordername
for (i in 1:nrow(biposplit)) {
  temp1<-biposplit[i,1]
  temp2<-biposplit[i,2]
  biposplit[i,1]<-paste0('',matchtable[temp1,2])
  biposplit[i,2]<-paste0('',matchtable[temp2,2])
  
}
edges2<-biposplit
nodes2<-data.frame(unique(c(edges2$source,edges2$target)))
colnames(nodes2)<-"label"

######################################################################################
##plots
pdf("results/figures/fig4_cross_chap.pdf", width  =8, height = 8)
par(mfrow=c(2,2), mar = c(0, 0, 0, 0) , oma = c(0, 0, 0, 0))

## uni-directional progress for female
circos.par(gap.after = 10)
edge41<-edges4[edges4$source==edges4$target,]
edge42<-edges4[edges4$source!=edges4$target,]
E4<-rbind(edge41,edge42)
chordDiagram(link.sort='asis',self.link=1,E4,order =sort(nodes4$label),col="#002fa7",grid.col <-'white',transparency = 0.25
             ,directional=1,direction.type="diffHeight+arrows", diffHeight = convert_height(0, "mm"),
             annotationTrack = c("grid", "axis"), link.arr.type = "triangle",link.arr.length=0.1,link.arr.col='gray'       
)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data(
    "xlim", sector.index = si, track.index = 1
  )
  ylim = get.cell.meta.data(
    "ylim", sector.index = si, track.index = 1
  )
  circos.text(
    mean(xlim), mean(ylim), si, sector.index = si,
    track.index = 1,  facing = "bending.inside",
    niceFacing = TRUE, col = "black",cex = 0.5
  )
}
circos.clear()
mtext('a',side=3,line=-2,adj=0,cex=1,font = 2)

## uni-directional progress for male
circos.par(gap.after = 10)
edge21<-edges2[edges2$source==edges2$target,]
edge22<-edges2[edges2$source!=edges2$target,]
E2<-rbind(edge21,edge22)
chordDiagram(link.sort='asis',self.link=1,E2,order =sort(nodes2$label),col="#002fa7",grid.col <-'white',transparency = 0.25
             ,directional=1,direction.type="diffHeight+arrows", diffHeight = convert_height(0, "mm"),
             annotationTrack = c("grid", "axis"), link.arr.type = "triangle",link.arr.length=0.1,link.arr.col='gray'    
)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data(
    "xlim", sector.index = si, track.index = 1
  )
  ylim = get.cell.meta.data(
    "ylim", sector.index = si, track.index = 1
  )
  circos.text(
    mean(xlim), mean(ylim), si, sector.index = si,
    track.index = 1,  facing = "bending.inside",
    niceFacing = TRUE, col = "black",cex = 0.5
  )
}
mtext('b',side=3,line=-2,adj=0,cex=1,font = 2)
circos.clear()

## bi-directional progress for female
circos.par(gap.after = 10)
edge31<-edges3[edges3$source==edges3$target,]
edge32<-edges3[edges3$source!=edges3$target,]
E3<-rbind(edge31,edge32)
chordDiagram(link.sort='asis',self.link=1,E3,order =sort(nodes3$label),col="#B3495C",grid.col <-'white',transparency = 0.25
             ,directional=2,direction.type="diffHeight+arrows", diffHeight = convert_height(0, "mm"),
             annotationTrack = c("grid", "axis"), link.arr.type = "triangle",link.arr.length=0.1,link.arr.col='gray'             
)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data(
    "xlim", sector.index = si, track.index = 1
  )
  ylim = get.cell.meta.data(
    "ylim", sector.index = si, track.index = 1
  )
  circos.text(
    mean(xlim), mean(ylim), si, sector.index = si,
    track.index = 1,  facing = "bending.inside",
    niceFacing = TRUE, col = "black",cex = 0.5
  )
}
mtext('c',side=3,line=-2,adj=0,cex=1,font = 2)
circos.clear()

## bi-directional progress for male
circos.par(gap.after = 10)
edge11<-edges1[edges1$source==edges1$target,]
edge12<-edges1[edges1$source!=edges1$target,]
E1<-rbind(edge11,edge12)
chordDiagram(link.sort='asis',self.link=1,E1,order =sort(nodes1$label),col="#B3495C",grid.col <-'white',transparency = 0.25
             ,directional=2,direction.type="diffHeight+arrows", diffHeight = convert_height(0, "mm"),
             annotationTrack = c("grid", "axis"), link.arr.type = "triangle",link.arr.length=0.1,link.arr.col='gray'    
)
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data(
    "xlim", sector.index = si, track.index = 1
  )
  ylim = get.cell.meta.data(
    "ylim", sector.index = si, track.index = 1
  )
  circos.text(
    mean(xlim), mean(ylim), si, sector.index = si,
    track.index = 1,  facing = "bending.inside",
    niceFacing = TRUE, col = "black",cex = 0.5
  )
}
mtext('d',side=3,line=-2,adj=0,cex=1,font = 2)
circos.clear()
dev.off()
