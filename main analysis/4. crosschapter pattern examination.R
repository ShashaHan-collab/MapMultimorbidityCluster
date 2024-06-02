rm(list = ls())
library(ggplot2)
library(reshape2)
library(devtools)
library(dplyr)
library(grid)
library(cowplot)
library(igraph)
library(mlbench)  
library(psych)
library(RColorBrewer)
library(circlize)
library(dendextend)
library(openxlsx)
library(readr)

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

## multimorbidity progress across the current disease classifications
EstBiPosCross <- EstBiPos %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

EstBiPosCross_sh<-EstBiPosCross%>%
  mutate(name1 = substr(pair, 1, 3),name2 = substr(pair, 5, 7))

EstBiPosNegCross <- EstBiPosNeg %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

EstBiPosNeg_sh<-EstBiPosNegCross%>%
  mutate(name1 = substr(pair, 1, 3),name2 = substr(pair, 5, 7))

##bi-direction
g<-EstBiPosCross_sh[,c(7,8)]
g<-graph_from_data_frame(g,directed = FALSE)

obs<-EstBiPosCross_sh%>%
  group_by(code1,code2)%>%
  summarize(count = n())
cross_obs<-sum(obs[obs$code1!=obs$code2,3])

x<-data.frame(colnames(m2))
colnames(x)<-'code'
y<-x%>%mutate(clu = ifelse(substr(code, 1, 2) %in% c('D0','D1','D2'), 'C', 
                           ifelse(substr(code, 1, 2) %in% c('H9'), 'H9', 
                                  ifelse(substr(code, 1, 2) %in% c('A0','A4'), 'B', substr(code,1,1)))))

rownames(y)<-y$code
record_p<-matrix(0,length(unique(y$clu)),length(unique(y$clu)))
rownames(record_p)<-unique(y$clu)
colnames(record_p)<-unique(y$clu)
n<-10000
times<-0
p_cross<-0
while(times<n){
  ##shuffle
  g1<-rewire(g,with = keeping_degseq(niter = vcount(g) * 10))
  g2<- igraph::as_data_frame(g1)
  colnames(g2)<-c('newname1','newname2')
  g3<-g2%>%
    mutate(newcode1 = substr(newname1, 1, 1),newcode2 = substr(newname2, 1, 1))%>%
    mutate(newcode1 = ifelse(substr(newname1, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname1, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname1, 1, 2) %in% c('A0','A4'), 'B', newcode1))),
           newcode2 = ifelse(substr(newname2, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname2, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname2, 1, 2) %in% c('A0','A4'), 'B', newcode2))))%>%
    mutate(chapter = paste0(newcode1,newcode2))%>%
    group_by(newcode1,newcode2)%>%
    summarize(count = n())
  
  cross_randomitem<-sum(g3[g3$newcode1!=g3$newcode2,3])
  p_cross<-ifelse(cross_randomitem>=cross_obs,p_cross+1,p_cross)
  for (i in 1:length(unique(y$clu))) {
    for (j in i:length(unique(y$clu))){
      code1<-rownames(record_p)[i]
      code2<-colnames(record_p)[j]
      ob<-obs$count[obs$code1==code1&obs$code2==code2]
      it1<-g3$count[g3$newcode1==code1&g3$newcode2==code2]
      it2<-g3$count[g3$newcode1==code2&g3$newcode2==code1]
      ob<-ifelse(length(ob)!=0,ob,0)
      if(ob==0){
        next
      }
      it1<-ifelse(length(it1)!=0,it1,0)
      it2<-ifelse(length(it2)!=0,it2,0)
      ifelse(i!=j,it<-it1+it2,it<-it1)
      if(it>=ob){
        record_p[i,j]<-record_p[i,j]+1
      }
    }
  }
  times<-times+1
}
record_p<-record_p/n
record_pp<-melt(record_p)
colnames(record_pp)<-colnames(obs)
contr<-merge(obs,record_pp,by.x = c("code1","code2"),by.y = c("code1","code2"))

##uni-direction
G<-EstBiPosNeg_sh[,c(7,8)]
G<-graph_from_data_frame(G,directed = TRUE)

obs2<-EstBiPosNeg_sh%>%
  group_by(code1,code2)%>%
  summarize(count = n())

cross_obs2<-sum(obs2[obs2$code1!=obs2$code2,3])
times=0
p_cross2<-0

record_p2<-matrix(0,length(unique(y$clu)),length(unique(y$clu)))
rownames(record_p2)<-unique(y$clu)
colnames(record_p2)<-unique(y$clu)

while(times<n){
  g1<-rewire(G,with = keeping_degseq(niter = vcount(G) * 10))
  g2<- igraph::as_data_frame(g1)
  colnames(g2)<-c('newname1','newname2')
  g3<-g2%>%
    mutate(newcode1 = substr(newname1, 1, 1),newcode2 = substr(newname2, 1, 1))%>%
    mutate(newcode1 = ifelse(substr(newname1, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname1, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname1, 1, 2) %in% c('A0','A4'), 'B', newcode1))),
           newcode2 = ifelse(substr(newname2, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname2, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname2, 1, 2) %in% c('A0','A4'), 'B', newcode2))))%>%
    mutate(chapter = paste0(newcode1,newcode2))%>%
    group_by(newcode1,newcode2)%>%
    summarize(count = n())
  
  cross_randomitem<-sum(g3[g3$newcode1!=g3$newcode2,3])
  p_cross2<-ifelse(cross_randomitem>=cross_obs2,p_cross2+1,p_cross2)
  for (i in 1:length(unique(y$clu))) {
    for (j in 1:length(unique(y$clu))){
      code1<-rownames(record_p2)[i]
      code2<-colnames(record_p2)[j]
      ob<-obs2$count[obs2$code1==code1&obs2$code2==code2]
      it1<-g3$count[g3$newcode1==code1&g3$newcode2==code2]
      ob<-ifelse(length(ob)!=0,ob,0)
      if(ob==0){
        next
      }
      it1<-ifelse(length(it1)!=0,it1,0)
      it<-it1
      if(it>=ob){
        record_p2[i,j]<-record_p2[i,j]+1
      }
    }
  }
  times<-times+1
}
record_p2<-record_p2/n
record_pp2<-melt(record_p2)
colnames(record_pp2)<-colnames(obs2)
contr2<-merge(obs2,record_pp2,by.x = c("code1","code2"),by.y = c("code1","code2"))

## for male
rm(list = ls())
dat<-read.csv('./results_male/pairwise_estimates.csv')
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

## multimorbidity progress across the current disease classifications
EstBiPosCross <- EstBiPos %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

EstBiPosCross_sh<-EstBiPosCross%>%
  mutate(name1 = substr(pair, 1, 3),name2 = substr(pair, 5, 7))

EstBiPosNegCross <- EstBiPosNeg %>%
  mutate(code1 = substr(pair, 1, 1),code2 = substr(pair, 5, 5))%>%
  mutate(code1 = ifelse(substr(pair, 1, 2) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 1, 2) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 1, 2) %in% c('A0','A4'), 'B', code1))),
         code2 = ifelse(substr(pair, 5, 6) %in% c('D0','D1','D2'), 'C', 
                        ifelse(substr(pair, 5, 6) %in% c('H9'), 'H9', 
                               ifelse(substr(pair, 5, 6) %in% c('A0','A4'), 'B', code2))))%>%
  mutate(chapter = paste0(code1,code2))

EstBiPosNeg_sh<-EstBiPosNegCross%>%
  mutate(name1 = substr(pair, 1, 3),name2 = substr(pair, 5, 7))

##bi-direction
g<-EstBiPosCross_sh[,c(7,8)]
g<-graph_from_data_frame(g,directed = FALSE)

obs<-EstBiPosCross_sh%>%
  group_by(code1,code2)%>%
  summarize(count = n())
cross_obs<-sum(obs[obs$code1!=obs$code2,3])

x<-data.frame(colnames(m2))
colnames(x)<-'code'
y<-x%>%mutate(clu = ifelse(substr(code, 1, 2) %in% c('D0','D1','D2'), 'C', 
                           ifelse(substr(code, 1, 2) %in% c('H9'), 'H9', 
                                  ifelse(substr(code, 1, 2) %in% c('A0','A4'), 'B', substr(code,1,1)))))

rownames(y)<-y$code
record_p<-matrix(0,length(unique(y$clu)),length(unique(y$clu)))
rownames(record_p)<-unique(y$clu)
colnames(record_p)<-unique(y$clu)
n<-10000
times<-0
p_cross<-0
while(times<n){
  ##shuffle
  g1<-rewire(g,with = keeping_degseq(niter = vcount(g) * 10))
  g2<- igraph::as_data_frame(g1)
  colnames(g2)<-c('newname1','newname2')
  g3<-g2%>%
    mutate(newcode1 = substr(newname1, 1, 1),newcode2 = substr(newname2, 1, 1))%>%
    mutate(newcode1 = ifelse(substr(newname1, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname1, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname1, 1, 2) %in% c('A0','A4'), 'B', newcode1))),
           newcode2 = ifelse(substr(newname2, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname2, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname2, 1, 2) %in% c('A0','A4'), 'B', newcode2))))%>%
    mutate(chapter = paste0(newcode1,newcode2))%>%
    group_by(newcode1,newcode2)%>%
    summarize(count = n())
  
  cross_randomitem<-sum(g3[g3$newcode1!=g3$newcode2,3])
  p_cross<-ifelse(cross_randomitem>=cross_obs,p_cross+1,p_cross)
  for (i in 1:length(unique(y$clu))) {
    for (j in i:length(unique(y$clu))){
      code1<-rownames(record_p)[i]
      code2<-colnames(record_p)[j]
      ob<-obs$count[obs$code1==code1&obs$code2==code2]
      it1<-g3$count[g3$newcode1==code1&g3$newcode2==code2]
      it2<-g3$count[g3$newcode1==code2&g3$newcode2==code1]
      ob<-ifelse(length(ob)!=0,ob,0)
      if(ob==0){
        next
      }
      it1<-ifelse(length(it1)!=0,it1,0)
      it2<-ifelse(length(it2)!=0,it2,0)
      ifelse(i!=j,it<-it1+it2,it<-it1)
      if(it>=ob){
        record_p[i,j]<-record_p[i,j]+1
      }
    }
  }
  times<-times+1
}
record_p<-record_p/n
record_pp<-melt(record_p)
colnames(record_pp)<-colnames(obs)
contr<-merge(obs,record_pp,by.x = c("code1","code2"),by.y = c("code1","code2"))

##uni-direction
G<-EstBiPosNeg_sh[,c(7,8)]
G<-graph_from_data_frame(G,directed = TRUE)

obs2<-EstBiPosNeg_sh%>%
  group_by(code1,code2)%>%
  summarize(count = n())

cross_obs2<-sum(obs2[obs2$code1!=obs2$code2,3])
times=0
p_cross2<-0

record_p2<-matrix(0,length(unique(y$clu)),length(unique(y$clu)))
rownames(record_p2)<-unique(y$clu)
colnames(record_p2)<-unique(y$clu)

while(times<n){
  g1<-rewire(G,with = keeping_degseq(niter = vcount(G) * 10))
  g2<- igraph::as_data_frame(g1)
  colnames(g2)<-c('newname1','newname2')
  g3<-g2%>%
    mutate(newcode1 = substr(newname1, 1, 1),newcode2 = substr(newname2, 1, 1))%>%
    mutate(newcode1 = ifelse(substr(newname1, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname1, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname1, 1, 2) %in% c('A0','A4'), 'B', newcode1))),
           newcode2 = ifelse(substr(newname2, 1, 2) %in% c('D0','D1','D2'), 'C', 
                             ifelse(substr(newname2, 1, 2) %in% c('H9'), 'H9', 
                                    ifelse(substr(newname2, 1, 2) %in% c('A0','A4'), 'B', newcode2))))%>%
    mutate(chapter = paste0(newcode1,newcode2))%>%
    group_by(newcode1,newcode2)%>%
    summarize(count = n())
  
  cross_randomitem<-sum(g3[g3$newcode1!=g3$newcode2,3])
  p_cross2<-ifelse(cross_randomitem>=cross_obs2,p_cross2+1,p_cross2)
  for (i in 1:length(unique(y$clu))) {
    for (j in 1:length(unique(y$clu))){
      code1<-rownames(record_p2)[i]
      code2<-colnames(record_p2)[j]
      ob<-obs2$count[obs2$code1==code1&obs2$code2==code2]
      it1<-g3$count[g3$newcode1==code1&g3$newcode2==code2]
      ob<-ifelse(length(ob)!=0,ob,0)
      if(ob==0){
        next
      }
      it1<-ifelse(length(it1)!=0,it1,0)
      it<-it1
      if(it>=ob){
        record_p2[i,j]<-record_p2[i,j]+1
      }
    }
  }
  times<-times+1
}
record_p2<-record_p2/n
record_pp2<-melt(record_p2)
colnames(record_pp2)<-colnames(obs2)
contr2<-merge(obs2,record_pp2,by.x = c("code1","code2"),by.y = c("code1","code2"))
