rm(list = ls())
library(ggplot2)
library(reshape2)
library(devtools)
library(dplyr)
library(grid)
library(cowplot)
library(stringr)
library(igraph)
library(mlbench)    
library(psych)
library(RColorBrewer)
library(circlize)
library(dendextend)
library(openxlsx)
library(readr)
library(ggraph)
library(tidygraph)
library(RcmdrMisc)
setwd("D:/RA/heatmap")
dat<-read.csv('main_analysis_Aug14/main_analysis_Aug14/results_female/pairwise_estimates_imp1.csv')
dat1<-read.csv('main_analysis_Aug14/main_analysis_Aug14/results_female/pairwise_estimates_imp1_alpha_0.05_2.csv')
dat <- anti_join(x = dat, y = dat1, by = c('expo','outcome'))
dat <- rbind(dat, dat1)
dat<-dat[,c(2:3,5:7)]
dat<-dat[order(dat[,2]),]
dat<-dat[order(dat[,1]),]
dat[is.na(dat)]<-0
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
dat$mean[(dat$low95ci<=0 & dat$up95ci>=0)]=0
dat<-dat[,1:3]
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
dat_z<-melt(m2,na.rm = TRUE)
name<-colnames(m2)
data <- read_tsv("ICD10_coding.tsv")
rownames(data)<-data$coding
name2<-data[name,]$meaning
M2<-ifelse(m2<=0,0,1)
diag(M2)<-0
MM<-M2
for (i in 1:nrow(M2)) {
  for (j in i+1:nrow(M2)) {
    if(j>nrow(M2)){break}
    MM[i,j]<-M2[i,j]+M2[j,i]
    MM[j,i]<-M2[i,j]+M2[j,i]
    
  }
}
diag(MM)<-2
k<-pamk(MM,criterion="multiasw",krange=6:20,critout=FALSE,seed = 1e+08)$nc
set.seed(1)
R<-KMeans(MM,centers=k,iter.max=1000, num.seeds=1000)
ata<-R$cluster
positive<-dat_z[dat_z$value>0,][,1:3]
positive$cl=0
for (i in 1:nrow(positive)) {
  if(ata[positive[i,1]]==ata[positive[i,2]]){positive$cl[i]=ata[positive[i,1]]}
}
edgess<-positive[positive$cl>0,]
nodes<-data.frame(names(ata),ata)
colnames(nodes)<-c("label",'group')
colnames(edgess)<-c('source','target','weight','group')
label2=substr(name2,5,200)
nodes$label2=str_wrap(label2,width = 20)
nodes$label3=as.factor(substr(name2,1,1))
nodes$label3[nodes$label3=='B']='A'
nodes$label3[nodes$label %in% c('D05','D12','D17','D22','D23','D25')]='C'
nodes$label3[nodes$label=='H91']='B'



isolate<-function(x){
  g <- graph_from_data_frame(x, directed = FALSE)
  components <- clusters(g)
  isolated_pairs <- which(components$csize == 2)
  isolated_nodes <- V(g)[components$membership %in% isolated_pairs]
  return(isolated_nodes)
}
modified<-c()
for (i in 1:9) {
  x<-edgess[edgess$group==i,1:2]
  modified<-c(modified,names(isolate(x)))
}
rows_with_specific_elements <- apply(edgess, 1, function(row) {
  any(sapply(modified, function(element) {
    grepl(pattern = element, x = row)
  }))
})

edgess<-edgess[!rows_with_specific_elements,]


for (i in 1:nrow(nodes)) {
  com<-nodes[i,2]
  if(sum(nodes[i,1]==edgess[edgess$group==com,1])+sum(nodes[i,1]==edgess[edgess$group==com,2])==0 | nodes[i,1]%in% modified )
  {addgroup<-c(positive[positive$Var1==nodes[i,1],2],positive[positive$Var2==nodes[i,1],1])
  a=(nodes[nodes$label %in% addgroup,2])
  nodes[i,2]<-as.numeric(names(table(a))[table(a) == max(table(a))])[1]
  gn<-nodes$label[nodes$group==nodes[i,2]]
  newedge1<-positive[positive$Var1==nodes[i,1]&positive$Var2 %in% gn,]
  newedge2<-positive[positive$Var2==nodes[i,1]&positive$Var1 %in% gn,]
  newedge<-rbind(newedge1,newedge2)
  colnames(newedge)<-c('source','target','weight','group')
  newedge$group=nodes[i,2]
  edgess<-rbind(edgess,newedge)
  }
}


save.image('female_constellations.RData')

rm(list = ls())
setwd("D:/RA/heatmap")
dat<-read.csv('main_analysis_Aug14/main_analysis_Aug14/results_male/pairwise_estimates_imp1.csv')
dat1<-read.csv('main_analysis_Aug14/main_analysis_Aug14/results_male/pairwise_estimates_imp1_alpha_0.05_2.csv')
dat <- anti_join(x = dat, y = dat1, by = c('expo','outcome'))
dat <- rbind(dat, dat1)
dat<-dat[,c(2:3,5:7)]
dat<-dat[order(dat[,2]),]
dat<-dat[order(dat[,1]),]
dat[is.na(dat)]<-0
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
dat$mean[(dat$low95ci<=0 & dat$up95ci>=0)]=0
dat<-dat[,1:3]
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
dat_z<-melt(m2,na.rm = TRUE)

name<-colnames(m2)
data <- read_tsv("D:/RA/heatmap/ICD10_coding.tsv")
rownames(data)<-data$coding
#intersect(colnames(m2),data$coding)
name2<-data[name,]$meaning
name2[160]<-'U07 New diseases of uncertain etiology or emergency use'

M2<-ifelse(m2<=0,0,1)
diag(M2)<-0
MM<-M2
for (i in 1:nrow(M2)) {
  for (j in i+1:nrow(M2)) {
    if(j>nrow(M2)){break}
    MM[i,j]<-M2[i,j]+M2[j,i]
    MM[j,i]<-M2[i,j]+M2[j,i]
    
  }
}
diag(MM)<-2
k<-pamk(MM,criterion="multiasw",krange=6:20,critout=FALSE,seed = 1e+08)$nc
set.seed(1)
R<-KMeans(MM,centers=9,iter.max=1000, num.seeds=1000)
ata<-R$cluster
positive<-dat_z[dat_z$value>0,][,1:3]
positive$cl=0
#ata[c('J33','J34')]=8

for (i in 1:nrow(positive)) {
  if(ata[positive[i,1]]==ata[positive[i,2]]){positive$cl[i]=ata[positive[i,1]]}
  
}

edgess<-positive[positive$cl>0,]

nodes<-data.frame(names(ata),ata)
colnames(nodes)<-c("label",'group')
colnames(edgess)<-c('source','target','weight','group')
label2=substr(name2,5,200)
nodes$label2=str_wrap(label2,width = 20)
nodes$label3=as.factor(substr(name2,1,1))
nodes$label3[nodes$label3=='B']='A'
nodes$label3[nodes$label %in% c('D05','D12','D17','D22','D23','D25')]='C'
nodes$label3[nodes$label=='H91']='B'
isolate<-function(x){
  g <- graph_from_data_frame(x, directed = FALSE)
  components <- clusters(g)
  isolated_pairs <- which(components$csize == 2)
  isolated_nodes <- V(g)[components$membership %in% isolated_pairs]
  return(isolated_nodes)
}
modified<-c()
for (i in 1:k) {
  x<-edgess[edgess$group==i,1:2]
  modified<-c(modified,names(isolate(x)))
}
rows_with_specific_elements <- apply(edgess, 1, function(row) {
  any(sapply(modified, function(element) {
    grepl(pattern = element, x = row)
  }))
})

edgess<-edgess[!rows_with_specific_elements,]


for (i in 1:nrow(nodes)) {
  com<-nodes[i,2]
  if(sum(nodes[i,1]==edgess[edgess$group==com,1])+sum(nodes[i,1]==edgess[edgess$group==com,2])==0 | nodes[i,1]%in% modified )
  {addgroup<-c(positive[positive$Var1==nodes[i,1],2],positive[positive$Var2==nodes[i,1],1])
  a=(nodes[nodes$label %in% addgroup,2])
  nodes[i,2]<-as.numeric(names(table(a))[table(a) == max(table(a))])[1]
  gn<-nodes$label[nodes$group==nodes[i,2]]
  newedge1<-positive[positive$Var1==nodes[i,1]&positive$Var2 %in% gn,]
  newedge2<-positive[positive$Var2==nodes[i,1]&positive$Var1 %in% gn,]
  newedge<-rbind(newedge1,newedge2)
  colnames(newedge)<-c('source','target','weight','group')
  newedge$group=nodes[i,2]
  edgess<-rbind(edgess,newedge)
  }
}
