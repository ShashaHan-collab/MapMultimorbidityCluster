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
library(pvclust)
setwd("~/.")
###############################################
###female
#ray
dat<-read.csv('main_analysis_Aug14/results_female/pairwise_estimates_imp1.csv')
dat1<-read.csv('main_analysis_Aug14/results_female/pairwise_estimates_imp1_alpha_0.05_2.csv')
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
## remove icd-10 codes starting with z
label <- !stringr::str_detect(colnames(m2), '[P]|[R-T]|[V-Z]\\d{2}')
m2 <- m2[, label]
m2 <- m2[label, ]

name<-colnames(m2)
data <- read_tsv("ICD10_coding.tsv")
rownames(data)<-data$coding
name2<-data[name,]$meaning


fun_rxy = function(x,y){
  x_bar = mean(x)
  y_bar = mean(y)
  fenzi = sum((x-x_bar)*(y-y_bar))
  fenmu = sqrt(sum((x-x_bar)^2)*sum((y-y_bar)^2))
  r_xy = fenzi/fenmu
  return(r_xy)
}

m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-t(m4)
t1<-Sys.time()
pv1<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
t2<-Sys.time()
#load('1014_1.RData')
clusters<-(cutree(pv1,h=0.661))
k<-max(clusters)
pau<-pv1[["edges"]][["au"]]

hclust_res <- as.hclust(pv1$hclust)
# 获取合并矩阵和高度
merge <- hclust_res$merge
height <- hclust_res$height

# 检查函数
check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  # 找到值所在的行
    return(row)
  } else {
    return(x)  # 如果不在矩阵中，返回原值
  }
}


last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  # 找出属于该类的所有样本
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
  # 找出所有合并事件中涉及这些成员的最大高度
  involved_merges <- which(apply(-(merge), 1, function(x) any(x %in% members)))
  historys<-c()
  ori<-involved_merges
  while (1) {
    
    for(me in involved_merges){
      if( me %in% merge[involved_merges,]| me %in% historys[1:(length(historys)-1)]) involved_merges<-involved_merges[involved_merges!=me]
      if(length(involved_merges)==1) {
        return((involved_merges))
        break
      }
    }
    if(length(involved_merges)==1) break
    historys<-c(historys,involved_merges,ori)
    historys<-order(unique(historys))
    involved_merges<- unique(sapply(involved_merges, check_value, merge))
    involved_merges<-involved_merges[order(involved_merges)]
    #involved_merges<-unique(which(matrix(merge %in% involved_merges,ncol=2), arr.ind = TRUE)[,1])
    
  }
  
}

# 对每个类应用函数
p<-c()
edges <- unlist(lapply(unique(clusters), function(x) last_merge(clusters, x)))
for ( j in 1:length(edges)) {
  if(edges[j]=='NULL'){
    p[j]<-1
  }else{
    p[j]<-pau[as.integer(edges[j])]
  }
}
names(p)<-edges
modi<-which(p<0.40)
modi_edge<-as.integer(names(modi))
modi_clu<-(clusters[clusters==modi])
child_cluster<-merge[modi_edge,]

find_nodes<-function(e){
  subs<-merge[e,]
  record<-c()
  if(all(subs<0)){
    return(subs)
  }
  else if(all(subs>0)){
    a1<-find_nodes(subs[1])
    a2<-find_nodes(subs[2])
    return(c(a1,a2))
  }
  else{
    b1<-subs[subs<0]
    b2<-subs[subs>0]
    b2<-find_nodes(b2)
    return(c(b1,b2))
  }
}
labels<-pv1[["hclust"]][["labels"]]

get_new_clu<-function(item,pv,modi_clu,child_cluster,k){
  if(all(pau[child_cluster]>=0.4)) {
    la<-labels[-find_nodes(item)]
    vec<-clusters[la]
    vec1=labels[-find_nodes(child_cluster[1])]
    vec2=labels[-find_nodes(child_cluster[2])]
    newc<-modi_clu
    newc[vec2]<-k+1
    return(newc)
  }
  else if(any(pau[child_cluster]>=0.4)){
    yesorno<-pau[child_cluster]!=0
    names(yesorno)<-child_cluster
    for (i in 1:2) {
      sub<-yesorno[i]
      nam<-names(yesorno)[i]
      if(sub) {
        e<-as.integer(nam)
        ord=-find_nodes(e)
        newc1<-labels[ord]
        newc1<-modi_clu[newc1]
        newc1[]<-k+1
        k=k+1
      }else{
        e<-as.integer(nam)
        child_cluster2<-merge[e,]
        modi_clu2<- -find_nodes(e)
        modi_clu2<-labels[modi_clu2]
        modi_clu2<-modi_clu[modi_clu2]
        newc2<-get_new_clu(e,pv,modi_clu2,child_cluster2,k)
        k=k+1
        
      }
    }
    newc<-c(newc1,newc2)
    return(newc)
  }
  
}
newc<-get_new_clu(modi_edge,pv1,modi_clu,child_cluster,k)
clusters[names(modi_clu)]<-newc
save(list ='clusters', file = "female_ray.RData")





#ring

m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-(m4)
pv3<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
#load('1014_3.RData')
clusters<-cutree(pv3,h=1.260)
k<-max(clusters)
pau<-pv3[["edges"]][["au"]]
hclust_res <- as.hclust(pv3$hclust)
# 获取合并矩阵和高度
merge <- hclust_res$merge
height <- hclust_res$height

# 检查函数
check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  # 找到值所在的行
    return(row)
  } else {
    return(x)  # 如果不在矩阵中，返回原值
  }
}

last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  # 找出属于该类的所有样本
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
  # 找出所有合并事件中涉及这些成员的最大高度
  involved_merges <- which(apply(-(merge), 1, function(x) any(x %in% members)))
  historys<-c()
  ori<-involved_merges
  while (1) {
    
    for(me in involved_merges){
      if( me %in% merge[involved_merges,]| me %in% historys[1:(length(historys)-1)]) involved_merges<-involved_merges[involved_merges!=me]
      if(length(involved_merges)==1) {
        return((involved_merges))
        break
      }
    }
    if(length(involved_merges)==1) break
    historys<-c(historys,involved_merges,ori)
    historys<-order(unique(historys))
    involved_merges<- unique(sapply(involved_merges, check_value, merge))
    involved_merges<-involved_merges[order(involved_merges)]
    #involved_merges<-unique(which(matrix(merge %in% involved_merges,ncol=2), arr.ind = TRUE)[,1])
    
  }
  
}

# 对每个类应用函数
p<-c()
edges <- unlist(lapply(unique(clusters), function(x) last_merge(clusters, x)))
for ( j in 1:length(edges)) {
  if(edges[j]=='NULL'){
    p[j]<-1
  }else{
    p[j]<-pau[as.integer(edges[j])]
  }
}
names(p)<-edges
modi<-which(p<0.40)
modi_edge<-as.integer(names(modi))
modi_clu<-(clusters[clusters==modi])
child_cluster<-merge[modi_edge,]


find_nodes<-function(e){
  subs<-merge[e,]
  record<-c()
  if(all(subs<0)){
    return(subs)
  }
  else if(all(subs>0)){
    a1<-find_nodes(subs[1])
    a2<-find_nodes(subs[2])
    return(c(a1,a2))
  }
  else{
    b1<-subs[subs<0]
    b2<-subs[subs>0]
    b2<-find_nodes(b2)
    return(c(b1,b2))
  }
}
labels<-pv3[["hclust"]][["labels"]]
get_new_clu<-function(item,pv,modi_clu,child_cluster,k){
  if(all(pau[child_cluster]>=0.4)) {
    la<-labels[-find_nodes(item)]
    vec<-clusters[la]
    vec1=labels[-find_nodes(child_cluster[1])]
    vec2=labels[-find_nodes(child_cluster[2])]
    newc<-modi_clu
    newc[vec2]<-k+1
    return(newc)
  }
  else if(any(pau[child_cluster]>=0.4)){
    yesorno<-pau[child_cluster]!=0
    names(yesorno)<-child_cluster
    for (i in 1:2) {
      sub<-yesorno[i]
      nam<-names(yesorno)[i]
      if(sub) {
        e<-as.integer(nam)
        ord=-find_nodes(e)
        newc1<-labels[ord]
        newc1<-modi_clu[newc1]
        newc1[]<-k+1
        k=k+1
      }else{
        e<-as.integer(nam)
        child_cluster2<-merge[e,]
        modi_clu2<- -find_nodes(e)
        modi_clu2<-labels[modi_clu2]
        modi_clu2<-modi_clu[modi_clu2]
        newc2<-get_new_clu(e,pv,modi_clu2,child_cluster2,k)
        k=k+1
        
      }
    }
    newc<-c(newc1,newc2)
    return(newc)
  }
  
}
newc<-get_new_clu(modi_edge,pv3,modi_clu,child_cluster,k)
clusters[names(modi_clu)]<-newc
save(list ='clusters', file = "female_ring.RData")

###############################################                       
###male
#ray
rm(list = ls())
setwd("D:/RA/heatmap")
dat<-read.csv('main_analysis_Aug14/results_male/pairwise_estimates_imp1.csv')
dat1<-read.csv('main_analysis_Aug14/results_male/pairwise_estimates_imp1_alpha_0.05_2.csv')
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
label <- !stringr::str_detect(colnames(m2), '[P]|[R-T]|[V-Z]\\d{2}')
m2 <- m2[, label]
m2 <- m2[label, ]

fun_rxy = function(x,y){
  x_bar = mean(x)
  y_bar = mean(y)
  fenzi = sum((x-x_bar)*(y-y_bar))
  fenmu = sqrt(sum((x-x_bar)^2)*sum((y-y_bar)^2))
  r_xy = fenzi/fenmu
  return(r_xy)
}
m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-t(m4)
t1<-Sys.time()
pv2<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
t2<-Sys.time()
t2-t1
load('1014_2.RData')

clusters<-(cutree(pv2,h=0.661))
k<-max(clusters)
pau<-pv2[["edges"]][["au"]]

hclust_res <- as.hclust(pv2$hclust)
# 获取合并矩阵和高度
merge <- hclust_res$merge
height <- hclust_res$height

# 检查函数
check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  # 找到值所在的行
    return(row)
  } else {
    return(x)  # 如果不在矩阵中，返回原值
  }
}


last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  # 找出属于该类的所有样本
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
  # 找出所有合并事件中涉及这些成员的最大高度
  involved_merges <- which(apply(-(merge), 1, function(x) any(x %in% members)))
  historys<-c()
  ori<-involved_merges
  while (1) {
    
    for(me in involved_merges){
      if( me %in% merge[involved_merges,]| me %in% historys[1:(length(historys)-1)]) involved_merges<-involved_merges[involved_merges!=me]
      if(length(involved_merges)==1) {
        return((involved_merges))
        break
      }
    }
    if(length(involved_merges)==1) break
    historys<-c(historys,involved_merges,ori)
    historys<-order(unique(historys))
    involved_merges<- unique(sapply(involved_merges, check_value, merge))
    involved_merges<-involved_merges[order(involved_merges)]
    #involved_merges<-unique(which(matrix(merge %in% involved_merges,ncol=2), arr.ind = TRUE)[,1])
    
  }
  
}

# 对每个类应用函数
p<-c()
edges <- unlist(lapply(unique(clusters), function(x) last_merge(clusters, x)))
for ( j in 1:length(edges)) {
  if(edges[j]=='NULL'){
    p[j]<-1
  }else{
    p[j]<-pau[as.integer(edges[j])]
  }
}
names(p)<-edges
modis<-which(p<0.40)
for (nums in 1:length(modis)) {
  modi<-modis[nums]
  modi_edge<-as.integer(names(modi))
  modi_clu<-(clusters[clusters==modi])
  child_cluster<-merge[modi_edge,]
  
  find_nodes<-function(e){
    subs<-merge[e,]
    record<-c()
    if(all(subs<0)){
      return(subs)
    }
    else if(all(subs>0)){
      a1<-find_nodes(subs[1])
      a2<-find_nodes(subs[2])
      return(c(a1,a2))
    }
    else{
      b1<-subs[subs<0]
      b2<-subs[subs>0]
      b2<-find_nodes(b2)
      return(c(b1,b2))
    }
  }
  labels<-pv2[["hclust"]][["labels"]]
  
  get_new_clu<-function(item,pv,modi_clu,child_cluster,k){
    if(all(pau[child_cluster]>=0.4)) {
      la<-labels[-find_nodes(item)]
      vec<-clusters[la]
      vec1=labels[-find_nodes(child_cluster[1])]
      vec2=labels[-find_nodes(child_cluster[2])]
      newc<-modi_clu
      newc[vec2]<-k+1
      return(newc)
    }
    else if(any(pau[child_cluster]>=0.4)){
      yesorno<-pau[child_cluster]!=0
      names(yesorno)<-child_cluster
      for (i in 1:2) {
        sub<-yesorno[i]
        nam<-names(yesorno)[i]
        if(sub) {
          e<-as.integer(nam)
          ord=-find_nodes(e)
          newc1<-labels[ord]
          newc1<-modi_clu[newc1]
          newc1[]<-k+1
          k=k+1
        }else{
          e<-as.integer(nam)
          child_cluster2<-merge[e,]
          modi_clu2<- -find_nodes(e)
          modi_clu2<-labels[modi_clu2]
          modi_clu2<-modi_clu[modi_clu2]
          newc2<-get_new_clu(e,pv,modi_clu2,child_cluster2,k)
          k=k+1
          
        }
      }
      newc<-c(newc1,newc2)
      return(newc)
    }
    
  }
  newc<-get_new_clu(modi_edge,pv2,modi_clu,child_cluster,k)
  clusters[names(modi_clu)]<-newc
  k<-k+1
}

save(list ='clusters', file = "male_ray.RData")


#ring
m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-(m4)
pv4<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
#load('1014_4.RData')
clusters<-cutree(pv4,h=1.260)
k<-max(clusters)
pau<-pv4[["edges"]][["au"]]
hclust_res <- as.hclust(pv4$hclust)
# 获取合并矩阵和高度
merge <- hclust_res$merge
height <- hclust_res$height

# 检查函数
check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  # 找到值所在的行
    return(row)
  } else {
    return(x)  # 如果不在矩阵中，返回原值
  }
}


last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  # 找出属于该类的所有样本
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
  # 找出所有合并事件中涉及这些成员的最大高度
  involved_merges <- which(apply(-(merge), 1, function(x) any(x %in% members)))
  historys<-c()
  ori<-involved_merges
  while (1) {
    
    for(me in involved_merges){
      if( me %in% merge[involved_merges,]| me %in% historys[1:(length(historys)-1)]) involved_merges<-involved_merges[involved_merges!=me]
      if(length(involved_merges)==1) {
        return((involved_merges))
        break
      }
    }
    if(length(involved_merges)==1) break
    historys<-c(historys,involved_merges,ori)
    historys<-order(unique(historys))
    involved_merges<- unique(sapply(involved_merges, check_value, merge))
    involved_merges<-involved_merges[order(involved_merges)]
    #involved_merges<-unique(which(matrix(merge %in% involved_merges,ncol=2), arr.ind = TRUE)[,1])
    
  }
  
}

# 对每个类应用函数
p<-c()
edges <- unlist(lapply(unique(clusters), function(x) last_merge(clusters, x)))
for ( j in 1:length(edges)) {
  if(edges[j]=='NULL'){
    p[j]<-1
  }else{
    p[j]<-pau[as.integer(edges[j])]
  }
}
names(p)<-edges
modi<-which(p<0.40)
modi_edge<-as.integer(names(modi))
modi_clu<-(clusters[clusters==modi])
child_cluster<-merge[modi_edge,]


find_nodes<-function(e){
  subs<-merge[e,]
  record<-c()
  if(all(subs<0)){
    return(subs)
  }
  else if(all(subs>0)){
    a1<-find_nodes(subs[1])
    a2<-find_nodes(subs[2])
    return(c(a1,a2))
  }
  else{
    b1<-subs[subs<0]
    b2<-subs[subs>0]
    b2<-find_nodes(b2)
    return(c(b1,b2))
  }
}
labels<-pv4[["hclust"]][["labels"]]
get_new_clu<-function(item,pv,modi_clu,child_cluster,k){
  if(all(pau[child_cluster]>=0.4)) {
    la<-labels[-find_nodes(item)]
    vec<-clusters[la]
    vec1=labels[-find_nodes(child_cluster[1])]
    vec2=labels[-find_nodes(child_cluster[2])]
    newc<-modi_clu
    newc[vec2]<-k+1
    return(newc)
  }
  else if(any(pau[child_cluster]>=0.4)){
    yesorno<-pau[child_cluster]!=0
    names(yesorno)<-child_cluster
    for (i in 1:2) {
      sub<-yesorno[i]
      nam<-names(yesorno)[i]
      if(sub) {
        e<-as.integer(nam)
        ord=-find_nodes(e)
        newc1<-labels[ord]
        newc1<-modi_clu[newc1]
        newc1[]<-k+1
        k=k+1
      }else{
        e<-as.integer(nam)
        child_cluster2<-merge[e,]
        modi_clu2<- -find_nodes(e)
        modi_clu2<-labels[modi_clu2]
        modi_clu2<-modi_clu[modi_clu2]
        newc2<-get_new_clu(e,pv,modi_clu2,child_cluster2,k)
        k=k+1
        
      }
    }
    newc<-c(newc1,newc2)
    return(newc)
  }
  
}
newc<-get_new_clu(modi_edge,pv4,modi_clu,child_cluster,k)
clusters[names(modi_clu)]<-newc
save(list ='clusters', file = "male_ring.RData")


