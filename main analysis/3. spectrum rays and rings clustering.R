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
## for females

##ray
dat<-read.csv('./results_female/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
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
merge <- hclust_res$merge
height <- hclust_res$height

check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  
    return(row)
  } else {
    return(x)  
  }
}

last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
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
  }
}

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

#rings
m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-(m4)
pv3<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
clusters<-cutree(pv3,h=1.260)
k<-max(clusters)
pau<-pv3[["edges"]][["au"]]
hclust_res <- as.hclust(pv3$hclust)
merge <- hclust_res$merge
height <- hclust_res$height


check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  
    return(row)
  } else {
    return(x) 
  }
}

last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
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
    
  }
  
}

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
### for males
#ray
rm(list = ls())

dat<-read.csv('./results_male/pairwise_estimates.csv')
colnames(dat)<-c('expo','outcome','mean','low95ci','up95ci')
m2 <- acast(dat, expo ~ outcome, value.var = "mean")
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
pv2<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
t2<-Sys.time()
t2-t1
load('1014_2.RData')

clusters<-(cutree(pv2,h=0.661))
k<-max(clusters)
pau<-pv2[["edges"]][["au"]]

hclust_res <- as.hclust(pv2$hclust)
merge <- hclust_res$merge
height <- hclust_res$height


check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  
    return(row)
  } else {
    return(x)  
  }
}


last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
 
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
    
  }
  
}


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


#rings
m4<-ifelse(m2>0,1,ifelse(m2<0,-1,0))
diag(m4)<-1
m5<-(m4)
pv4<-pvclust(m5,method.hclust = 'ward.D2',nboot = 2000,iseed = 1)
#load('1014_4.RData')
clusters<-cutree(pv4,h=1.260)
k<-max(clusters)
pau<-pv4[["edges"]][["au"]]
hclust_res <- as.hclust(pv4$hclust)

merge <- hclust_res$merge
height <- hclust_res$height

check_value <- function(x, data) {
  if (x %in% data) {
    row <- which(data == x, arr.ind = TRUE)[,1]  
    return(row)
  } else {
    return(x)  
  }
}


last_merge <- function(cluster_assignment, cluster_number) {
  merge<-merge[1:(length(cluster_assignment)-k),]
  members <- which(cluster_assignment == cluster_number)
  if(length(members)==1){
    return('NULL')
    
  }
  
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


