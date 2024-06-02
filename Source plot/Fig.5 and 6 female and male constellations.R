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

setwd("~/.")

##for female
load('female_constellations.RData')

## set label colors
l3<-sort(unique(nodes$label3))
set.seed(1)
rainb<-sample(rainbow(15),15)
rainb[4]<-'darkred'
rainb[7]<-'skyblue3'
rainb[11]<-'yellowgreen'
rainb[12]<-'#CDBE70'
centralcol<-'seagreen4'
for (i in 1:15) {
  nodes$label4[nodes$label3==l3[i]]=rainb[i]
}

## plots
core<-c()
n_max <- 10  
for (n in 1:n_max) {
  var_name <- paste("p", n, sep="")  
  e<-edgess[edgess$group==n,1:2]
  e2<-merge(e,e,by.x = c("source","target"),by.y = c("target","source")) 
  e1<-setdiff(e,e2) ##uni_direct
  e2<- e2[!duplicated(t(apply(e2, 1, sort))),] ##bi_direct
  ## weighted edges
  for (i in 1:nrow(e2)) {
    e21<-edgess[edgess$source==e2[i,1]&edgess$target==e2[i,2],3]
    e22<-edgess[edgess$source==e2[i,2]&edgess$target==e2[i,1],3]
    e2$wei[i]<-ifelse(e21>=0.01&e22>=0.01,3,ifelse(e21<0.01&e22<0.01,1,2))
  }
  e2$directed<-FALSE
  if(nrow(e1)!=0) {
    for (i in 1:nrow(e1)) {
      e1$wei[i]<-edgess[edgess$source==e1[i,1]&edgess$target==e1[i,2],3]
      e1$wei[i]<-ifelse(e1$wei[i]>=0.01,2,1)
    }
    e1$directed<-TRUE
  }
  ed<-rbind(e1,e2)
  net_pc<-graph_from_data_frame(
    d=ed,vertices=nodes[nodes$group==n,])
  graph_pc<-as_tbl_graph(net_pc)
  no<- nodes[nodes$group==n,]
  subm<-m2[no$label,no$label]
  subd<-melt(subm,na.rm = TRUE)
  subd<-subd[subd$value>0,]
  colnames(subd)<-c('source','target','weight')
  subg<-graph_from_data_frame(subd,directed = TRUE)
  ## page_rank to get core
  core[n]<-names(sort(page_rank(subg)$vector,decreasing = TRUE)[1])
  
  ##get network
  p<-ggraph(graph_pc, layout = 'kk')+
    geom_edge_link(aes(filter= directed==TRUE,edge_width=wei),color="gray",
                   arrow = arrow(length = unit(3.5, 'mm')),
                   end_cap=circle(0.35,'inches'),start_cap =circle(0.35,'inches'),show.legend = FALSE)+
    geom_edge_link(color=centralcol,aes(filter= directed==FALSE,edge_width=wei),
                   show.legend = FALSE)+
    geom_node_label(aes(label = name), family = "serif", fontface = "bold",show.legend = FALSE,cex=5,color = nodes$label4[nodes$group==n])+
    geom_node_text(aes(label = (label2)), family = "sans",cex=2, fontface = "italic",repel = TRUE,box.padding=0.8,segment.color=NA)+
    theme(panel.background = element_rect(fill='white'), plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))+
    scale_edge_width(range=c(0,max(ed$wei)/3))+
    scale_x_continuous(expand = c(0, 0.15))+
    scale_y_continuous(expand = c(0, 0.15))
  if(n==9){
    p<-ggraph(graph_pc, layout = 'kk')+
      geom_edge_link(aes(filter= directed==TRUE,edge_width=wei),color="gray",
                     arrow = arrow(length = unit(3.5, 'mm')),
                     end_cap=circle(0.35,'inches'),start_cap =circle(0.35,'inches'),show.legend = FALSE)+
      geom_edge_link(color=centralcol,aes(filter= directed==FALSE,edge_width=wei),
                     show.legend = FALSE)+
      geom_node_label(aes(label = name), family = "serif", fontface = "bold",show.legend = FALSE,cex=5,color = nodes$label4[nodes$group==n])+
      geom_node_text(aes(label = (label2)), family = "sans",cex=2, fontface = "italic",repel = TRUE,box.padding=0.8,segment.color=NA)+
      theme(legend.position = c(0.5,0.075),legend.key = element_rect(fill='white'),panel.background = element_rect(fill='white'),
            plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))+
      guides(color=guide_legend(override.aes = list(color='black')))+
      scale_edge_width(name='',range=c(0,max(ed$wei)/3))+
      scale_x_continuous(expand = c(0, 0.3))+
      scale_y_continuous(expand = c(0.2, 0.4))
  }
  assign(var_name, p)  
}
## combined sub figures
prow1<-plot_grid(p8,p10,p4,rel_widths =c(1,1.5,1.5), labels = "auto", hjust = 0, vjust = 1,label_size = 20,nrow=1)
prow2<-plot_grid(p5,p1,p7,p2,labels = c("d","e",'f','g'), hjust = 0, vjust = 1,label_size = 20,nrow=1)
prow3<-plot_grid(p6,p3,p9,rel_widths = c(1,2.5,0.5), nrow = 1,labels = c('h','i','j'),label_size = 20)
prowsum<-plot_grid(prow1, prow2,prow3, nrow=3,rel_widths =c(1.5,1,1))
ggsave(filename = "results/figures/fig5_female constellations.pdf", plot = prowsum, width = 25, height =20)

rm(list = ls())
##for male
load('male_constellations.RData')

## set label colors
l3<-sort(unique(nodes$label3))
set.seed(1)
rainb<-sample(rainbow(15),15)
rainb[4]<-'darkred'
rainb[7]<-'skyblue3'
rainb[11]<-'yellowgreen'
rainb[12]<-'#CDBE70'
centralcol<-'seagreen4'
for (i in 1:15) {
  nodes$label4[nodes$label3==l3[i]]=rainb[i]
}

## plots
core<-c()
n_max <- 9  
for (n in 1:n_max) {
  var_name <- paste("p", n, sep="")  
  e<-edgess[edgess$group==n,1:2]
  e2<-merge(e,e,by.x = c("source","target"),by.y = c("target","source")) 
  e1<-setdiff(e,e2) ##uni_direct
  e2<- e2[!duplicated(t(apply(e2, 1, sort))),] ##bi_direct
  ## weighted edges
  for (i in 1:nrow(e2)) {
    e21<-edgess[edgess$source==e2[i,1]&edgess$target==e2[i,2],3]
    e22<-edgess[edgess$source==e2[i,2]&edgess$target==e2[i,1],3]
    e2$wei[i]<-ifelse(e21>=0.01&e22>=0.01,3,ifelse(e21<0.01&e22<0.01,1,2))
  }
  e2$directed<-FALSE
  if(nrow(e1)!=0) {
    for (i in 1:nrow(e1)) {
      e1$wei[i]<-edgess[edgess$source==e1[i,1]&edgess$target==e1[i,2],3]
      e1$wei[i]<-ifelse(e1$wei[i]>=0.01,2,1)
    }
    e1$directed<-TRUE
  }
  ed<-rbind(e1,e2)
  net_pc<-graph_from_data_frame(
    d=ed,vertices=nodes[nodes$group==n,])
  graph_pc<-as_tbl_graph(net_pc)
  no<- nodes[nodes$group==n,]
  subm<-m2[no$label,no$label]
  subd<-melt(subm,na.rm = TRUE)
  subd<-subd[subd$value>0,]
  colnames(subd)<-c('source','target','weight')
  subg<-graph_from_data_frame(subd,directed = TRUE)
  ## page_rank to get core
  core[n]<-names(sort(page_rank(subg)$vector,decreasing = TRUE)[1])
  
  ##get network
  p<-ggraph(graph_pc, layout = 'kk')+
    geom_edge_link(aes(filter= directed==TRUE,edge_width=wei),color="gray",
                   arrow = arrow(length = unit(3.5, 'mm')),
                   end_cap=circle(0.35,'inches'),start_cap =circle(0.35,'inches'),show.legend = FALSE)+
    geom_edge_link(color=centralcol,aes(filter= directed==FALSE,edge_width=wei),
                   show.legend = FALSE)+
    geom_node_label(aes(label = name), family = "serif", fontface = "bold",show.legend = FALSE,cex=5,color = nodes$label4[nodes$group==n])+
    geom_node_text(aes(label = (label2)), family = "sans",cex=2, fontface = "italic",repel = TRUE,box.padding=0.8,segment.color=NA)+
    theme(panel.background = element_rect(fill='white'), plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))+
    scale_edge_width(range=c(0,max(ed$wei)/3))+
    scale_x_continuous(expand = c(0, 0.15))+
    scale_y_continuous(expand = c(0, 0.15))
  if(n==6){
    p<-ggraph(graph_pc, layout = 'kk')+
      geom_edge_link(aes(filter= directed==TRUE,edge_width=wei),color="gray",
                     arrow = arrow(length = unit(3.5, 'mm')),
                     end_cap=circle(0.35,'inches'),start_cap =circle(0.35,'inches'),show.legend = FALSE)+
      geom_edge_link(color=centralcol,aes(filter= directed==FALSE,edge_width=wei),
                     show.legend = FALSE)+
      geom_node_label(aes(label = name), family = "serif", fontface = "bold",show.legend = FALSE,cex=5,color = nodes$label4[nodes$group==n])+
      geom_node_text(aes(label = (label2)), family = "sans",cex=2, fontface = "italic",repel = TRUE,box.padding=0.8,segment.color=NA)+
      theme(legend.position = c(0.5,0.075),legend.key = element_rect(fill='white'),panel.background = element_rect(fill='white'),
            plot.margin = margin(t = 0, r = 0, b = 0,l = 0,unit = "cm"))+
      guides(color=guide_legend(override.aes = list(color='black')))+
      scale_edge_width(name='',range=c(0,max(ed$wei)/3))+
      scale_x_continuous(expand = c(0, 0.3))+
      scale_y_continuous(expand = c(0.2, 0.4))
  }
  assign(var_name, p)  
}
## combined sub figures
prow1<-plot_grid(p3,p8,p7,rel_widths =c(1,1.5,1.5), labels = "auto", hjust = 0, vjust = 1,label_size = 20,nrow=1)
prow2<-plot_grid(p1,p2,p5,labels = c("d","e",'f'),rel_widths = c(1.5,1,1), hjust = 0, vjust = 1,label_size = 20,nrow=1)
prow3<-plot_grid(p9,p4,p6,rel_widths = c(1,2.5,1), nrow = 1,labels = c('g','h','i'),label_size = 20)
prowsum<-plot_grid(prow1, prow2,prow3, nrow=3,rel_widths =c(1.5,1,1))
ggsave(filename = "results/figures/fig6_male constellations.pdf", plot = prowsum, width = 25, height =20)
