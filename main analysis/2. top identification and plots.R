library(data.table)
library(corrplot)
library(factoextra)
library(tidyr)
library(RColorBrewer)
library(rcartocolor)
library(showtext)
library(broom)
library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(dplyr)

library(scales) 
library("ggsci")
library("gridExtra")

const <- list(); const <- within(const, {
  color_vline = 'grey50'
  font_plot = 'roboto'
  vertical_gap = 0.3
  fill_strong_pos = '#BB0021FF'
  fill_weak_pos = '#F39B7FFF'
  fill_strong_neg = '#008280FF'
  fill_weak_neg = '#91D1C2FF'
  text_x_position5 = 1.0
  font_table = 'Roboto Condensed'
  font_countries = 'Roboto Condensed'
  font_xaxis = 'Roboto Condensed'
  fontsize_table = 3.3
})


setwd("~/.")

## for females
## Step 0 : Get results 
Results <- read.csv("./results_female/pairwise_correlation_matrix.csv",header=TRUE)[-1]
## expo
dim = dim(Results)[1]

expo_strong_pos <- Results %>%
  mutate(across(everything(), ~ ifelse(. >= 0.01, 1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_weak_pos <- Results %>%
  mutate(across(everything(), ~ ifelse((. <= 0.01)&(. > 0), 1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_strong_neg <- Results %>%
  mutate(across(everything(), ~ ifelse(. <= -0.01, -1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_weak_neg <- Results %>%
  mutate(across(everything(), ~ ifelse((. >= -0.01)&(. < 0), -1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

Expo <- cbind(expo_strong_pos,expo_weak_pos,expo_strong_neg,expo_weak_neg) %>%
  data.frame(.)%>%
  mutate("Anterior" = rownames(.))%>%
  pivot_longer(., cols = starts_with("exp"), names_to = 'EffectDirection', values_to = 'Proportion')%>%
  mutate(EffectDirection = factor(EffectDirection, levels=c('expo_strong_pos', 'expo_weak_pos',
                                                            'expo_strong_neg','expo_weak_neg')))

order.diseases <- Expo %>%
  filter(EffectDirection == "expo_strong_pos")%>%
  arrange(Proportion) %>%
  mutate(Anterior=factor(x=Anterior, levels=Anterior)) %>%
  select(Anterior) %>% unlist() %>% levels()

order.diseases_expo_f_pos <- order.diseases[(length(order.diseases)-9):length(order.diseases)]

order.diseases <- Expo %>%
  filter(EffectDirection == "expo_strong_neg")%>%
  arrange(-Proportion) %>%
  mutate(Anterior=factor(x=Anterior, levels=Anterior)) %>%
  select(Anterior) %>% unlist() %>% levels()

order.diseases_expo_f_neg <- order.diseases[(length(order.diseases)-9):length(order.diseases)]


## outcome 
out_strong_pos <- Results %>%
  mutate(across(everything(), ~ ifelse(. >= 0.01, 1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_weak_pos <- Results %>%
  mutate(across(everything(), ~ ifelse((. <= 0.01)&(. > 0), 1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_strong_neg <- Results %>%
  mutate(across(everything(), ~ ifelse(. <= -0.01, -1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_weak_neg <- Results %>%
  mutate(across(everything(), ~ ifelse((. >= -0.01)&(. < 0), -1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

Out <- cbind(out_strong_pos,out_weak_pos,out_strong_neg,out_weak_neg) %>%
  data.frame(.)%>%
  mutate("Posterior" = rownames(.))%>%
  pivot_longer(., cols = starts_with("out"), names_to = 'EffectDirection', values_to = 'Proportion')%>%
  mutate(EffectDirection = factor(EffectDirection, levels=c('out_strong_pos', 'out_weak_pos',
                                                            'out_strong_neg','out_weak_neg')))

order.diseases <- Out %>%
  filter(EffectDirection == "out_strong_pos")%>%
  arrange(Proportion) %>%
  mutate(Posterior=factor(x=Posterior, levels=Posterior)) %>%
  select(Posterior) %>% unlist() %>% levels()

order.diseases_out_f_pos <- order.diseases[(length(order.diseases)-9):length(order.diseases)]

order.diseases <- Out %>%
  filter(EffectDirection == "out_strong_neg")%>%
  arrange(-Proportion) %>%
  mutate(Posterior=factor(x=Posterior, levels=Posterior)) %>%
  select(Posterior) %>% unlist() %>% levels()

order.diseases_out_f_neg <- order.diseases[(length(order.diseases)-9):length(order.diseases)]

###############################################
## for males 

## Step 0 : Get results 
Results <- read.csv("./results_male/pairwise_correlation_matrix.csv",header=TRUE)[-1]

## expo
dim = dim(Results)[1]

expo_strong_pos <- Results %>%
  mutate(across(everything(), ~ ifelse(. >= 0.01, 1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_weak_pos <- Results %>%
  mutate(across(everything(), ~ ifelse((. <= 0.01)&(. > 0), 1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_strong_neg <- Results %>%
  mutate(across(everything(), ~ ifelse(. <= -0.01, -1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

expo_weak_neg <- Results %>%
  mutate(across(everything(), ~ ifelse((. >= -0.01)&(. < 0), -1, 0)))%>%
  rowSums(.,na.rm = TRUE)/dim*100

Expo <- cbind(expo_strong_pos,expo_weak_pos,expo_strong_neg,expo_weak_neg) %>%
  data.frame(.)%>%
  mutate("Anterior" = rownames(.))%>%
  pivot_longer(., cols = starts_with("exp"), names_to = 'EffectDirection', values_to = 'Proportion')%>%
  mutate(EffectDirection = factor(EffectDirection, levels=c('expo_strong_pos', 'expo_weak_pos',
                                                            'expo_strong_neg','expo_weak_neg')))


order.diseases <- Expo %>%
  filter(EffectDirection == "expo_strong_pos")%>%
  arrange(Proportion) %>%
  mutate(Anterior=factor(x=Anterior, levels=Anterior)) %>%
  select(Anterior) %>% unlist() %>% levels()

order.diseases_expo_m_pos <- order.diseases[(length(order.diseases)-9):length(order.diseases)]

order.diseases <- Expo %>%
  filter(EffectDirection == "expo_strong_neg")%>%
  arrange(-Proportion) %>%
  mutate(Anterior=factor(x=Anterior, levels=Anterior)) %>%
  select(Anterior) %>% unlist() %>% levels()

order.diseases_expo_m_neg <- order.diseases[(length(order.diseases)-9):length(order.diseases)]


## outcome 
out_strong_pos <- Results %>%
  mutate(across(everything(), ~ ifelse(. >= 0.01, 1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_weak_pos <- Results %>%
  mutate(across(everything(), ~ ifelse((. <= 0.01)&(. > 0), 1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_strong_neg <- Results %>%
  mutate(across(everything(), ~ ifelse(. <= -0.01, -1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

out_weak_neg <- Results %>%
  mutate(across(everything(), ~ ifelse((. >= -0.01)&(. < 0), -1, 0)))%>%
  colSums(.,na.rm = TRUE)/dim*100

Out <- cbind(out_strong_pos,out_weak_pos,out_strong_neg,out_weak_neg) %>%
  data.frame(.)%>%
  mutate("Posterior" = rownames(.))%>%
  pivot_longer(., cols = starts_with("out"), names_to = 'EffectDirection', values_to = 'Proportion')%>%
  mutate(EffectDirection = factor(EffectDirection, levels=c('out_strong_pos', 'out_weak_pos',
                                                            'out_strong_neg','out_weak_neg')))



order.diseases <- Out %>%
  filter(EffectDirection == "out_strong_pos")%>%
  arrange(Proportion) %>%
  mutate(Posterior=factor(x=Posterior, levels=Posterior)) %>%
  select(Posterior) %>% unlist() %>% levels()

order.diseases_out_m_pos <- order.diseases[(length(order.diseases)-9):length(order.diseases)]

order.diseases <- Out %>%
  filter(EffectDirection == "out_strong_neg")%>%
  arrange(-Proportion) %>%
  mutate(Posterior=factor(x=Posterior, levels=Posterior)) %>%
  select(Posterior) %>% unlist() %>% levels()

order.diseases_out_m_neg <- order.diseases[(length(order.diseases)-9):length(order.diseases)]


