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


barplot_exp_f_pos <-
  Expo %>%
  filter(Anterior %in% order.diseases_expo_f_pos)%>%
  mutate(Anterior=factor(Anterior, levels=order.diseases_expo_f_pos))%>%
  ggplot(aes(y = Anterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5, fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Anterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'expo_strong_pos' = const$fill_strong_pos, 'expo_weak_pos' = const$fill_weak_pos,
    'expo_strong_neg' = const$fill_strong_neg, 'expo_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = c(0.08,0.20)
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Pattens of effects along consequential spectrum rays \n Proportions of effect types',
    y = 'ICD10 codes (females)',
    fill = 'EffectDirection'
  )

plot(barplot_exp_f_pos)

barplot_exp_f_neg <-
  Expo %>%
  filter(Anterior %in% order.diseases_expo_f_neg)%>%
  mutate(Anterior=factor(Anterior, levels=order.diseases_expo_f_neg))%>%
  ggplot(aes(y = Anterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5, fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Anterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'expo_strong_pos' = const$fill_strong_pos, 'expo_weak_pos' = const$fill_weak_pos,
    'expo_strong_neg' = const$fill_strong_neg, 'expo_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Pattens of effects along consequential spectrum rays \nProportions of effeccts types',
    y = 'ICD10 codes (females)',
    fill = 'Effects'
  )

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

barplot_out_f_pos <-
  Out %>%
  filter(Posterior %in% order.diseases_out_f_pos)%>%
  mutate(Posterior=factor(Posterior, levels=order.diseases_out_f_pos))%>%
  ggplot(aes(y = Posterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5,  fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Posterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'out_strong_pos' = const$fill_strong_pos, 'out_weak_pos' = const$fill_weak_pos,
    'out_strong_neg' = const$fill_strong_neg, 'out_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Pattens of effects along causal spectrum rings \n Proportions of causal types',
    y = 'ICD10 codes (females)',
    fill = 'Effects'
  )

barplot_out_f_neg <-
  Out %>%
  filter(Posterior %in% order.diseases_out_f_neg)%>%
  mutate(Posterior=factor(Posterior, levels=order.diseases_out_f_neg))%>%
  ggplot(aes(y = Posterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5,  fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Posterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'out_strong_pos' = const$fill_strong_pos, 'out_weak_pos' = const$fill_weak_pos,
    'out_strong_neg' = const$fill_strong_neg, 'out_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Pattens of effects along causal spectrum rings \n Proportions of causal types ',
    y = 'ICD10 codes (females)',
    fill = 'Effects'
  )


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

barplot_exp_m_pos <-
  Expo %>%
  filter(Anterior %in% order.diseases_expo_m_pos)%>%
  mutate(Anterior=factor(Anterior, levels=order.diseases_expo_m_pos))%>%
  ggplot(aes(y = Anterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5,  fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Anterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'expo_strong_pos' = const$fill_strong_pos, 'expo_weak_pos' = const$fill_weak_pos,
    'expo_strong_neg' = const$fill_strong_neg, 'expo_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Proportions of effect types',
    y = 'ICD10 codes (males)',
    fill = 'Effects'
  )

barplot_exp_m_neg <-
  Expo %>%
  filter(Anterior %in% order.diseases_expo_m_neg)%>%
  mutate(Anterior=factor(Anterior, levels=order.diseases_expo_m_neg))%>%
  ggplot(aes(y = Anterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5,  fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Anterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'expo_strong_pos' = const$fill_strong_pos, 'expo_weak_pos' = const$fill_weak_pos,
    'expo_strong_neg' = const$fill_strong_neg, 'expo_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Proportions of effect types',
    y = 'ICD10 codes (males)',
    fill = 'Effects'
  )


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

barplot_out_m_pos <-
  Out %>%
  filter(Posterior %in% order.diseases_out_m_pos)%>%
  mutate(Posterior=factor(Posterior, levels=order.diseases_out_m_pos))%>%
  ggplot(aes(y = Posterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5,  fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Posterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'out_strong_pos' = const$fill_strong_pos, 'out_weak_pos' = const$fill_weak_pos,
    'out_strong_neg' = const$fill_strong_neg, 'out_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=14),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Proportions of causal types',
    y = 'ICD10 codes (males)',
    fill = 'Effects'
  )

barplot_out_m_neg <-
  Out %>%
  filter(Posterior %in% order.diseases_out_m_neg)%>%
  mutate(Posterior=factor(Posterior, levels=order.diseases_out_m_neg))%>%
  ggplot(aes(y = Posterior)) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 10.2, ymax = 10.5, fill = 'grey90') +
  geom_vline(
    xintercept = seq(-60, 60, 20),
    color = '#FFFFFF', size = 0.2
  ) +
  geom_col(
    aes(fill = EffectDirection, x = Proportion, group = Posterior),
    position = position_stack(), 
    data = .%>% data_frame(),
    width = 0.6
    #,show.legend =TRUE
  ) +
  geom_vline(xintercept = 0, color = 'grey50') +
  scale_fill_manual(values = c(
    'out_strong_pos' = const$fill_strong_pos, 'out_weak_pos' = const$fill_weak_pos,
    'out_strong_neg' = const$fill_strong_neg, 'out_weak_neg' = const$fill_weak_neg
  ),labels = c('Strong Increase','Weak Increase', 'Strong Reduction','Weak Reduction')
  ) +
  scale_x_continuous(breaks = seq(-100, 100, 20),
                     labels = c('100','80', '60', '40','20','0','20','40','60','80','100'), position = 'top') +
  scale_y_discrete(expand = expansion(mult = .01))+
  #scale_y_discrete(expand = expansion(add = 1.5))+
  theme_bw()+
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y=element_text(size=20),
    axis.title.x=element_text(size=20),
    axis.text=element_text(size=16),
    legend.text=element_text(size=20),
    #panel.background = element_rect(fill = 'grey95', color = 'grey95'),
    strip.background = element_rect(fill = 'grey90', color = 'grey90'),
    legend.position = "none"
  ) +
  guides(fill=guide_legend(title=" "))+
  labs(
    x = 'Proportions of causal types',
    y = 'ICD10 codes (males)',
    fill = 'Effects'
  )

plot(barplot_out_m_neg)


prow <- plot_grid(barplot_exp_f_pos, 
                  barplot_exp_f_neg,
                  barplot_exp_m_pos,
                  barplot_exp_m_neg,
                  barplot_out_f_pos, 
                  barplot_out_f_neg, 
                  barplot_out_m_pos,
                  barplot_out_m_neg,
                  rel_widths = 1, 
                  nrow = 4,labels = c("a", "b" ,"","","c", "d","",""),label_size = 30)
prow

## Patterns of effects along the consequential spectrum rays and causal spectrum rings
pdf("results/figures/fig3_effects_rays_pos_neg_10ICD.pdf", width  = 36, height = 18)
plot(prow)
dev.off()

