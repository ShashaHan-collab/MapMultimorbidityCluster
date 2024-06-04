# This is an example for a female disease.
# The code for male diseases should be modified with the corresponding data.
library(dplyr)
library(glmnet)
library(doParallel)
setwd("~/")
glmnet.control(epsnr = 1e-5, mxitnr = 25)

lr_tmle <- function(dcode2, diffAY = 0){
  diseaseAdate <- paste0(dcode1, "date")
  diseaseYdate <- paste0(dcode2, "date")
  # select the subpopulation for studying the effects of diseaseA on disease Y
  ukbdata_pairwise <- ukbdata_i %>% 
    mutate(diff_datesAY = as.numeric(difftime(get(diseaseYdate), get(diseaseAdate), units = "days")),
           diff_datesAI = as.numeric(difftime(get(diseaseAdate), ini_date, units = "days")),
           diff_datesYI = as.numeric(difftime(get(diseaseYdate), ini_date, units = "days")))%>%
    # exclude the people who had diseaseA prior to the baseline 
    filter(is.na(diff_datesAI) | diff_datesAI>=0 )%>%
    # exclude the people who had diseaseY prior to the baseline 
    filter(is.na(diff_datesYI) | diff_datesYI>=0 )%>%
    # include the people who did not have diseaseY prior to diseaseA after baseline
    filter(is.na(diff_datesAY) | diff_datesAY>=diffAY )%>%
    select(all_of(c(cols.x,dcode1,dcode2)), diff_datesAY)
  
  ## Step 0. Prefit the propensity score model and exclude individuals with extreme values
  A <- ukbdata_pairwise[,dcode1] %>% as.numeric()
  X <- ukbdata_pairwise %>% select(all_of(cols.x))
  X_dummy = model.matrix(~., X)[,-1]
  X_dummy = as(X_dummy, "dgCMatrix")
  g <- glmnet(X_dummy, A, family = binomial(), lambda = 0, lower.limits = -8, upper.limits = 8, maxit = 200, thresh = 1e-6)
  if (g$jerr<0) return(c(dcode1, dcode2, rep(NA, 9)))
  g_w = as.vector(predict(g, s = 0, newx = X_dummy, type = 'response'))
  g_min = min(g_w[A == 1])
  g_max = max(g_w[A == 0])
  ukbdata_pairwise = ukbdata_pairwise %>% filter(g_w >= g_min & g_w <= g_max)
  rm(X_dummy, X); gc()
  
  minyears <- round(min(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  maxyears <- round(max(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  medyears <- round(median(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  sp_num <- nrow(ukbdata_pairwise)
  
  if(sp_num <= 3){
    cat("Number of sample size is smaller than 3 in the scenario of",dcode1,"to",dcode2)
    return(c(dcode1, dcode2, sp_num, rep(NA, 8)))
  }
  
  ## Step 1. Learn the propensity score model
  Y <- ukbdata_pairwise[,dcode2]%>% as.numeric()
  A <- ukbdata_pairwise[,dcode1]%>% as.numeric()
  X <- ukbdata_pairwise %>% select(all_of(cols.x))
  rm(ukbdata_pairwise); gc()
  
  X_dummy = model.matrix(~., X)[,-1]
  X_dummy = as(X_dummy, "dgCMatrix")
  g <- glmnet(X_dummy, A, family = binomial(), lambda = 0, lower.limits = -8, upper.limits = 8, maxit = 200, thresh = 1e-6)
  if (g$jerr<0) return(c(dcode1, dcode2, sp_num, rep(NA, 8)))
  g_w = as.vector(predict(g, s = 0, newx = X_dummy, type = 'response'))
  H_1 <- 1/g_w
  H_0 <- -1/(1-g_w) # Pr(A=0|X) is 1-Pr(A=1|X)
  H_A <- ifelse(A == 1, H_1, H_0) # if A is 1 (treated), assign H_1, else H_0
  rm(X_dummy); gc()
  
  ## Step 2. Learn the outcome model
  X_A_dummy = model.matrix(~., cbind(A, X))[,-1]
  X_A_dummy = as(X_A_dummy, "dgCMatrix")
  Q <- glmnet(X_A_dummy, Y, family = binomial(), lambda = 0, lower.limits = -8, upper.limits = 8, thresh = 1e-6)
  if (Q$jerr<0) return(c(dcode1, dcode2, sp_num, rep(NA, 8)))
  coefA = coef(Q, s = 0)['A',1]
  Q_A_link = as.vector(predict(Q, s = 0, newx = X_A_dummy))
  Q_A = pmax(pmin(plogis(Q_A_link), 1 - 1e-6), 1e-6)
  Q_1 = pmax(pmin(plogis(Q_A_link + (1-A)*coefA), 1 - 1e-6), 1e-6)
  Q_0 = pmax(pmin(plogis(Q_A_link - A*coefA), 1 - 1e-6), 1e-6)
  rm(X_A_dummy); gc()
  
  ## Step 3. Estimate the flunctuation parameter model
  glm_fit <- glm(Y ~ -1 + offset(qlogis(Q_A)) + H_A, family=binomial)
  eps <- coef(glm_fit)
  minmax = max(abs(glm_fit$linear.predictors))
  extreme_value = ifelse(minmax>30, 1, 0)
  ## Step 4. Estimate the initial parameters
  Q_A_update <- plogis(qlogis(Q_A) + eps*H_A)
  Q_1_update <- plogis(qlogis(Q_1) + eps*H_1)
  Q_0_update <- plogis(qlogis(Q_0) + eps*H_0)
  
  ## Step 5. Compute the statistics parameters
  tmle_ate <- mean(Q_1_update - Q_0_update)
  
  ## Step 6. Calculate the Standard Errors for Confidence Intervals and P-values
  infl_fn <- (Y - Q_A_update) * H_A + Q_1_update - Q_0_update - tmle_ate
  tmle_se <- sqrt(var(infl_fn)/sp_num)
  conf_low <- tmle_ate - 1.96*tmle_se
  conf_high <- tmle_ate + 1.96*tmle_se
  pval <- 2 * (1 - pnorm(abs(tmle_ate / tmle_se)))
  
  estimates <- round(c(tmle_ate, conf_low, conf_high,pval),5)
  update.estimates <- c(dcode1, dcode2, sp_num, estimates,
                        minyears, medyears, maxyears, extreme_value)
  return(update.estimates)
}

############################################
## set disease states at baseline 
icd10file <- readxl::read_excel('disease_prevalence_sex_v4.xlsx')
selectlabel <- icd10file %>% 
  filter(!stringr::str_detect(ICD10, '[Q-T]|[V-Y]\\d{2}')) %>% 
  filter(prevalence_female >= 0.01) %>% pull(ICD10)

############################################
cols.var <- c("age", "tdi","edu","ethnic",
              "BMI","IPAQ", "smoke","drink","sleep","phone_use",
              "DBP", "SBP","HbA1c", "HDL", "LDL", "TG","calc","CRP" , 
              "gas", "leisure",
              "agedeathmoth", "agedeathfath","longill","illfather","illmother")
cols.disease.base <- paste0(selectlabel, "_base")
cols.x <- c(cols.var, cols.disease.base)

selectlabel <- selectlabel[!stringr::str_detect(selectlabel, 'Z\\d{2}')] # remove ZXX 
dnum <- length(selectlabel) # number of diseases under consideration

## do for each of the five imputed datasets
imp_num = 5
path = 'main_analysis/'
for (i in 1:imp_num){
  ukbdata_i = readRDS(paste0('ukbdata_analysis_comp/ukbdata_imputed_female_', i, '.rds'))
  
  ## format output
  CORRmat <- data.frame(matrix(NA,dnum,dnum+1))
  TMLEmat <- data.frame(matrix(NA,0, 11))
  names(CORRmat) <- c("expo", selectlabel)
  CORRmat$expo <- selectlabel
  names(TMLEmat) <- c("expo","outcome","sample_size","mean","low95ci","up95ci","pvalue","mininterval","medinterval","maxinterval","extreme_value")
  
  ## estimation
  cl <- makeCluster(getOption("cl.cores", 10)) # memory limit
  registerDoParallel(cl)
  for (dcode1 in selectlabel){
    cat(sprintf("imp:%d treat:(%d / %d)", i, match(dcode1, selectlabel), dnum), as.character(Sys.time()), '\n')
    selectlabel2 <- selectlabel[!selectlabel %in% c(dcode1)]
    tmle_rlt <- foreach(dcode2=selectlabel2, .packages=c('glmnet','dplyr'),
                       .combine = rbind) %dopar% lr_tmle(dcode2, diffAY = 0)
    colnames(tmle_rlt) <- c("expo","outcome","sample_size","mean","low95ci","up95ci","pvalue","mininterval","medinterval","maxinterval","extreme_value")
    TMLEmat <- rbind(TMLEmat, tmle_rlt, make.row.names = FALSE)
    CORRmat[CORRmat$expo == dcode1, selectlabel2] <- tmle_rlt[,'mean']
    
    write.csv(TMLEmat, file = paste0(path,'results_female/pairwise_estimates_imp',i,'.csv'))
    write.csv(CORRmat, file = paste0(path,'results_female/pairwise_correlation_matrix_imp',i,'.csv'))
    gc()
  }
  stopCluster(cl)
}
