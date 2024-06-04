library(dplyr)
library(glmnet)
library(doParallel)

glmnet.control(epsnr = 1e-5, mxitnr = 25)
detectCores(logical = F)
lr_tmle <- function(dcode2){
  diseaseAdate <- paste0(dcode1, "date")
  diseaseYdate <- paste0(dcode2, "date")
  # select the subpopulation for studying the effects of diseaseA on disease Y
  ukbdata_pairwise <- ukbdata_i %>% 
    # mutate_at(c(eval(diseaseAdate),eval(diseaseYdate),"ini_date"), as.Date, format = "%Y-%m-%d")%>%
    mutate(diff_datesAY = as.numeric(difftime(get(diseaseYdate), get(diseaseAdate), units = "days")),
           diff_datesAI = as.numeric(difftime(get(diseaseAdate), ini_date, units = "days")),
           diff_datesYI = as.numeric(difftime(get(diseaseYdate), ini_date, units = "days")))%>%
    # exclude the people who had diseaseA prior to the baseline 
    filter(is.na(diff_datesAI) | diff_datesAI>=0 )%>%
    # exclude the people who had diseaseY prior to the baseline 
    filter(is.na(diff_datesYI) | diff_datesYI>=0 )%>%
    # include the people who did not have diseaseY prior to diseaseA after baseline
    filter(is.na(diff_datesAY) | diff_datesAY>=0 )%>%
    select(all_of(c(cols.x,dcode1,dcode2)), diff_datesAY)
  
  minyears <- round(min(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  maxyears <- round(max(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  medyears <- round(median(ukbdata_pairwise$diff_datesAY,na.rm = TRUE)/365.2425,3)
  sp_num <- nrow(ukbdata_pairwise)
  
  if(sp_num <= 3){
    cat("Number of sample size is samller than 3 in the scenario of",dcode1,"to",dcode2)
    continue
  }
  
  ## Step 1. Learn the outcome model
  col.treat_out <- c(dcode1,dcode2)
  Y <- ukbdata_pairwise[,col.treat_out[2]]%>% as.numeric()
  A <- ukbdata_pairwise[,col.treat_out[1]]%>% as.numeric()
  X <- ukbdata_pairwise %>% select(all_of(cols.x)) %>% select(-all_of(paste0(col.treat_out, '_base')))
  rm(ukbdata_pairwise); gc()

  
  X_A_dummy = model.matrix(~., cbind(A, X))[,-1]
  Q <- glmnet(X_A_dummy, Y, family = binomial(), lambda = 0, maxit = 50, thresh = 1e-6)
  if (Q$jerr<0) return(rep(NA, 9))
  coefA = coef(Q, s = 0)['A',1]
  Q_A_link = as.vector(predict(Q, s = 0, newx = X_A_dummy))
  Q_A = plogis(Q_A_link)
  Q_1 = plogis(Q_A_link + (1-A)*coefA)
  Q_0 = plogis(Q_A_link - A*coefA)
  rm(X_A_dummy); gc()
  
  ## Step 2. Learn the propensity score model
  X_dummy = model.matrix(~., X)[,-1]
  g <- glmnet(X_dummy, A, family = binomial(), lambda = 0, maxit = 50, thresh = 1e-6)
  if (g$jerr<0) return(rep(NA, 9))
  g_w = as.vector(predict(g, s = 0, newx = X_dummy, type = 'response'))
  H_1 <- 1/g_w
  H_0 <- -1/(1-g_w) # Pr(A=0|X) is 1-Pr(A=1|X)
  H_A <- ifelse(A == 1, H_1, H_0) # if A is 1 (treated), assign H_1, else H_0
  rm(X_dummy); gc()
  
  ## Step 3. Estimate the flunctuation parameter model
  glm_fit <- glm(Y ~ -1 + offset(qlogis(Q_A)) + H_A, family=binomial)
  eps <- coef(glm_fit)
  
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
  
  estimates <- round(c(tmle_ate, conf_low, conf_high,pval),3)
  update.estimates <- c(dcode1, dcode2,estimates,minyears,medyears,maxyears)
  return(update.estimates)
}

############################################
## set disease states at baseline 
icd10file <- readxl::read_excel('disease_prevalence_v2.xlsx')
selectlabel <- icd10file %>% 
  filter(!stringr::str_detect(ICD10, '[O-P]|[R-T]|[V-Z]\\d{2}')) %>% 
  filter(prevalence >= 0.01) %>% pull(ICD10)

############################################
## format output
dnum = length(selectlabel); #number of diseases under consideration
CORRmat <- data.frame(matrix(NA,dnum,dnum+1))
TMLEmat <- data.frame(matrix(NA,0, 9))

names(CORRmat) <- c("expo",selectlabel)
CORRmat$expo <- selectlabel
names(TMLEmat) <- c("expo","outcome","mean","low95ci","up95ci","pvalue","mininterval","medinterval","maxinterval")

############################################
## estimation
cols.var <- c("sex", "age", "tdi","edu","ethnic",
              "BMI","IPAQ", "smoke","drink","sleep","phone_use",
              "DBP", "SBP","HbA1c", "HDL", "LDL", "TG","calc","CRP" , 
              "gas",
              "leisure",
              "agedeathmoth", "agedeathfath","longill","illfather","illmother")

cols.disease.base <- paste0(selectlabel, "_base")
cols.x <- c(cols.var,cols.disease.base)

## do for each of the five imputed datasets
imp_num = 5

for (i in 1:imp_num){
  ## Uncomment below for the final analysis in the server
  ukbdata_i = readRDS(paste0('ukbdata_analysis_comp/ukbdata_imputed_', i, '.rds'))
  
  cl <- makeCluster(getOption("cl.cores", 4))
  registerDoParallel(cl)
  
  for (dcode1 in selectlabel){
    cat(sprintf("imp:%d treat:(%d / %d)", i, match(dcode1, selectlabel), dnum), as.character(Sys.time()), '\n')
    
    selectlabel2 <-  selectlabel[! selectlabel %in% c(dcode1)]
    tmle_rlt = foreach(dcode2=selectlabel2, .packages=c('glmnet','dplyr'),
                       .combine = rbind) %dopar% lr_tmle(dcode2)
    colnames(tmle_rlt) = c("expo","outcome","mean","low95ci","up95ci","pvalue","mininterval","medinterval","maxinterval")
    TMLEmat = rbind(TMLEmat, tmle_rlt, make.row.names = FALSE)
    CORRmat[CORRmat$expo == dcode1, selectlabel2] <- tmle_rlt[,'mean']
    
    write.csv(TMLEmat, file = paste0('results/pairwise_estimates_imp',i,'.csv'))
    write.csv(CORRmat, file = paste0('results/pairwise_correlation_matrix_imp',i,'.csv'))
  }
  stopCluster(cl)
}
