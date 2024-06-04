library(dplyr)
library(data.table)
library(stringr)
library(parallel)
library(readxl)
library(openxlsx)
library(mice)

setwd("~/")

## auxliary functions
isTRUE <- function(x){
  !is.na(x) & x
}

## extract ICD-10 disease and diagnosis dates
alldatac = readRDS('ukb_variables.rds')
ICD10 = alldatac %>% select(contains('41270'))
ICD10_date = alldatac %>% select(contains('41280'))
extract_date = alldatac %>% select(eid, sex)
rm(alldatac); gc()

# get disease frequencies stratified by sex
n_male = sum(extract_date$sex)
n_female = sum(1 - extract_date$sex)
regx1 = '^[A-Z]\\d{2}'
ICD10_2 = ICD10 %>% mutate(across(, ~str_extract(.x, regx1)))
rm(ICD10); gc()
alldata = cbind(extract_date, ICD10_2)
code_male = alldata %>% filter(sex == 1) %>% select(-eid, -sex) %>%
  unlist() %>% table() %>% data.frame() %>% arrange(desc(Freq))
code_female = alldata %>% filter(sex == 0) %>% select(-eid, -sex) %>%
  unlist() %>% table() %>% data.frame() %>% arrange(desc(Freq))
write.csv(code_male, 'ICD10_code_freq_male.csv', row.names = F)
write.csv(code_female, 'ICD10_code_freq_female.csv', row.names = F)
rm(alldata); gc()

# extract diseases with prevalence < 1%
code_male = code_male[code_male$Freq >= n_male*0.01, 1] %>% as.character()
code_female = code_female[code_female$Freq >= n_female*0.01, 1] %>% as.character()
code2_list = union(code_male, code_female)
cl <- makeCluster(getOption("cl.cores", 6))
extract_date = data.frame(extract_date)
for (k in 1:length(code2_list2)){
  codei = code2_list2[k]
  d10_date = ICD10_date
  d10_date[isTRUE(ICD10_2!=codei)] = NA
  d10_date = parApply(cl, d10_date, 1, min, na.rm = TRUE)
  d10_date = as.IDate(d10_date)
  extract_date[codei] = d10_date
  rm(d10_date); gc()
  cat(k,'/ ', length(code2_list2), as.character(Sys.time()), '\n')
  if (k %% 10 == 0){
    saveRDS(extract_date, file = 'ukb_icd10_disease_sex.rds')
  }
}
saveRDS(extract_date, file = 'ukb_icd10_disease_sex.rds')
stopCluster(cl)

# calculate disease prevalence
ukb_icd10 = readRDS('ukb_icd10_disease_sex.rds')
ukb_cov = readRDS('ukbdata_v2.rds')
ukb_icd10 = ukb_cov %>% select(eid, sex) %>% left_join(ukb_icd10)
icd10 = data.frame(ICD10 = names(ukb_icd10), samples = colSums(!is.na(ukb_icd10)))
icd10 = icd10[3:nrow(icd10),]
row.names(icd10) = NULL
icd10$prevalence = icd10$samples/nrow(ukb_icd10)
icd10$sample_male = colSums(!is.na(ukb_icd10 %>% filter(sex == 1) %>% select(-eid,-sex)))
icd10$sample_female = colSums(!is.na(ukb_icd10 %>% filter(sex == 0) %>% select(-eid,-sex)))
icd10$prevalence_male = icd10$sample_male/sum(ukb_icd10$sex)
icd10$prevalence_female = icd10$sample_female/sum(1-ukb_icd10$sex)
icd10 = icd10 %>% filter(prevalence_male>=0.01 | prevalence_female>=0.01) %>% arrange(desc(prevalence))
write.xlsx(icd10, file = 'disease_prevalence_sex_v4.xlsx')

## merge ICD-10 diagnosis with other variables
disease_sel = icd10$ICD10
ukbdata_v3 = ukb_cov %>% select(eid:illmother)
ukb_icd10 = ukb_icd10 %>% select(eid, all_of(disease_sel))
ukb_icd10_binary = ukb_icd10 %>% mutate(across(, ~as.numeric(!is.na(.x)))) %>% select(-eid)
names(ukb_icd10) = c('eid', paste0(disease_sel, 'date'))
ukb_icd10 = cbind(ukb_icd10, ukb_icd10_binary)
ukbdata_v3 = ukbdata_v3 %>% left_join(ukb_icd10)
saveRDS(ukbdata_v3, file = 'ukbdata_v3.rds')

## multiple imputation using mice
ukbdata = readRDS("ukbdata_v3.rds")
ukbdata = ukbdata %>% mutate(death_date = unlist(death_date))

# calculate missing percentage
num_list = c( "age", "tdi", "BMI", "sleep", "phone_use",
               "DBP", "SBP","HbA1c", "HDL", "LDL", "TG", "calc", "CRP" , 
               "agedeathmoth","agedeathfath")
bin_list = c("ethnic", "smoke", "leisure",
              "longill","illfather","illmother")
poly_list = c("edu","IPAQ", "drink", "gas")
map(ukbdata[,c(num_list, bin_list,poly_list)], ~mean(is.na(.)))
map(ukbdata[,c('ini_date',"death_date")], ~mean(is.na(.)))

# impute missing data
for (sex0 in 0:1){
  sex_name = ifelse(sex0==0, 'female', 'male')
  ukbdata_subset = ukbdata %>% filter(sex == sex0)
  eid = ukbdata_subset$eid
  # we do not consider disease at baseline due to computational issues.
  ukbdata_to_impute = ukbdata_subset %>% select(c(num_list, bin_list,poly_list))
  ukbdata_to_impute = ukbdata_to_impute %>% 
    mutate_at(num_list, as.numeric)%>%
    mutate_at(c(bin_list,poly_list), as.factor)

  init = mice(ukbdata_to_impute, maxit=0) 
  meth = init$method
  meth[num_list] = "pmm" 
  meth[bin_list] = "logreg" 
  meth[poly_list] = "polyreg"
  
  set.seed(123)
  ukbdata_imputed = mice(ukbdata_to_impute, method = meth, m = 5, maxit = 3)
  ukbdata_imputed_comp = complete(ukbdata_imputed, action = "long")
  ukbdata_imputed_comp = ukbdata_imputed_comp %>% mutate(eid = rep(eid,5)) %>% select(-.id)
  saveRDS(ukbdata_imputed_comp, file = paste0('ukbdata_analysis_comp/ukbdata_imputed_',sex_name,'.rds'))
}
rm(list = ls());gc()

# combine imputed datasets with remaining datasets
ukbdata = readRDS("ukbdata_v3.rds")
ukbdata = ukbdata %>% mutate(death_date = unlist(death_date))

cols.var = c("sex", "age", "tdi","edu","ethnic",
              "BMI","IPAQ", "smoke","drink","sleep","phone_use",
              "DBP", "SBP","HbA1c", "TC", "HDL", "LDL", "TG","calc","chol","CRP" , 
              "gas", "leisure",
              "agedeathmoth", "agedeathfath","longill","illfather","illmother")
ukbdata = ukbdata %>% select(-cols.var)
for (sex0 in 0:1){
  sex_name = ifelse(sex0==0, 'female', 'male')
  ukbvar_comp = readRDS(paste0('ukbdata_analysis_comp/ukbdata_imputed_',sex_name,'.rds'))
  ukbdata_subset = ukbvar_comp %>% left_join(.,ukbdata,by = c("eid"))
  
  ## set disease states at baseline 
  icd10file = read_excel('disease_prevalence_sex_v4.xlsx')
  selectlabel = icd10file$ICD10
  for(dcode in selectlabel){
    diseaseDdate = paste0(dcode, "date")
    ukbdata_subset = ukbdata_subset %>% 
      mutate_at(c(eval(diseaseDdate),"ini_date"), as.Date, format = "%Y-%m-%d")%>%
      mutate(diff_datesDI = as.numeric(difftime(get(diseaseDdate), ini_date, units = "days")))%>%
      mutate(!!paste0(dcode,'_base') := ifelse((is.na(diff_datesDI) | diff_datesDI>=0 ),0,1))
  }    
  ukbdata_subset = ukbdata_subset %>% select(-diff_datesDI)
  labels.base = paste0(selectlabel, "_base")
  
  # save the data set for each imputation
  imp_num = 5
  for (i in 1:imp_num){
    ukbdata_subset_i = ukbdata_subset %>% filter(.imp == i)
    saveRDS(ukbdata_subset_i, file = paste0('ukbdata_analysis_comp/ukbdata_imputed_',sex_name,'_', i, '.rds'))
  }
  gc()
}







