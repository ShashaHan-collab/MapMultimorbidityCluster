
library(dplyr)
library(data.table)
library(stringr)
library(openxlsx)
library("readxl")

setwd("~/")

## auxliary functions
std_date <- function(x){
  if (is.na(x)) return(NA)
  if (x<1000) return(NA)
  year = str_extract(as.character(x), '^\\d{4}')
  if (is.na(year)) {cat(x,' '); return(NA)}
  mon = str_extract(as.character(x), '(?<=\\d{4}\\.)\\d{1,2}')
  if (is.na(mon)) {cat(x,' '); return(paste0(year,'-01-01'))}
  if (nchar(mon)==1) mon = paste0('0', mon)
  return(paste0(year,'-',mon,'-','01'))
}
replaceNaN <- function(x){
  x[is.nan(x)] = NA
  x
}

## merge separate data files
alldata = fread('./data/ukb_data/field_select.csv')
addf2 = fread('./data/ukb_data/field_select2.csv')
addf3 = fread('./data/ukb_data/field_select3.csv')
addf4 = fread('./data/ukb_data/field_select4.csv')
addf5 = fread('./data/ukb_data/field_select5.csv')
addf6 = fread('./data/ukb_data/field_select6.csv')

alldatac = alldata %>%# left_join(eyeimg) %>% 
  left_join(addf2) %>%  left_join(addf3) %>%  left_join(addf4) %>%  
  left_join(addf5) %>% left_join(addf6) 

rm(alldata); gc()
## 40000: date of death;
## 40001: primary cause of death;
## 40002: secondary cause of death.




rm(list = ls()); gc()

## clean data
alldatac = read.csv('data/data_clean_v1.csv')
colname = names(alldatac)

varset = c('eid', 'X31.0.0', 'X34.0.0', 'X52.0.0', 
           'X53.0.0', #'X53.1.0',#'X53.2.0','X53.3.0',
           'X4079.0.0', 'X4079.0.1', #'X4079.1.0', 'X4079.1.1', 'X4079.2.0', 'X4079.2.1', 'X4079.3.0', 'X4079.3.1', 
           'X4080.0.0', 'X4080.0.1', 
           'X21001.0.0',
           'X22032.0.0')
ukbdata_cleaned = alldatac %>% select(all_of(varset)) 
names(ukbdata_cleaned) = c('eid', 'sex', 'birth_y', 'birth_m', 
                'ini_date', #'Rep_date',
                'DBP1', 'DBP2', 
                'SBP1', 'SBP2', 
                'BMI',
                'IPAQ')

## 31:Sex
## 34:Year of birth
## 52:Month of birth
## 53:Date of attending assessment centre
## 4079:Diastolic blood pressure, automated reading
## 4080:Systolic blood pressure, automated reading
## 21001: BMI
## 22032: derived MET (Metabolic Equivalent Task) scores data based on IPAQ (International Physical Activity Questionnaire) guidelines
#age(y) at Ini_date,IPAQ: https://biobank.ndph.ox.ac.uk/showcase/ukb/docs/ipaq_analysis.pdf
ukbdata_cleaned = ukbdata_cleaned %>% mutate(birth_d = paste0(birth_y, '.', birth_m,'.1'), .after = birth_m) %>% 
  mutate(birth_d = as.Date(birth_d, format = '%Y.%m.%d')) %>% 
  mutate(age = (as.Date(ini_date) - birth_d)/365.25, .after = sex) %>% 
  mutate(age = as.numeric(age))
mean(is.na(ukbdata_cleaned$age))
mean(is.na(ukbdata_cleaned$sex))
#SBP/DBP at Ini_date
ukbdata_cleaned$DBP = rowMeans(ukbdata_cleaned[,c('DBP1','DBP2')], na.rm = T) %>% 
  tidyr::replace_na()
ukbdata_cleaned$SBP = rowMeans(ukbdata_cleaned[,c('SBP1','SBP2')], na.rm = T) %>% 
  tidyr::replace_na()

#HbA1c
HbA1c = alldatac %>% select(contains('30750'))
ukbdata_cleaned$HbA1c = HbA1c[,1]
#TC/HDL/LDL/TG
TC = alldatac %>% select(contains('30690'))
ukbdata_cleaned$TC = TC[,1]
HDL = alldatac %>% select(contains('30760'))
ukbdata_cleaned$HDL = HDL[,1]
LDL = alldatac %>% select(contains('30780'))
ukbdata_cleaned$LDL = LDL[,1]
TG = alldatac %>% select(contains('30870'))
ukbdata_cleaned$TG = TG[,1]
mean(is.na(HbA1c[,1]))
mean(is.na(TC[,1]))
mean(is.na(HDL[,1]))
mean(is.na(LDL[,1]))
mean(is.na(TG[,1]))

ukbdata_cleaned <- ukbdata_cleaned %>% 
  select(-birth_y, -birth_m, -DBP1, -DBP2, -SBP1, -SBP2) 

#edu -3:NA, -7:0, 2-6:1, 
#1. college or university:2 
#2.to 6. From A levels to other professional qualifications
#-7: none of qualifications
edu = alldatac %>% select(contains('6138.0'))
eduf = c()
for (i in 1:nrow(edu)){
  edui = edu[i, ]
  edui = edui[!is.na(edui)]
  if (length(edui)==0 | edui[1] == -3) {eduf[i] = NA; next}
  edui = sort(edui)
  if (edui[1] == -7) eduf[i] = 0
  if (edui[1] == 1) eduf[i] = 2
  if (edui[1] > 1) eduf[i]=1
}
table(eduf, useNA = 'ifany')
ukbdata_cleaned$edu = eduf

#Ethnic background -3/-1:NA  white(1):0 non-white(>1):1
ethnic = alldatac %>% select(contains('21000.0'))
ethnic = ethnic[,1]
ethnic[which(ethnic == -3 | ethnic == -1)] = NA
ethnic0 = str_extract(ethnic, '^\\d') %>% as.numeric()
table(ethnic0, useNA = 'ifany')
ethnic0[ethnic0 == 1] = 0
ethnic0[ethnic0 > 1] = 1
ukbdata_cleaned = ukbdata_cleaned %>% mutate(ethnic = ethnic0)

#smoke  -3:NA 0:0 >0:1
smoke = alldatac %>% select(contains('1239.0'))
smoke = smoke[,1]; smoke[which(smoke == -3)] = NA
table(smoke, useNA = 'ifany')
smoke[smoke>0] = 1
ukbdata_cleaned = ukbdata_cleaned %>% mutate(smoke = smoke)

#drink -3:NA 6:0 2-5:1, 1:2
drink = alldatac %>% select(contains('1558.0'))
drink = drink[,1]; drink[which(drink == -3)] = NA
table(drink, useNA = 'ifany')
drink = factor(drink, levels = 1:6, labels = c(2,1,1,1,1,0))
ukbdata_cleaned = ukbdata_cleaned %>% mutate(drink = drink)

#sleep -3/-1:NA
sleep = alldatac %>% select(contains('1160.0'))
sleep = sleep[,1]; sleep[which(sleep == -3 | sleep == -1)] = NA
table(sleep, useNA = 'ifany')
ukbdata_cleaned = ukbdata_cleaned %>% mutate(sleep = sleep)

#phone_use -3/-1:NA
phone_use = alldatac %>% select(contains('1120.0'))
phone_use = phone_use[,1]; phone_use[which(phone_use == -3 | phone_use == -1)] = NA
table(phone_use, useNA = 'ifany')
ukbdata_cleaned = ukbdata_cleaned %>% mutate(phone_use = phone_use)

# death date
death_date = alldatac %>% select(contains('40000.0.0'))
hist(as.Date(death_date$X40000.0.0), breaks = 'year')
ukbdata_cleaned <- ukbdata_cleaned %>% mutate(death_date = death_date$X40000.0.0)%>% unlist()

death_cause = alldatac %>% select(contains(c('40001.0.0','40002.0')))
death_cause[death_cause==""]<-NA
not_all_na <- function(x) any(!is.na(x))
not_any_na <- function(x) all(!is.na(x))
death_cause <- death_cause %>% sapply(., substring, 1, 3) %>% unlist()%>%
          data.frame %>%
          select(where(not_all_na))
  
names(death_cause) <- paste0("deathcause",1:15)

death.date = alldatac %>% select(contains('40000.0.0'))%>% unlist()
ukbdeath = alldatac %>% select("eid") %>%
  mutate(death_date = death.date)

ukbdeath = cbind(ukbdeath, death_cause)
saveRDS(ukbdeath, file = "data/ukbdeath_v1.rds")


## extract disease and diagnosis dates
ukbdisease = alldatac %>% select("eid") 
asc <- function(x){strtoi(charToRaw(x),16L)}
## 65 for the asc II code. 
## 2600 for the maximal number of ICD10. 
code2num <- function(x){
  x1 = substring(x,first=1L,last=1L)
  x2 = substring(x,first=2L,last=10000L)
  return((asc(x1)-65)*100+as.numeric(x2))
}
num2code <- function(x){
  a = floor(x/100)
  x1 = rawToChar(as.raw(a+65))
  x2 = x-100*a
  if(x2 >= 10){
    x2 = as.character(x-100*a)
    return(paste0(x1,x2))
  }
  else{
    x2 = as.character(x-100*a)
    return(paste0(x1,'0',x2))
  }
}
ICD10 <- alldatac %>% select(contains('41270'))
ICD10_date <- alldatac %>% select(contains('41280'))
N = 2600 # the total number of disease types coded in ICD10
selectlabel <- c("I10","Z86","Z92","E78","K57","Z87","Z88","K29","R10","R07",
                 "K21","Z85","K44","Z96","I25","Z53","J45","M19","Z51","Z90",
                 "E11","H26","K62","M17","N39","I48","E66","Y83","K63","I20")


icd10file <- read_excel('data/disease_prevalence.xlsx')
selectlabel100 <- icd10file$ICD10
selectlabel2 <- selectlabel100[!selectlabel100 %in% selectlabel]

for(dcode in selectlabel2){
  d10 = apply(ICD10, 1, function(x) which(str_detect(x, pattern = dcode)))
  d10_date = c()
  for (n in 1:nrow(ICD10)){
    d10_date[n] = ifelse(length(d10[[n]]) == 0, NA, ICD10_date[n, d10[[n]]] %>% unlist() %>% min())
  }
  d10[sapply(d10, length) > 0 ] <- 1
  d10[sapply(d10, length) ==0 ] <- 0
  temp <- as.data.frame(cbind(unlist(d10),d10_date))
  names(temp) <- c(dcode, paste0(dcode, "date"))
  ukbdisease <- cbind(ukbdisease,temp)
}

ukbdata <- merge(ukbdata_cleaned,ukbdisease,by="eid")

saveRDS(ukbdata, file = "data/ukbdata_v1.rds")



