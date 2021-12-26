# library(tibble)
library(readr)
###
##
#
# File to analyse and study the data to derive the final data
# to predict (diagnose) AD among people with cognitive impairment
# set wd and read data
wd <- getwd()
setwd("..")

raw_data <- read_csv("investigator_nacc55.csv") # compare to read.csv - maybe factors are being coherced to NA
setwd(wd)
nobs <- nrow(raw_data)

# select columns
base <- c("NACCID","BIRTHYR","SEX","NACCAGE","RACE","EDUC","NACCAPOE","NACCNE4S")
# "HISPANIC"
# NACCAPOE is the APOE mutation, NACCNE4S number of e4 allelles.
# NACCAGE - age at visit

rm_apoe <- which(raw_data[,"NACCNE4S"]==9) # 9= not known not assessed

#visit_t <- c("VISITYR","NACCAGE")

dxs <- c("DEMENTED","COGFLAGO","NACCCOGF","COGMODE","NACCALZD",
         "FDGAD","TAUPETAD","CSFTAU","OTHBIOM")
# NACCALZD - Alzheimer as etiologic cause of cognitive disorder;
# ""NACCALZP" - Primary cause of cognitive disorder
# "DXMETHOD"; COGFLAGO = age begin symptoms; 
# NACCCOGF - predominant symptom that was first recognized as decline in cognition
# COGMODE - mode of onset of cognitive symptoms
# DECAGE - What age did decline begin
# FDGAD - FDG pattern of AD
# CSFTAU - elevated tau or ptau
# TAUPETAD - tau PET evidence of AD
# OTHBIOM - other biomarker found

scores <- c("NACCMMSE","PENTAGON","LOGIPREV","BOSTON","DIGIF","MOCATOTS")
# LOGIPREV = previous [MMSE] score. MOCATOTS = total MOCA raw scores
# DIGIF = Digit span forward trials correct

# comorbidities and health
comorbid <- c("CVHATT","HATTYEAR","CVCHF","CBSTROKE","NACCSTYR","PD","PDYR","TBI",
              "DIABETES","B12DEF","ARTHRIT","THYDIS","DEP2YRS","DEPOTHR","ANXIETY",
              "OCD")
# cvhatt - heart attack, hattyear - year of most recent; CVCHF - congestive heart failure
# NACCSTYR - most recent year of stroke. PD - parkinson's disease. TBI traumatic brain injury


# lifestyle (EXERCISE - DIET?)
lifestyle <- c("SMOKYRS","QUITSMOK","ALCOCCAS","ALCOHOL","ABUSOTHR","NACCBMI","ANYMEDS")

mind <- c("MEMORY","ORIENT","PERSCARE","SATIS","BORED")


all <- c(base,dxs,scores,comorbid,lifestyle,mind)
all %in% colnames(raw_data)
#View(raw_data[-rm_apoe,all])

process_data <- raw_data[-rm_apoe,all]


### how missing values are coded
# in general -4 is not available
# 9, 99, 999 is unknown
# 8, 88, 888 is no cognitive impairment

missing <- apply(process_data, 2, function(x) x==-4 | is.na(x))

# number of missing values by category
mis_base <- apply(missing[,base], 1, sum)
mis_dx <- apply(missing[,dxs], 1, sum)
mis_scores <- apply(missing[,scores], 1, sum)
mis_comorbid <- apply(missing[,comorbid], 1, sum)
mis_lifestyle <- apply(missing[,lifestyle], 1, sum)  
mis_mind <- apply(missing[,mind], 1, sum)  

num_vars <- c(length(dxs),length(scores),length(comorbid),length(lifestyle),
              length(mind))

## number of observations with all variabLes
sum(mis_mind==0)/nrow(process_data)

dxs2 <- setdiff(dxs,c("COGFLAGO","FDGAD","TAUPETAD","CSFTAU","OTHBIOM","OTHBIOMX"))
scores2 <- setdiff(scores,c("MOCATOTS","LOGIPREV"))
comorbid2 <- setdiff(comorbid,c("TBI","ARTHRIT","ANXIETY","OCD","THYDIS"))
lifestyle2 <- setdiff(lifestyle,c("ALCOCCAS"))
mind2 <- mind

all2 <- c(base,dxs2,scores2,comorbid2,lifestyle2,mind2)
process_data2 <- process_data[,all2]

## see proportion of AZ
vars43 <-c("CVHATT","CVCHF","CBSTROKE","PD","DIABETES","B12DEF","DEP2YRS",
           "DEPOTHR","SMOKYRS","ALCOHOL","ABUSOTHR")
n_subs <- length((process_data2$NACCID))
ad_tot <- length((process_data2$NACCID[process_data2$NACCALZD==1]))
ad_tot/n_subs

n_missing <- numeric(length(vars43))
ad_missing <-numeric(length(vars43))
prop_missing <- numeric(length(vars43))
for (i in 1:length(vars43)) {
  n_missing[i] <- length((process_data2$NACCID[missing[,vars43[i]]]))
  ad_missing[i] <-length((process_data2$NACCID[missing[,vars43[i]] & process_data2$NACCALZD==1]))
  prop_missing[i] <- ad_missing[i]/n_missing[i]

}
obs_missing <- apply(missing[,vars43], 1, sum)
# unique(obs_missing) = 0 11
## Same proportions and Missingness comes from the same observation
data_fullvars <- process_data2[obs_missing==0,]
data_reduced <- process_data2[,setdiff(all2,vars43)]

# > length(unique(data_reduced$NACCID))
# [1] 33816
# > length(unique(data_fullvars$NACCID))
# [1] 33816

### Outcome of script is processed data
## The missingness is at some visits, removing the missing observations preserves the subjects
# set no cognitive impairment to no AD, so remove observations
finalData <- data_fullvars[data_fullvars$NACCALZD!=8,]
vars <- list(base=base,dx = dxs2, scores=scores2,comorbid=comorbid2,
             lifestyle=lifestyle2, mind=mind2)


# rm(list=setdiff(ls(),c("raw_data","finalData","vars")))


