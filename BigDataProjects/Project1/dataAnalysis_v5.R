
### Set directory and read tables ####
setwd("~/Berkeley Classes/2021Spring/Big data A public health perspective/Project I")
source('auxiliar_functions.R', echo=F)
library(future)
plan(multicore) #multiprocess

# load("~/Berkeley Classes/2021Spring/Big data A public health perspective/Project I/ProcessedTables.Rdata")
wd <- getwd()
setwd('Data')

# read data tables
data_claims  <- read.csv("Claims.csv")
data_daysY2  <- read.csv("DaysInHospital_Y2.csv")
data_daysY3  <- read.csv("DaysInHospital_Y3.csv")
data_drug    <- read.csv("DrugCount.csv")
data_lab     <- read.csv("LabCount.csv")
data_members <- read.csv("Members.csv")
data_target  <- read.csv("Target.csv")

# read look-up tables
#lookup_condition <- read.csv("Lookup PrimaryConditionGroup.csv")
#lookup_procedure <- read.csv("Lookup ProcedureGroup.csv")

setwd(wd)

### Data cleaning and preprocessing ####

#
## see missing data and NAs

# 20329 missing either age (5,753) or sex (17,552). Sex imputed based on likelihood given age
# then age imputed as average given sex

# First, set age to numeric to the beginning of the category interval
data_members$AgeAtFirstClaim <- as.numeric(gsub("([0-9]+).*$","\\1",data_members$AgeAtFirstClaim))

# input sex to likelihood-based general or on age 
bool_men <- data_members$Sex=="M"
bool_fem <- data_members$Sex=="F"
p <- sum(bool_men) /(sum(bool_men)+sum(bool_fem))

n <- sum(is.na(data_members$AgeAtFirstClaim))
men_no_age <- rbinom(n,1,p)

data_members$Sex[which(is.na(data_members$AgeAtFirstClaim))*men_no_age] <- "M"
data_members$Sex[which(is.na(data_members$AgeAtFirstClaim))*(1-men_no_age)] <- "F"

bool_id <- data_members$Sex==""
data_members$AgeAtFirstClaim[is.na(data_members$AgeAtFirstClaim)] <- -1
ages <- unique(data_members$AgeAtFirstClaim)
ages <- ages[ages!=-1]
for(age in ages){
  
  bool_age <- data_members$AgeAtFirstClaim==age
  men <- sum(bool_men & bool_age)
  fem <- sum(bool_fem & bool_age)
  p   <- men / (men+fem)
  
  n <- sum(bool_age)
  men_no_age <- rbinom(n,1,p)
  
  data_members$Sex[which(data_members$AgeAtFirstClaim==age)*men_no_age] <- "M"
  data_members$Sex[which(data_members$AgeAtFirstClaim==age)*(1-men_no_age)] <- "F"
}
# Input empty ages to the average given sex
m_avg <- mean(data_members$AgeAtFirstClaim[bool_men],na.rm = T)
f_avg <- mean(data_members$AgeAtFirstClaim[bool_fem],na.rm = T)

data_members$AgeAtFirstClaim[data_members$AgeAtFirstClaim==-1&bool_men] <- m_avg
data_members$AgeAtFirstClaim[data_members$AgeAtFirstClaim==-1&bool_fem] <- f_avg

#remove dummy variables for this segment
rm(list = c("bool_age","bool_fem","bool_id","bool_men","men","fem",
            "men_no_age","m_avg","f_avg","p","n"))

#
## Clean data encoding and coercing to numeric 
data_claims$LengthOfStay[data_claims$LengthOfStay==""] <- "0"
weeks_id <- as.numeric(grepl("*weeks*",data_claims$LengthOfStay))*7
weeks_id[weeks_id==0] <- 1
data_claims$LengthOfStay <- substr(data_claims$LengthOfStay,1,1)
data_claims$LengthOfStay <- as.numeric(data_claims$LengthOfStay)*weeks_id

data_claims$PayDelay[which(data_claims$PayDelay=='162+')] <- '162'
data_claims$PayDelay <- as.numeric(data_claims$PayDelay)

data_drug$DrugCount[which(data_drug$DrugCount=='7+')] <- '7'
data_drug$DrugCount <- as.numeric(data_drug$DrugCount)

data_lab$LabCount[which(data_lab$LabCount=='10+')] <- '10'
data_lab$LabCount <- as.numeric(data_lab$LabCount)


# set "" data as 0 in DSFS
data_claims$DSFS[data_claims$DSFS==""] <- "0"
data_drug$DSFS[data_drug$DSFS==""]     <- "0"
data_lab$DSFS[data_lab$DSFS==""]       <- "0"

# Clean DSFS code to the first month in interval
data_claims$DSFS <- as.numeric(gsub("([0-9]+).*$","\\1",data_claims$DSFS))
data_drug$DSFS   <- as.numeric(gsub("([0-9]+).*$","\\1",data_drug$DSFS))
data_lab$DSFS    <- as.numeric(gsub("([0-9]+).*$","\\1",data_lab$DSFS))


##
## informative: compute how many clients in training data have two consecutive
##  years of claims (and total) and for the data to predict

## over 71435 clients with claims on Y3
#two_years <- intersect(data_claims$MemberID[data_claims$Year=="Y1"],data_claims$MemberID[data_claims$Year=="Y2"])

## over 70942 members in target data Y4
#two_years_pred <- intersect(data_claims$MemberID[data_claims$Year=="Y2"],data_claims$MemberID[data_claims$Year=="Y3"])


#
## Create a wide data table(rows=MemberIDxYear,columns=All vars)
# table ordering is: MemberID,data_member,(DaysInHospital+target),data_drug,data_lab,claims

# build MemberID-Year combination. DIH Y_i are to be predicted by data Y_(i-1)
claims_yr2 <- cbind.data.frame(MemberID = data_daysY2$MemberID,
                               Truncated = data_daysY2$ClaimsTruncated,
                               DIH = data_daysY2$DaysInHospital,Year = "Y1")

claims_yr3 <- cbind.data.frame(MemberID = data_daysY3$MemberID,
                               Truncated = data_daysY3$ClaimsTruncated,
                               DIH = data_daysY3$DaysInHospital,Year = "Y2")

claims_yrt <- cbind.data.frame(MemberID = data_target$MemberID,
                               Truncated = data_target$ClaimsTruncated,
                               DIH = data_target$DaysInHospital,Year = "Y3")

claims_byIDYR <- rbind.data.frame(claims_yr2,claims_yr3,claims_yrt)

#change sex to wide in demographic table and add it to big table
wide_members <- num_claims_by_category(cbind(data_members,Year = "AUX"),"Sex",T)[-2]
final_table  <- merge(claims_byIDYR,wide_members,by = "MemberID", all.y = T)
final_table  <- final_table[,c(1,4,5,2,3)]

# add drug and laboratory to the table, label columns properly
colnames(data_drug)[3] <- "DrugDSFS"
colnames(data_lab)[3]  <- "LabDSFS"

# Aggregate numeric vectors by their mean
drug_num <- aggregate_by_vec(data_drug,c("DrugDSFS","DrugCount"))
lab_num <- aggregate_by_vec(data_lab,c("LabDSFS","LabCount"))

final_table <- merge(final_table,drug_num,by = c("MemberID","Year"), all = T)
final_table <- merge(final_table,lab_num,by = c("MemberID","Year"), all = T)

#NAs are induced by IDxYear that have no drug or lab, which really means zero
final_table[is.na(final_table)] <- 0

#
## process and add claims data table
# Sum numerical values describing each claim
vec_claim <- c("PayDelay","LengthOfStay","SupLOS","DSFS")
claims_num <- aggregate_by_vec(data_claims,vec_claim)


# Count total of claims
claims_total <- data_claims %>% count(MemberID, Year)
colnames(claims_total)[3] <- "NumClaims"

claims_total <- merge(claims_total,claims_num,by=c("MemberID","Year"),all = T)

###
# table of claims by specialty
categories <- c("Specialty","PlaceSvc","CharlsonIndex","ProcedureGroup","PrimaryConditionGroup")
table_claims_list <- build_wide_format_table(data_claims,categories,drop_var = T)

#change list of data tables to big table of data
claim_big_table <- list_to_data_frame(table_claims_list)
claim_big_table <- merge(claims_total,claim_big_table)

final_table <-  merge(final_table,claim_big_table,by = c("MemberID","Year"), all = T)
#save.image("DataProcess.Rdata")



### Variable selection and Modelling ####
library(rpart)
library(rattle)
library(rpart.plot)
library(RColorBrewer)

#clear environment
 rm(list = setdiff(ls(),"final_table"))
source('auxiliar_functions.R', echo=F)

# collapse white spaces for formula calling
colnames(final_table) <- gsub(" ","",colnames(final_table))

#
## separate into the modeling data and the data to predict, remove ID and year
ids_pred <- final_table$Year=="Y3"
modeling_data <- final_table[!ids_pred,-c(1,2)]
predict_data <- final_table[ids_pred,-c(1,2)]

# order DIH to be the first column for easier ID 
modeling_data <- modeling_data[,c(3,1,2,4:ncol(modeling_data))]
predict_data <- predict_data[,c(3,1,2,4:ncol(predict_data))]

#
## split model into training and testing samples (~70%/30%)  separated on DIH 0 or >0
# split into members with 0 and positive DIH
n_model <- nrow(modeling_data)
dih_zer <- modeling_data$DIH==0
dih_pos <- modeling_data$DIH>0

# Draw a IID sample of ~70% size of data and split in train and test set
train_zer <- rbinom(sum(dih_zer),1,0.7)
train_pos <- rbinom(sum(dih_pos),1,0.7)

train_ids <- numeric(n_model)

train_ids[dih_zer] <- dih_zer[dih_zer]*train_zer
train_ids[dih_pos] <- dih_pos[dih_pos]*train_pos

data_train <- modeling_data[as.logical(train_ids),]
data_test <-  modeling_data[as.logical(1-train_ids),]

#
## Variable screening considering only main and square terms

# create a big table with the main and quadratic terms
col_rmids <- c(1,2,3) # indicator variables are equal to the square value and remove outcome
data_train_sq <- data_train[,-col_rmids]^2 #remove indicator columns 
#give names of square variable
colnames(data_train_sq) <- paste(colnames(data_train)[-col_rmids],"^2",sep = "")

#create big table with main and square of variables
id_lins <- rep(T,ncol(data_train))
id_sqrs  <- rep(T,ncol(data_train_sq))

data_train_mix <- cbind.data.frame(data_train,data_train_sq)

id_lin_mix  <- c(id_lins,rep(F,length(id_sqrs)))
id_sqrs_mix <- c(rep(F,length(id_lins)),id_sqrs)
id_sqrs_mix[1] <- T #response var

# screen for linear and quadratic terms only and a 10% FDR
screen_linvals <- screen_lm(data_train_mix[,id_lin_mix],.1)
screen_sqvals  <- screen_lm(data_train_mix[,id_sqrs_mix],.1)

# retrieve ids for both runs and create vector of mixed terms (with heredity principle)
id_lins <- c(T,screen_linvals$SelVars) # add response
id_sqrs <- screen_sqvals$SelVars

id_lins[-col_rmids] <- (id_lins[-col_rmids] | id_sqrs) # heredity

# screen for final significant variables
screen_mix <- screen_lm(data_train_mix[,c(id_lins,id_sqrs)],cutoff = 0.10,FDR = T)


# retrieve variables into lin and square component
final_vars <- rep(F,length(c(id_lins,id_sqrs)))
final_vars[c(id_lins,id_sqrs)] <- c(T,screen_mix$SelVars)

id_lins <- final_vars[c(1:length(id_lins))]
id_sqrs <- final_vars[length(id_lins) + c(1:length(id_sqrs))]

#heredity principle, all square terms also include the main component
id_lins[-col_rmids] <- (id_lins[-col_rmids] | id_sqrs) 


# build new data bases for training, testing and predicting
training_final <- convert2screened(data_train,id_lins,id_sqrs,col_rmids)
testing_final  <- convert2screened(data_test,id_lins,id_sqrs,col_rmids)

predict_final <- convert2screened(predict_data,id_lins,id_sqrs,col_rmids)

# ## liner model with selected vars
fin_mod_lin <- lm(DIH ~.,data = training_final)

summa_lm <- get_rmse(predict(fin_mod_lin,testing_final),testing_final)

predDIH_lm <- predict(fin_mod_lin,predict_final)
sum(predDIH_lm)

# ## Build a tree with the selected variables
tree_mod <- rpart(DIH ~ .,data = training_final, control = rpart.control(cp=0.005))
fancyRpartPlot(tree_mod, main = "Days as categories", caption = NULL)

preds_trday <- predict(tree_mod,testing_final)

predDIH_tree <- predict(tree_mod,predict_final)

summa_tree <- get_rmse(preds_trday,testing_final)

### SECOND INVESTIGATION QUESTION ####
# Determine whether hospitalization days changed for men and women
m_hosp <- 
m_Nohosp <- sum(modeling_data$M == 1 & modeling_data$DIH==0)

f_hosp <- sum(modeling_data$M == 0 & modeling_data$DIH>0)
f_Nohosp <- sum(modeling_data$M == 0 & modeling_data$DIH==0)


expected_mat <- matrix(NA,2,2)
colnames(expected_mat) <- c("M","F")
rownames(expected_mat) <- c("Hosp","NoHosp")

observed_mat <- expected_mat
observed_mat[1,1] <- sum(modeling_data$M == 1 & modeling_data$DIH>0)
observed_mat[1,2] <- sum(modeling_data$M == 0 & modeling_data$DIH>0)
observed_mat[2,1] <- sum(modeling_data$M == 1 & modeling_data$DIH==0)
observed_mat[2,2] <- sum(modeling_data$M == 1 & modeling_data$DIH==0)

n_all <- nrow(modeling_data)/1000 #to avoid integer overflow
n_men <- sum(modeling_data$M == 1)/1000
n_hosp <- sum(modeling_data$DIH>0)/1000

expected_mat[1,1] <- n_men*n_hosp/n_all
expected_mat[1,2] <- (n_all-n_men)*n_hosp/n_all
expected_mat[2,1] <- n_men*(n_all-n_hosp)/n_all
expected_mat[2,2] <- (n_all-n_men)*(n_all-n_hosp)/n_all
expected_mat <- expected_mat*1000

chi_stats <- sum((expected_mat-observed_mat)^2/expected_mat)

p_val <- pchisq(chi_stats,1,lower.tail = F)
