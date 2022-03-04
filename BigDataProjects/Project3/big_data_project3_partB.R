
### Quick peak at all datasets ####
# Look into all variables per subfolder to pre-select project
# rm(list = ls())
# 
# #directory with the folders with the cvs by source
 setwd("~/Berkeley Classes/2021Spring/Big data A public health perspective/Project 3")
 setwd("data-proj-3b/RDSC-07-30-Update")
# crtwd <- getwd()
# 
# # iterate over all folders and read only the header
# folders <- list.dirs()
# n_folders <- length(folders)
# # list with the name of the folder
# databases_folders <- vector("list",length = n_folders)
# names(databases_folders) <- sapply(c(1:n_folders),
#                                    function(x){strsplit(folders[x],"./")[[1]][2]})
# for(i in c(2:n_folders)){
# 
#   # read inside the folder
#   setwd(folders[i])
#   files <- list.files()
#   n_files <- length(files)
# 
#   # add all column names
#   for(j in c(1:n_files)){
#     cols_in_file <- colnames(read.csv(files[j],header = T,nrows = 1))
#     databases_folders[[i]] <-append(databases_folders[[i]],cols_in_file)
#   }
#   setwd(crtwd)
# }

# Potential projects to work on:
# - governement_of_mexico: Analyze the spread of the disease by date of entry to de country
# - nextstrain: Analyze how new strains are being spread and compare to the first one (and part a))
# - our_world_in_data: See new cases and deaths as function of tests, stringency?, gdp, cardiovascular disease, poverty,
# diabetes, hospitals_beds_per_thousand


### Better peak at the 3 potential problems to address
crtwd <- getwd()
#
# governement_of_mexico

# setwd("government_of_mexico")
# data_mexico <- read.csv("covid-19-mexico-daily-cases.csv",header = T,encoding = "UTF-8") # data with accents
# setwd(crtwd)
# 10K registries with date_of_arrival_in_mexico. All are confirmed or suspicious cases
# other information is age, sex, state and start of symptoms. 61 countries of origin, maybe less for different writting

# nextstrain
# setwd("nextstrain")
# data_strain <- read.csv("covid-19-genetic-phylogeny.csv",header = T)
# it just have all the 16284 strains found in 9 different species. it has dates, so I could study if it has sped up
# still couldnt address the confounding of simply having more tests. other would be country and country_exposure
# but hard to say what do this ones mean

setwd(crtwd)

# our_world_in_data
setwd("our_world_in_data")
data_world1 <- read.csv("coronavirus-disease-covid-19-statistics-and-research.csv",header = T)
# data_world2 <- read.csv("covid-19-testing-all-observations.csv",header = T)
# data_world3 <- read.csv("covid-19-testing-latest-data-and-source-details.csv",header = T)
setwd(crtwd)

#data_world1 seems to have the most complete info on cases, deaths, stringency index, tests, gdp and poverty...
# rm(list = setdiff(ls(),"data_world1"))

### Answer research question ####

# Impact of interventions on curving the impact of covid
# create database
library(ltmle)
#remove NAs in intervention

ids_AY <- (is.na(data_world1$new_cases_per_million) |
             is.na(data_world1$stringency_index))
            
data_clean <- data_world1[!ids_AY,]

# define variables
data_models <- data.frame("Country" = data_clean$location)
data_models$date <- data_clean$date

data_models$W1  <- as.numeric(as.POSIXct(data_clean$date,format="%Y-%m-%d"))
data_models$W1  <- floor(data_models$W1/60/60/24) # from seconds to days

data_models$W2  <- data_clean$new_tests_smoothed_per_thousand
data_models$W3  <- data_clean$population_density
data_models$W4  <- data_clean$median_age
data_models$W5  <- data_clean$aged_65_older
data_models$W6  <- data_clean$aged_70_older
data_models$W7  <- data_clean$gdp_per_capita
data_models$W8  <- data_clean$extreme_poverty
data_models$W9  <- data_clean$cardiovasc_death_rate
data_models$W10 <- data_clean$diabetes_prevalence
data_models$W11 <- (data_clean$female_smokers + data_clean$male_smokers)/2 #assumes 50%F/M
data_models$W12 <- data_clean$hospital_beds_per_thousand
data_models$W13 <- data_clean$life_expectancy
data_models$W14 <- data_clean$handwashing_facilities

data_models$stringency_id <- data_clean$stringency_index


# define intervention and outcome node
med_strin <- median(data_clean$stringency_index,na.rm = T)
data_models$A <- as.numeric(data_clean$stringency_index>=med_strin)
data_models$Y <- data_clean$new_cases_per_million

# quick function to clean NAs and remove variables with var = 0
clean_nas <- function(data_table){
  m <- ncol(data_table)
  ids <- NULL
  for(i in c(1:m)){
    ids_nas <- is.na(data_table[,i])
    drop_col <- sum(ids_nas)>=(nrow(data_table)/2)
    if(var(data_table[,i],na.rm = T)!=0||drop_col){
      
      in_med <- median(data_table[!ids_nas,i])
      if(in_med!=""){
        data_table[ids_nas,i] <- in_med
      }else{
        data_table[ids_nas,i] <- mode(data_table[!ids_nas,i])
      }
      
    }else{
      ids <- c(ids,i)
    }
    drop_col <- FALSE
  }
  if(length(ids)>0) {data_table <- data_table[,-ids]}
  return(data_table)
}


data_models <- cbind(data_models[,c(1,2)],clean_nas(data_models[,-c(1,2)]))
# remove the outlier and negative values
data_models <- data_models[data_models$Y!=max(data_models$Y),]
data_models$Y[data_models$Y<0] <- 0 

# add stringency index from previous week and 2 weeks
# compute week
data_models$week <- ceiling(data_models$W1/7)%%52 #I know it is less than a year
data_models$week <- data_models$week - min(data_models$week)+1
data_models$CntrWeek <- paste(data_models$Country,data_models$week)

data_reg <- data_models

# for weeks before
week_1 <- data_models %>% group_by(Country,week) %>% 
            summarise(WeekT1 = mean(stringency_id))

for (i in c(1,2)) {
  # add 1 week (and 2 weeks next cycle)
  week_1$week <- week_1$week+1 # to match with previous week
  week_1$CntrWeek <- paste(week_1$Country,week_1$week)
  
  data_reg <- merge(data_reg,week_1,by.x = "CntrWeek",
                       by.y="CntrWeek",all.x = T)
  
  colnames(data_reg)[which(colnames(data_reg)=="week.x")] <- "week"
  data_reg <- data_reg[,-which(colnames(data_reg)=="week.y")]
  colnames(week_1)[3] <- "weekT2"
}
#
#remove week 1 and 2 that dont have info on previous two weeks
data_reg <- data_reg[data_reg$week>2,]

# delete doubled up variables from the merging with the stringency averages
# and other character (eg country nAME)
ids_rep <- c(1,2,3,19,22,21,24) 
data_reg <- data_reg[,-ids_rep]
colnames(data_reg)

# NAs from missing data in original set (doubled check reasons for missigness before automating)
id_debug <- rep(FALSE,nrow(data_reg))
for(i in c(1:ncol(data_reg))){
  id_debug <- id_debug | is.na(data_reg[,i])
}
data_reg <- data_reg[!id_debug,]

# linear model
linear_model1 <- lm(Y~. ,data = data_reg)
mean((linear_model1$fitted.values-data_reg$Y)^2)/var(data_reg$Y)
#fit seems rather week ,try interactions or other models

tree_mod <- rpart(Y ~ .,data = data_reg)# control = rpart.control(cp=0.005))
fancyRpartPlot(tree_mod, main = "Cases/million", caption = NULL)

y_predict <- predict(tree_mod, data_reg)
mean((y_predict-data_reg$Y)^2)/var(data_reg$Y)

#biggest predictor is W2,W1, W3, W9, weekT2 (interact with W9,W2), stringency_id
colnames(data_reg)
ids_filtered <- c(1,2,3,9,18,15,16)
data_fil <- data_reg[,ids_filtered]
colnames(data_fil)

reg_form <- paste(colnames(data_fil)[-7],collapse = "+")
reg_form <- paste(c(reg_form,I("W2*W1"),I("stringency_id*W1"),I("stringency_id*W2"),
              I("W2*W9"),I("weekT2*W9"),I("W2*W9")),collapse = "+")
reg_form <- paste("Y~",reg_form,collapse = "")

data_strds <- scale(data_fil)

linear_model2 <- lm(reg_form,as.data.frame(data_strds))
linear_model2$coefficients

mean((linear_model2$fitted.values-data_reg$Y)^2)/var(data_reg$Y)



