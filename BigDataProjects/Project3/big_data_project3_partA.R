
setwd("~/Berkeley Classes/2021Spring/Big data A public health perspective/Project 3")

rm(list=ls())

### read data and put in format for estimation ####
database <- read.csv("data-proj-3a.csv")

names(database)
summary(database)
head(database)

# change date to numeric
database$DateNum <- as.numeric(as.POSIXct(database$Date,format="%Y-%m-%d"))
# numeric date has seconds, change to daily units
database$DateNum <- floor(database$DateNum/60/60/24) 

# remove (set to 0) negative values
database$TargetValue[database$TargetValue<0] <- 0

#dates in DB
# as.Date(min(database$DateNum)) # 19-January-2020
# as.Date(max(database$DateNum)) # 10-June-2020, 140 days


#separate into confirmed cases and deaths
cases  <- database[database$Target=="ConfirmedCases",] #484820 cases, including 0
deaths <- database[database$Target=="Fatalities",] #484820 deaths, including 0
rm("database")

# Add cases by date and country (187 countries * 140 days)
total_cases_by_country <- cases %>%
  group_by(Country_Region,DateNum) %>% 
  summarise(Population=mean(Population),TotalCases = sum(TargetValue),Weight=mean(Weight))

total_cases <- cases %>%
  group_by(DateNum) %>% 
  summarise(Population=mean(Population),TotalCases = sum(TargetValue),Weight=mean(Weight))

# Add deaths by date and country (187 countries * 140 days)
total_deaths_by_country <- deaths %>%
  group_by(Country_Region,DateNum) %>% 
  summarise(Population=mean(Population),TotalCases = sum(TargetValue),Weight=mean(Weight))

total_deaths <- deaths %>%
  group_by(DateNum) %>% 
  summarise(Population=mean(Population),TotalCases = sum(TargetValue),Weight=mean(Weight))

### exploratory analysis and graphs ####
library(ggplot2)
library(zoo)

# total global cases
plot_cases <- ggplot(data=total_cases_by_country,
                     mapping = aes(x= as.Date(DateNum),y=TotalCases/Population, colour = factor(Country_Region))) +
              geom_line() + xlab("Date") + ylab("Total Cases (K)") + theme(legend.position = "none")
print(plot_cases)
# total global deaths
plot_deaths <- ggplot(data=total_deaths_by_country,
                     mapping = aes(x= as.Date(DateNum),y=TotalCases/1e3, colour = factor(Country_Region))) +
  geom_line() + xlab("Date") + ylab("Total deaths (K)") + theme(legend.position = "none")
print(plot_deaths)



# it seems there are outliers we can analyze with tsoutliers
library(tsoutliers)
casests <- ts(total_cases$TotalCases) #creates time series objects
plot.ts(casests)

outliers_cases <-tso(casests) #automatic procedure for outliers detection
plot(outliers_cases)
# an arima model might not seem appropiate


### data base format and apply SL ####
source("superlearner_fun.R")
total_cases$Interaction <- total_cases$DateNum*total_cases$Population/1e8
total_cases$Sin3 <- sin((total_cases$DateNum-18284)/(2*pi)*3)
total_cases$Sin5 <- sin((total_cases$DateNum-18284)/(2*pi)*5)
total_cases$Sin7 <- sin((total_cases$DateNum-18284)/(2*pi)*7)

total_deaths$Interaction <- total_deaths$DateNum*total_deaths$Population/1e8
total_deaths$Sin3 <- sin((total_deaths$DateNum-18284)/(2*pi)*3)
total_deaths$Sin5 <- sin((total_deaths$DateNum-18284)/(2*pi)*5)
total_deaths$Sin7 <- sin((total_deaths$DateNum-18284)/(2*pi)*7)

varnames <- c("DateNum","Population","Sin3","Sin5","Sin7")
outname  <- "TotalCases"

# Super Learner (commented out the "bad" CV schemes)
# cases_sl1 <- prediction_sl_origin(total_cases,outname = outname,varnames = varnames)
deaths_sl <- prediction_sl_origin(total_deaths,outname = outname,varnames = varnames)

cases_sl <- prediction_sl_window(total_cases,outname = outname,varnames = varnames)
# deaths_sl2 <- prediction_sl_window(total_deaths,outname = outname,varnames = varnames)


# cross validate the Super Learner
cases_sl_cv <- cross_validate_SL(total_cases,cases_sl$SuperLearner,
                                 outname,varnames)
death_sl_cv <- cross_validate_SL(total_deaths,deaths_sl$SuperLearner,
                                 outname,varnames)


# Proportion of variance explained (although no F theory to say "significance")
## MSE is different than the risk estimate
mse_cases <- 1-mean((cases_sl$Preds-total_cases$TotalCases)^2)/
              var(total_cases$TotalCases)
mse_deaths <- 1-mean((deaths_sl$Preds-total_deaths$TotalCases)^2)/
  var(total_deaths$TotalCases)


# Plots the fits
sl_cases_plot  <- data.frame(Prediction = cases_sl$Preds,Real=total_cases$TotalCases,
                             Date = as.Date(total_cases$DateNum))
sl_deaths_plot <- data.frame(Prediction = deaths_sl$Preds,Real=total_deaths$TotalCases,
                             Date = as.Date(total_deaths$DateNum))

ggplot(data=sl_cases_plot,aes(x=Date)) + geom_point(aes(y=Real)) +
  geom_line(aes(y=Prediction)) + ylab("Predicted vs. real cases")

ggplot(data=sl_deaths_plot,aes(x=Date)) + geom_point(aes(y=Real)) +
  geom_line(aes(y=Prediction)) + ylab("Predicted vs. real deaths")



