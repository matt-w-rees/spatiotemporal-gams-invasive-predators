# Identify the minimum interval to approximate independence
# code and methods adapated from Iannarilli et al. (2019).

#' Load lorelogram package
if (!require("devtools")) install.packages("devtools", 
                                           repos = "http://cran.us.r-project.org", 
                                           dependencies = "Imports")
install.packages("usethis")
library(usethis)
devtools::install_github("FabiolaIannarilli/lorelogram") 
library(lorelogram)
library(camtrapR)
library(dplyr)

# load updated functions from Iannarilli et al. (2019).
source("r_scripts/Modified_function_for_detection_histories_by_minute.R")

# load data
records <- read.csv("derived_data/records_solar.csv")
camdata <- read.csv("raw_data/camdata.csv")


## FILTER TO 1 GRID DEPLOYMENT AT A TIME  (memory limit exhausts / R crashes when you try to do any more)
records <- filter(records, data_source == "matt" & block == "annya")
camdata <- filter(camdata, data_source == "matt" & block == "annya")

# make camera operation matrix using camtrapR package
camop <- cameraOperation(CTtable = camdata, stationCol = "station_year", setupCol = "date_start", retrievalCol = "date_end", writecsv = FALSE, hasProblems = FALSE, dateFormat = "%Y-%m-%d")
write.csv(camop, "derived_data/camop.csv") # have to save and reload to get it in the right format for the function below
camop_problem <- read.csv("derived_data/camop.csv", header=TRUE, check.names = FALSE) 
row.names(camop_problem) <- camop_problem[,1]
camop_problem <- camop_problem[,-1] #to remove first column with station name 

# run modified dethist function - written by Iannarilli et al. (2019)
DetHist_Minute_cat <- myfunc_detectionHistory_minute(recordTable = records,
                                                     species = "cat", #change species name
                                                     camOp = camop_problem,
                                                     stationCol =  "station_year",
                                                     speciesCol = "species",
                                                     recordDateTimeCol = "date_time", 
                                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                                     occasionLength = 10, 
                                                     includeEffort = FALSE,
                                                     scaleEffort = FALSE, 
                                                     day1 = "survey", 
                                                     timeZone="Australia/Brisbane") #adjust time zone

DetHist_Minute_fox <- myfunc_detectionHistory_minute(recordTable = records,
                                                     species = "fox", #change species name
                                                     camOp = camop_problem,
                                                     stationCol =  "station_year",
                                                     speciesCol = "species",
                                                     recordDateTimeCol = "date_time", 
                                                     recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                                     occasionLength = 10, 
                                                     includeEffort = FALSE,
                                                     scaleEffort = FALSE, 
                                                     day1 = "survey", 
                                                     timeZone="Australia/Brisbane") #adjust time zone

#write.csv(x = DetHist_Minute_cat$detection_history, file = "derived_data/DetHist_Minute.csv")
#DetHist_Minute <- read.csv("derived_data/DetHist_Minute.csv")

LOR_Stat_data_cat <- lorelogram(as.data.frame(DetHist_Minute_cat$detection_history), max_lag = 60)
lor_lag_to_indep(LOR_Stat_data_cat)

LOR_Stat_data_fox <- lorelogram(as.data.frame(DetHist_Minute_fox$detection_history), max_lag = 60)
lor_lag_to_indep(LOR_Stat_data_fox)

# 30 minutes is sufficient

# END