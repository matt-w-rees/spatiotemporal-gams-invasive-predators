## FORMAT RECORDS INTO COUNTS PER HOUR FOR EACH CAMERA-TRAP DEPLOYOMENT 

rm(list = ls())
library(readr)
library(dplyr)
library(reshape2)
library(lubridate)
library(maptools)     

# load function to change clock time to radial time for overlap package
source("r_scripts/function_clock_time_to_radian.R")

# load data
records <- read.csv("derived_data/records_solar_ind.csv")
camdata <- read.csv("raw_data/camdata.csv")


# DATAFRAME 1:  HOUR ----------------------------------------------------
# need to do this for summing
records$occ <- 1

# make solar hour an integer, of character class
records$hour <- as.integer(records$aa_clock)

# melt the dataframe
gam_data <- melt(records, id.var = colnames(records), measure.var = "occ")

# cast to derive count per year per station_year long format
gam_data = dcast(gam_data, station_year + hour ~ species, fill = 0, fun = sum)

# make a table with all unique combinations - (use camdata so we don't just take the sites with detections)
df <- expand.grid(station_year = unique(camdata$station_year), hour = 0:23)

# merge into gam_data so every hour for every station is provided (even with no predator detections) - fills missing values with NA
gam_data <- left_join(df, gam_data)
# change NA's to 0's
gam_data$fox <- ifelse(is.na(gam_data$fox), 0, gam_data$fox)
gam_data$cat <- ifelse(is.na(gam_data$cat), 0, gam_data$cat)
unique(is.na(gam_data)) # check it worked

# sort by station, hour
gam_data <- arrange(gam_data, station_year, hour)
head(gam_data)

## merge in site covariates
gam_data <- left_join(gam_data, camdata, by = "station_year")
head(gam_data)


# ADD SUNSET SUNRISE TIMES FOR STATIONS ------------------------------------------
# make date class
gam_data$date_start <- ymd(gam_data$date_start, tz = "Australia/Brisbane")
gam_data$date_end <- ymd(gam_data$date_end, tz = "Australia/Brisbane")
# take a matrix of the coordinates
coords <- matrix(c(gam_data$longitude, gam_data$latitude), nrow = length(gam_data$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 
# calculate the time of sunrise given the location and timing of each species detection - midway point for survey
sunriset_df <- as.data.frame(sunriset(coords, gam_data$date_start + days(as.integer(gam_data$survey_duration / 2)), direction = "sunrise", POSIXct.out=TRUE))
# calculate the time of sunset given the location and timing of each species detection - midway point for survey
sunsett_df <- as.data.frame(sunriset(coords, gam_data$date_start + days(as.integer(gam_data$survey_duration / 2)), direction = "sunset", POSIXct.out=TRUE))
# load hms r package here so it does't conflict with lubridate^
library(hms)
# add time only to dataframe using hms package
gam_data$sunrise <- hms::as_hms(lubridate::ymd_hms(sunriset_df$time, tz = "Australia/Brisbane"))
gam_data$sunset <- hms::as_hms(lubridate::ymd_hms(sunsett_df$time, tz = "Australia/Brisbane"))

# looks good?
head(gam_data)
gam_data <- subset(gam_data, select = -c(X, X1))

# save
write.csv(gam_data, "derived_data/counts_hour.csv")

# END