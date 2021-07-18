## add average anchored solar time to species records
# using methods from Vazquez et al. 2019 https://doi.org/10.1111/2041-210X.13290

rm(list = ls())
library(maptools)        
library(dplyr)
library(activity) 
library(lubridate)

# load function to change clock time to radial time for overlap package
source("r_scripts/function_clock_time_to_radian.R")

# load data
records <- read.csv("raw_data/predator_records.csv")
camdata <- read.csv("raw_data/camdata.csv")

# remove stupid X's
camdata <- subset(camdata, select = -c(X, X1))
records <- subset(records, select = -c(X))

# add camdata info to records 
records <- left_join(records, camdata, by = "station_year")

# take a matrix of the coordinates
coords <- matrix(c(records$longitude, records$latitude), nrow = length(records$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 

# specify detection time as a date class, as well as timezone (use Brisbane because it doesn't have daylight savings)
records$date_time <- ymd_hms(records$date_time, tz = "Australia/Brisbane")

# calculate the time of sunrise given the location and timing of each species detection 
sunriset_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunrise", POSIXct.out=TRUE))
# convert to radian time 
sunrise_radian <- ClocktimeToTimeRad(sunriset_df$time) 
# calculate the time of sunset given the location and timing of each species detection 
sunsett_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunset", POSIXct.out=TRUE))
# convert to radian 
sunset_radian <- ClocktimeToTimeRad(sunsett_df$time) 

# save these as a separate matrix for the transtime function below 
anchors <- matrix(c(sunrise_radian, sunset_radian), nrow = length(sunset_radian))

# we also need to convert the detection time to radian time
det_radian <- ClocktimeToTimeRad(records$date_time)

# we can now calculate the average anchored time (time relative to both sunrise and sunset) for every detection
det_anchored <- transtime(det_radian, anchors)

# add to dataframe
records$aa_radian <- det_anchored
# as well as in clock time (new cols)
records$aa_clock <- records$aa_radian * 3.81971

# looks good?
head(records)
summary(head(records))

# save 
write.csv(records, "derived_data/records_solar.csv")

# END