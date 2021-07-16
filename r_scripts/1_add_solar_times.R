## add average anchored solar time to species records

library(maptools)        
library(dplyr)
library(activity) 
library(lubridate)

rm(list = ls())


# load function to change clock time to radial time for overlap package
ClocktimeToTimeRad <- function(Clocktime,
                               timeformat = "%Y-%m-%d %H:%M:%S")
{
  DateTime2 <- strptime(as.character(Clocktime), format = timeformat, tz = "UTC")
  Time2     <- format(DateTime2, format = "%H:%M:%S", usetz = FALSE)
  Time.rad  <- (as.numeric(as.POSIXct(strptime(Time2, format = "%H:%M:%S", tz = "UTC"))) -
                  as.numeric(as.POSIXct(strptime("0", format = "%S", tz = "UTC")))) / 3600 * (pi/12)
  return(Time.rad)
}


# load data
records <- read.csv("raw_data/predator_records.csv")
camdata <- read.csv("raw_data/camdata.csv")
# remove stupid X's
camdata <- subset(camdata, select = -c(X, X1))
records <- subset(records, select = -c(X))
# merge files
records <- left_join(records, camdata, by = "station_year")

# date class
records$date_time <- ymd_hms(records$date_time, tz = "Australia/Brisbane")


# ADD SOLAR TIMES FOR DETECTIONS ------------------------------------------
# take a matrix of the coordinates
coords <- matrix(c(records$longitude, records$latitude), nrow = length(records$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 

# calculate the time of sunrise given the location and timing of each species detection 
sunriset_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunrise", POSIXct.out=TRUE))
# convert to radian time and add to the dataframe
sunrise_radian <- ClocktimeToTimeRad(sunriset_df$time) 
# calculate the time of sunset given the location and timing of each species detection 
sunsett_df <- as.data.frame(sunriset(coords, records$date_time, direction = "sunset", POSIXct.out=TRUE))
# convert to radian time and add to the dataframe
sunset_radian <- ClocktimeToTimeRad(sunsett_df$time) 
# save these as a separate matrix for the transtime function below 
anchors <- matrix(c(sunrise_radian, sunset_radian), nrow = length(sunset_radian))

# we also need to convert the detection time to radian time
det_radian <- ClocktimeToTimeRad(records$date_time)

# we can now calculate the average anchored time (time relative to both sunrise and sunset) for every detection
det_anchored <- transtime(det_radian, anchors)

# add to dataframe
records$aa_radian <- det_anchored
records$sunrise_radian <- sunrise_radian
records$sunset_radian <- sunset_radian
# convert back to clock time (new cols)
records$aa_clock       <- records$aa_radian * 3.81971
records$sunrise_clock <- records$sunrise_radian * 3.81971
records$sunset_clock  <- records$sunset_radian * 3.81971
head(records)
summary(head(records))

# save 
write.csv(records, "derived_data/records_solar.csv")