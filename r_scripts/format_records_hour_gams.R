# load r packages
library(maptools)        
library(dplyr)
library(activity) 
library(lubridate)
library(reshape2)
library(sf)
# note: chron package is loaded further down

# load function to change clock time to radial time for overlap package
source("r_scripts/function_clock_time_to_radian.R")


# LOAD / PREPARE DATA -----------------------------------------------------

# load data
records <- read.csv("raw_data/predator_records.csv")
camdata <- read.csv("raw_data/camdata.csv")
# remove stupid X's
camdata <- subset(camdata, select = -c(X))
records <- subset(records, select = -c(X))
# add camdata info to records - NA for sites without pred detections?
records <- left_join(records, camdata, by = c("station", "station_year"))


# 1) ADD SOLAR TIMES -------------------------------------------------------
## add average anchored solar time to species records
# using methods from Vazquez et al. 2019 https://doi.org/10.1111/2041-210X.13290

# take a matrix of the coordinates
coords <- matrix(c(records$longitude, records$latitude), nrow = length(records$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 
# specify detection time as a date class, as well as timezone (use Brisbane because it doesn't have daylight savings)
records$date_time <- ymd_hms(records$date_time, tz = "Australia/Brisbane")
records$date_start <- ymd(records$date_start, tz = "Australia/Brisbane")
records$date_end <- ymd(records$date_end, tz = "Australia/Brisbane")


## ADD ADJUSTED HOUR FOR EACH SPECIES DETECTION
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


## ADD AVERAGE SUNRISE SUNSET TIMES FOR EACH DEPLOYMENT
# calculate the time of sunrise given the location and timing of each species detection - midway point for survey
sunriset_df <- as.data.frame(sunriset(coords, records$date_start + days(as.integer(records$survey_duration / 2)), direction = "sunrise", POSIXct.out=TRUE))
# calculate the time of sunset given the location and timing of each species detection - midway point for survey
sunsett_df <- as.data.frame(sunriset(coords, records$date_start + days(as.integer(records$survey_duration / 2)), direction = "sunset", POSIXct.out=TRUE))
# load hms r package here so it does't conflict with lubridate^
library(hms)
# add time only to dataframe using hms package
records$sunrise_survey <- hms::as_hms(lubridate::ymd_hms(sunriset_df$time, tz = "Australia/Brisbane"))
records$sunset_survey <- hms::as_hms(lubridate::ymd_hms(sunsett_df$time, tz = "Australia/Brisbane"))
# looks good?
head(records)


# 2) REMOVE REPEAT DETECTIONS ---------------------------------------------
# (used lorelograms [in a separate script] to decide on 30 minutes)

# date_time as date class
records$date_time <- ymd_hms(records$date_time, tz = "Australia/Brisbane") # brisbane because we ignored daylight savings
## remove records of the same species at a particular cam-trap within 30 minutes
time_to_independence <- 1800 # number of seconds in 30 minutes
# get time (in seconds) since last visit
records <- records %>%
  group_by(station_year, species) %>%
  arrange(date_time) %>%
  mutate(seconds_diff = date_time - lag(date_time))
# change NA's to a value greater than the time to independence so they get kept when we filter out 
records$seconds_diff[is.na(records$seconds_diff)] <- time_to_independence + 1
# filter out records less than 60 minutes
records <- filter(records, seconds_diff > time_to_independence)
# remove seconds_diff col
records <- subset(records, select = -c(seconds_diff))
# looks good?
head(records)


# 3) RESHAPE COUNTS PER HOUR  -------------------------------------------------

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
## merge in site covariates
gam_data <- left_join(gam_data, camdata, by = "station_year")
head(gam_data)



# MODIFY DATA / VARIABLES --------------------------------------------------------

## EXCLUDE DATA
# drop zoi's data 
gam_data <- filter(gam_data, data_source != "zoi")
# exclude stations left out for less than 2 weeks 
gam_data <- filter(gam_data, survey_duration >= 10)
# as we only had three "Riverine Grassy Woodlands or Forests" unique survey sites --> drop em
gam_data <- filter(gam_data, XGROUPNAME != "Riverine Grassy Woodlands or Forests")

## MODIFY EVC GROUPS
# as we only had 20 "Rainforests" unique survey sites, all of which interspersed in "Wet Forests" in the Otways, reclassify as "Wet Forests". 
gam_data$XGROUPNAME <- if_else(gam_data$XGROUPNAME == "Rainforests", "Wet or Damp Forests", as.character(gam_data$XGROUPNAME))
# abbreviate EVC group names for plotting
gam_data$XGROUPNAME <- if_else(gam_data$XGROUPNAME == "Wet or Damp Forests", "Wet Forests", gam_data$XGROUPNAME)
gam_data$XGROUPNAME <- if_else(gam_data$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", "Swampy Scrubs", gam_data$XGROUPNAME)
gam_data$XGROUPNAME <- substr(gam_data$XGROUPNAME, 1, nchar(gam_data$XGROUPNAME) - 1)


## ADD VARIABLES

# add xy coordinates
gam_data_sf <- st_as_sf(gam_data, coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 32754)
gam_data$x <- st_coordinates(gam_data_sf)[,1]
gam_data$y <- st_coordinates(gam_data_sf)[,2]

# single fox count for each cam-trap (and do the same for cats for good measure)
gam_data <- gam_data %>%
  group_by(station_year)  %>%
  mutate(fox_count_total = sum(fox),
         cat_count_total = sum(cat))

# fox presence / absence
gam_data$fox_pa <- if_else(gam_data$fox_count_total > 0, "present", "absent")
gam_data$cat_pa <- if_else(gam_data$cat_count_total > 0, "present", "absent")
xtabs(~gam_data$fox_pa + gam_data$cat_pa)

## Add habitat type cov (for models 2/3)
gam_data$habitat_type <- if_else(gam_data$XGROUPNAME == "Wet Forest", "wet_otway", "dry_otway")
gam_data$habitat_type <- if_else(gam_data$region == "glenelg", "dry_glenelg", gam_data$habitat_type)
gam_data$habitat_type <- if_else(gam_data$region == "glenelg", "dry_glenelg", gam_data$habitat_type)
gam_data$habitat_type <- if_else(is.na(gam_data$habitat_type), "dry_otway", gam_data$habitat_type)
unique(gam_data$habitat_type)



## SAVE
# looks good but?
head(gam_data)
# now save
write.csv(gam_data, "derived_data/predator_counts_hour.csv")


# END