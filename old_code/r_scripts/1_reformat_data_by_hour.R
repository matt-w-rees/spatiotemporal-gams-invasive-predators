# This script reformats data from a row per detection, to a row per hour at a site with counts of fox and cat detections
# Also converts detection times to solar time and adds human footprint index variable 


# load r packages
library(maptools)        
library(dplyr)
library(activity) 
library(lubridate)
library(reshape2)
library(sf)
library(hms)
library(terra)

# set terra tempdir (for cleanup) + raster progress bar
terraOptions(tempdir ="temp/", progress = 1, memfrac = 0.8)

# load function to change clock time to radial time for overlap package
source("r_scripts/function_clock_time_to_radian.R")



# LOAD / PREPARE DATA -----------------------------------------------------
# load human footprint index raster (1km)
hfi <- terra::rast("raw_data/wildareas-v3-2009-human-footprint-geotiff/wildareas-v3-2009-human-footprint.tif")

# records of species detections (1 row per detection)
records <- read.csv("raw_data/predator_records.csv")
# camera-trap deplopyment / site information 
camdata <- read.csv("raw_data/camdata.csv")
# stitch together
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

# add time only to dataframe using hms package
records$sunrise_survey <- hms::as_hms(lubridate::ymd_hms(sunriset_df$time, tz = "Australia/Brisbane"))
records$sunset_survey <- hms::as_hms(lubridate::ymd_hms(sunsett_df$time, tz = "Australia/Brisbane"))



# 2) REMOVE REPEAT DETECTIONS ---------------------------------------------
# (previously used lorelograms to decide on 30 minute window)

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



# ADD HUMAN FOOTPRINT INDEX TO DATASET -----------------------------------------------
# 128 == NA value: replace
hfi <- subst(hfi, 128, NA)

# make cam locations spatial (to relate to raster) and reproject to same crs as hfi 
gam_data_sp <- sf::st_as_sf(gam_data, coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = "ESRI:54009") 

# get hfi values for each camera-trap site - accounting for adjacent 1km cells too (due to cat homerange size and proximity to the ocean)
vals <- terra::extract(hfi, terra::vect(gam_data_sp), method = "bilinear")

# add these values to original gam_data
gam_data$hfi <- vals$`wildareas-v3-2009-human-footprint`



# MODIFY DATA / VARIABLES FOR MODELLING --------------------------------------------------------
## EXCLUDE DATA
# drop zoi's data, exclude stations left out for less than 2 weeks and as we only had three "Riverine Grassy Woodlands or Forests" unique survey sites --> drop em
gam_data <- filter(gam_data, data_source != "zoi" & survey_duration >= 10 & XGROUPNAME != "Riverine Grassy Woodlands or Forests")

## MODIFY EVC GROUPS
# as we only had 20 "Rainforests" unique survey sites, all of which interspersed in "Wet Forests" in the Otways, reclassify as "Wet Forests".
# abbreviate EVC group names for plotting
gam_data$XGROUPNAME <- if_else(gam_data$XGROUPNAME == "Rainforests" | gam_data$XGROUPNAME == "Wet or Damp Forests", "Wet Forests", 
                        if_else(gam_data$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", "Swampy Scrubs", gam_data$XGROUPNAME))
# and remove the plurarity (for plotting)
gam_data$XGROUPNAME <- substr(gam_data$XGROUPNAME, 1, nchar(gam_data$XGROUPNAME) - 1)


## ADD VARIABLES
# add habitat type cov (for models 2/3)
gam_data$habitat_type <- if_else(gam_data$XGROUPNAME == "Wet Forest", "Wet forest (Otway Ranges)", "Dry vegetation (Otway Ranges)")
gam_data$habitat_type <- if_else(gam_data$region == "glenelg", "Dry vegetation (Glenelg region)", gam_data$habitat_type)

# add xy coordinates -- vicgrid crs 
gam_data_sf <- st_as_sf(gam_data, coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 3111) 
gam_data$x <- st_coordinates(gam_data_sf)[,1]
gam_data$y <- st_coordinates(gam_data_sf)[,2]

# fox summaries for each cam-trap 
gam_data <- gam_data %>%
  group_by(station_year) %>%
  mutate(fox_count_total = sum(fox),
         fox_count_adj = log((fox_count_total / survey_duration) + 1), # new variable for fox counts adjusted for survey effort - cats more likely care about fox activity on the log scale
         fox_pa = if_else(fox_count_total > 0, "present", "absent")) %>%
  ungroup() 


# remove unneccesary cols - clean it up a bit 
gam_data_clean <- dplyr::select(gam_data, 
                             data_source, region, block, 
                             station, longitude, latitude, x, y, 
                             station_year, date_start, date_end, survey_duration, 
                             hour, cat, fox, fox_count_total, fox_count_adj, 
                             XGROUPNAME, habitat_type, foxbaits, hfi)


### SAVE DATA
write.csv(gam_data, "derived_data/predator_counts_hour.csv")


# END 