# REMOVE REPEAT OBSERVATIONS FOR SPECIES WIHIN 30 MINUTES AT A CAMERA FOR APPROXIMATE TEMPORAL INDEPENDENCE

rm(list = ls())
library(dplyr)
library(lubridate)

# load records
records <- read.csv("derived_data/records_solar.csv")

# date_time as date class
records$date_time <- ymd_hms(records$date_time, tz = "Australia/Brisbane")

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
records$seconds_diff <- NULL

# looks good?
head(records)

# save
write.csv(records, "derived_data/records_solar_ind.csv")

# END