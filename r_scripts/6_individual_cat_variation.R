## REFORMAT CAT IDENTIFIED RECORDS FOR DIEL ACTIVITY GAMS

rm(list = ls())
library(stringr)
library(maptools)        
library(dplyr)
library(activity) 
library(lubridate)
library(reshape2)
library(mgcv)
library(gratia)
library(ggplot2)
library(patchwork)

# load function to change clock time to radial time for overlap package
source("r_scripts/function_clock_time_to_radian.R")

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# load cat records
records <- read.csv("raw_data/cat_id_detections.csv")
records$X <- NULL

# load camdata - but just matts deployments
camdata <- read.csv("raw_data/camdata.csv")
camdata <- filter(camdata, data_source == "matt")
camdata$X <- NULL
camdata$X.1 <- NULL



# ADD SOLAR TIMES ---------------------------------------------------------
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


# REMOVE REPEAT DETECTIONS ------------------------------------------------

## remove records of the same species at a particular cam-trap within 30 minutes
time_to_independence <- 1800 # number of seconds in 30 minutes

# get time (in seconds) since last visit
records <- records %>%
  group_by(station_year, individual) %>%
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



# SPLIT COMBINED CATS ---------------------------------------------------------
# duplicate row and split individual name
records_double1 <- filter(records, individual == "otway_swirl_venus_&_otway_swirl_peter" | individual == "annya_swirl_ben_&_annya_swirl_motley")
records_double2 <- filter(records, individual == "otway_swirl_venus_&_otway_swirl_peter" | individual == "annya_swirl_ben_&_annya_swirl_motley")
records_double1$individual <- str_split_fixed(records_double1$individual, "_&_", 2)[,1]
records_double2$individual <- str_split_fixed(records_double2$individual, "_&_", 2)[,2]
records_double <- rbind(records_double1, records_double2)

# remove these rows from records
records <- filter(records, is.na(individual) | individual != "otway_swirl_venus_&_otway_swirl_peter")
records <- filter(records, is.na(individual) | individual != "annya_swirl_ben_&_annya_swirl_motley")

# now add the new ones in 
records <- rbind(records, records_double)


# RESHAPE BY HOUR ---------------------------------------------------------
# DATAFRAME 1:  HOUR ----------------------------------------------------
# need to do this for summing
records$occ <- 1

# make solar hour an integer, of character class
records$hour <- as.integer(records$aa_clock)
records$species <- "cat"

# split by region
records_g <- filter(records, region == "glenelg")
records_o <- filter(records, region == "otways")
# do same for camdata
camdata_g <- filter(camdata, region == "glenelg")
camdata_o <- filter(camdata, region == "otways")

# melt the dataframe
gam_data_g <- melt(records_g, id.var = colnames(records_g), measure.var = "occ")
gam_data_o <- melt(records_o, id.var = colnames(records_o), measure.var = "occ")

# cast to derive count per year per station_year long format
gam_data_g = dcast(gam_data_g,  individual + station_year + hour ~ species, fill = 0, fun = sum)
gam_data_o = dcast(gam_data_o,  individual + station_year + hour ~ species, fill = 0, fun = sum)

# make a table with all unique combinations - (use camdata so we don't just take the sites with detections)
df_g <- expand.grid(station_year = unique(camdata_g$station_year), individual = unique(records_g$individual), hour = 0:23)
df_o <- expand.grid(station_year = unique(camdata_o$station_year), individual = unique(records_o$individual), hour = 0:23)

# merge into gam_data so every hour for every station is provided (even with no predator detections) - fills missing values with NA
gam_data_g <- left_join(df_g, gam_data_g)
gam_data_o <- left_join(df_o, gam_data_o)

# change NA's to 0's
gam_data_g$cat <- ifelse(is.na(gam_data_g$cat), 0, gam_data_g$cat)
gam_data_o$cat <- ifelse(is.na(gam_data_o$cat), 0, gam_data_o$cat)

# sort by station, hour
gam_data_g <- arrange(gam_data_g, individual, station_year, hour)
gam_data_o <- arrange(gam_data_o, individual, station_year, hour)

## merge in site covariates
gam_data_g <- left_join(gam_data_g, camdata_g, by = "station_year")
gam_data_o <- left_join(gam_data_o, camdata_o, by = "station_year")





# GAMS --------------------------------------------------------------------

## OTWAYS
## Transform variable class  for GAMs
gam_data_o <- transform(gam_data_o,
                       hour = as.integer(hour),
                       individual = factor(individual, ordered = FALSE), 
                       station = factor(station, ordered = FALSE), 
                       survey_duration = as.integer(survey_duration)
)

gam_cat_ind_o <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                           s(hour, individual, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                           s(longitude, latitude, bs = "ds",  m = c(1, 0.5), k = 100) +
                           s(station, bs = "re") +  
                           offset(log(survey_duration)), 
                     data = gam_data_o, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_cat_ind_o)
plot(gam_cat_ind_o, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)
draw(gam_cat_ind_o)



## GLENELG
## Transform variable class  for GAMs
gam_data_g <- transform(gam_data_g,
                        hour = as.integer(hour),
                        individual = factor(individual, ordered = FALSE), 
                        station = factor(station, ordered = FALSE), 
                        survey_duration = as.integer(survey_duration)
)

gam_cat_ind_g <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                           s(hour, individual, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                           s(longitude, latitude, bs = "ds",  m = c(1, 0.5), k = 150) +
                           s(station, bs = "re") +  
                           offset(log(survey_duration)), 
                     data = gam_data_g, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_cat_ind_g)
plot(gam_cat_ind_g, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)
draw(gam_cat_ind_g)



# ADD SUNSET SUNRISE TIMES FOR STATIONS ------------------------------------------
## glenelg 
# make date class
gam_data_g$date_start <- ymd(gam_data_g$date_start, tz = "Australia/Brisbane")
gam_data_g$date_end <- ymd(gam_data_g$date_end, tz = "Australia/Brisbane")
# take a matrix of the coordinates
coords <- matrix(c(gam_data_g$longitude, gam_data_g$latitude), nrow = length(gam_data_g$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 
# calculate the time of sunrise given the location and timing of each species detection - midway point for survey
sunriset_df <- as.data.frame(sunriset(coords, gam_data_g$date_start + days(as.integer(gam_data_g$survey_duration / 2)), direction = "sunrise", POSIXct.out=TRUE))
# calculate the time of sunset given the location and timing of each species detection - midway point for survey
sunsett_df <- as.data.frame(sunriset(coords, gam_data_g$date_start + days(as.integer(gam_data_g$survey_duration / 2)), direction = "sunset", POSIXct.out=TRUE))
# load hms r package here so it does't conflict with lubridate^
library(hms)
# add time only to dataframe using hms package
gam_data_g$sunrise <- hms::as_hms(lubridate::ymd_hms(sunriset_df$time, tz = "Australia/Brisbane"))
gam_data_g$sunset <- hms::as_hms(lubridate::ymd_hms(sunsett_df$time, tz = "Australia/Brisbane"))

# otways
# make date class
gam_data_o$date_start <- ymd(gam_data_o$date_start, tz = "Australia/Brisbane")
gam_data_o$date_end <- ymd(gam_data_o$date_end, tz = "Australia/Brisbane")
# take a matrix of the coordinates
coords <- matrix(c(gam_data_o$longitude, gam_data_o$latitude), nrow = length(gam_data_o$latitude))
# add a coordinate reference system (can't be UTM)
coords <- sp::SpatialPoints(coords, proj4string=sp::CRS("+proj=longlat +datum=WGS84")) 
# calculate the time of sunrise given the location and timing of each species detection - midway point for survey
sunriset_df <- as.data.frame(sunriset(coords, gam_data_o$date_start + days(as.integer(gam_data_o$survey_duration / 2)), direction = "sunrise", POSIXct.out=TRUE))
# calculate the time of sunset given the location and timing of each species detection - midway point for survey
sunsett_df <- as.data.frame(sunriset(coords, gam_data_o$date_start + days(as.integer(gam_data_o$survey_duration / 2)), direction = "sunset", POSIXct.out=TRUE))
# add time only to dataframe using hms package
gam_data_o$sunrise <- hms::as_hms(lubridate::ymd_hms(sunriset_df$time, tz = "Australia/Brisbane"))
gam_data_o$sunset <- hms::as_hms(lubridate::ymd_hms(sunsett_df$time, tz = "Australia/Brisbane"))

# get averages (use in plot)
library(chron)
mean(times(gam_data_g$sunrise))
mean(times(gam_data_g$sunset))
mean(times(gam_data_o$sunrise))
mean(times(gam_data_o$sunset))


# PLOTS -------------------------------------------------------------------
## 1st get top 7 cats with most detections
top7_g <- slice_max(as.data.frame(table(records_g$individual)), order_by = Freq, n = 7, with_ties = FALSE)
top7_o <- slice_max(as.data.frame(table(records_o$individual)), order_by = Freq, n = 7, with_ties = FALSE)

## glenelg
## Plot global smooth (using the gratia package)
plot_global <- draw(evaluate_smooth(gam_cat_ind_g, overall_uncertainty = TRUE, "s(hour)"), xlab = "Hour",  title = "", subtitle = "i.   Global function") + 
  geom_vline(xintercept = c(7.01,18.12), colour = "black", size = 0.6, linetype="dotted") 
## Plot the difference smooths using the gratia package
plot_difference <- draw(evaluate_smooth(gam_cat_ind_g, "s(hour,individual)"), xlab = "Hour",  title = "", subtitle = "ii.   Individual deviations from the global function") + 
  geom_vline(xintercept = c(7.01,18.12), colour = "black", size = 0.6, linetype="dotted") 
## Plot the actual hourly activity for each vegetation type   
# make a new dataframe to predict into
df <- expand.grid(hour = 0:23,
                  individual = unique(top7_g$Var1),
                  survey_duration = mean(gam_data_g$survey_duration))
# get smooths to exlude
x = sapply(gam_cat_ind_g$smooth, "[[",  "label")
# predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
df <- cbind(df, predict.gam(gam_cat_ind_g, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:6], "s(survey_duration)")))
# shorten cat names
df$individual <- as.character(df$individual)
df$individual <- recode(df$individual, annya_swirl_lola = "faith", c_tabby_quasimodo = "brisket", m_swirl_socksy = "boots", annya_swirl_ben = "dolly", annya_swirl_motley = "motley", hotspur_tabby_1 = "sprout", m_tabby_murray = "murray")
# plot top 7 individuals
plot_ind <- ggplot(data=df, aes(x=hour, y=fit, group=individual)) +
  geom_line(aes(y=fit, x=hour), lwd = 0.7) +
  facet_wrap(~individual, nrow =1) +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
  geom_vline(xintercept = c(7.01,18.12), colour = "black", size = 0.6, linetype="dotted") + 
  labs(subtitle = "iii.   Predicted activity", x = "Hour", y = "log(count)") 

# plot / save
png("figs/cat_ind_g.png", width = 9, height = 5, res = 600, units = "in")
(plot_global | plot_difference) / plot_ind  + plot_annotation(title = "A.   Glenelg region")
dev.off()


## otways
## Plot global smooth (using the gratia package)
plot_global <- draw(evaluate_smooth(gam_cat_ind_o, overall_uncertainty = TRUE, "s(hour)"), xlab = "Hour",  title = "", subtitle = "i.   Global function") + 
  geom_vline(xintercept = c(7.27,17.38), colour = "black", size = 0.6, linetype="dotted") 
## Plot the difference smooths using the gratia package
plot_difference <- draw(evaluate_smooth(gam_cat_ind_o, "s(hour,individual)"), xlab = "Hour",  title = "", subtitle = "ii.   Individual deviations from the global function") + 
  geom_vline(xintercept = c(7.27,17.38), colour = "black", size = 0.6, linetype="dotted") 
## Plot the actual hourly activity for each vegetation type   
# make a new dataframe to predict into
df <- expand.grid(hour = 0:23,
                  individual = unique(top7_o$Var1),
                  survey_duration = mean(gam_data_o$survey_duration))
# get smooths to exlude
x = sapply(gam_cat_ind_o$smooth, "[[",  "label")
# predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
df <- cbind(df, predict.gam(gam_cat_ind_o, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:6], "s(survey_duration)")))
# shorten cat names
df$individual <- as.character(df$individual)
df$individual <- recode(df$individual, otway_tabby_javier = "pav", otway_swirl_skip = "skip", otway_tabby_crass = "crass", otway_tabby_meg = "pauline", otway_swirl_bluey = "bluey", otway_swirl_chowder = "chowder", otway_tabby_catsup = "chen_kenichi")
# plot top 7 individuals
plot_ind <- ggplot(data=df, aes(x=hour, y=fit, group=individual)) +
  geom_line(aes(y=fit, x=hour), lwd = 0.7) +
  facet_wrap(~individual, nrow =1) +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
  geom_vline(xintercept = c(7.27,17.38), colour = "black", size = 0.6, linetype="dotted") + 
  labs(subtitle = "iii.   Predicted activity", x = "Hour", y = "log(count)") 

# plot / save
png("figs/cat_ind_o.png", width = 9, height = 5, res = 600, units = "in")
(plot_global | plot_difference) / plot_ind  + plot_annotation(title = "B.   Western Otway Ranges")
dev.off()

# to combine figures, type in the terminal (using imagemagick): 
# convert figs/cat_ind_g.png figs/cat_ind_o.png  -append figs/cat_ind.png

# END