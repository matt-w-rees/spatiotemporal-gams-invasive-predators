## FIT GAMS TO MODEL DIEL ACTIVITY FOR FOXES AND FERAL CATS
# MATTHEW REES

# SETUP ------------------------------------------------------------
library(dplyr)
library(mgcv)
library(ggplot2)
library(patchwork)
library(viridis)
library(gratia)
library(chron)
library(sf)

rm(list = ls())
options(scipen = 999) 

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())


# PREPARE DATA ------------------------------------------------------------
## load gam data - 23 rows for each hour of every deployment
records <- read.csv("derived_data/counts_hour.csv")

# drop zoi's data 
records <- filter(records, data_source != "zoi")

# exclude stations left out for less than 7 days - something musta went wrong with those few cams
records <- filter(records, survey_duration >= 7)
hist(records$survey_duration / 23, breaks = 50)
summary(records$survey_duration)

# as we only had three "Riverine Grassy Woodlands or Forests" unique survey sites --> drop em
records <- filter(records, XGROUPNAME != "Riverine Grassy Woodlands or Forests")

# as we only had 20 "Rainforests" unique survey sites, all of which interspersed in "Wet Forests" in the Otways, reclassify as "Wet Forests". 
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Rainforests", "Wet or Damp Forests", as.character(records$XGROUPNAME))

# abbreviate EVC group names for plotting
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Wet or Damp Forests", "Wet Forests", records$XGROUPNAME)
records$XGROUPNAME <- if_else(records$XGROUPNAME == "Riparian Scrubs or Swampy Scrubs and Woodlands", "Swampy Scrubs", records$XGROUPNAME)
records$XGROUPNAME <- substr(records$XGROUPNAME, 1, nchar(records$XGROUPNAME) - 1)
table(records$XGROUPNAME)

# add xy coordinates
records_sf <- st_as_sf(records, coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = 32754)
records$x <- st_coordinates(records_sf)[,1]
records$y <- st_coordinates(records_sf)[,2]
head(records)


# ADD FOX COUNTS PER CAMERA -----------------------------------------------
## Single fox count for each cam-trap (and do the same for cats for good measure)
records <- records %>%
  group_by(station_year)  %>%
  mutate(fox_count = sum(fox),
         cat_count = sum(cat))

# fox presence / absence
records$fox_pa <- if_else(records$fox_count > 0, "present", "absent")
records$cat_pa <- if_else(records$cat_count > 0, "present", "absent")
xtabs(~records$fox_pa + records$cat_pa)

# adjust (positive only - we don't want negative counts) fox counts for survey effort 
records$fox_count_adj <- ifelse(records$fox_count > 0, records$fox_count/log(records$survey_duration), records$fox_count)

# and take the log (add 1 because you cant log 0)
records$fox_count_adj <- log(records$fox_count_adj + 1)
hist(records$fox_count_adj)

## Add habitat type cov (for models 2/3)
records$habitat_type <- if_else(records$XGROUPNAME == "Wet Forest", "wet_otway", "dry_otway")
records$habitat_type <- if_else(records$region == "glenelg", "dry_glenelg", records$habitat_type)
records$habitat_type <- if_else(records$region == "glenelg", "dry_glenelg", records$habitat_type)
records$habitat_type <- if_else(is.na(records$habitat_type), "dry_otway", records$habitat_type)
unique(records$habitat_type)

## Transform variable class  for GAMs
records <- transform(records,
                     hour = as.integer(hour),
                     cat = as.integer(cat), 
                     fox = as.integer(fox), 
                     foxbaits = as.integer(foxbaits), 
                     fox_pa = factor(fox_pa, ordered = FALSE), 
                     fox_count_adj = as.numeric(fox_count_adj), 
                     region = factor(region, ordered = FALSE), 
                     station = factor(station, ordered = FALSE), 
                     vegetation_group = factor(XGROUPNAME, ordered = FALSE), 
                     habitat_type = factor(habitat_type, ordered = FALSE), 
                     survey_duration = as.integer(survey_duration)
)
summary(records)

# save modified dataframe
write.csv(records, "derived_data/counts_hour_cleaned.csv")


# GAMS --------------------------------------------------------------------
# 1)  DO PREDATORS CHANGE DIEL ACTIVITY ACROSS VEGETATION TYPES?
# 1A) FOXES:
gam_fox <- bam(fox ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) + 
                     s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                     s(longitude, latitude, bs = "ds",  m = c(1, 0.5), k = 200) +
                     s(region, bs = "re") +  
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                    data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_fox)
plot(gam_fox, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)


# 1B) FERAL CATS:
gam_cat_1 <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                       s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                       s(longitude, latitude, bs = "ds",  m = c(1, 0.5), k = 200) +
                       s(region, bs = "re") +  
                       s(station, bs = "re") +  
                       offset(log(survey_duration)), 
                 data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_cat_1)
plot(gam_cat_1, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)


# 2b) fox counts
gam_cat_2b <- bam(cat ~ habitat_type + t2(hour, fox_count_adj, by = habitat_type, bs = c("cc", "ts"), k = c(8, 5), full = TRUE) +  
                    s(longitude, latitude, bs = "ds",  m = c(1, 0.5), k = 200) +
                    s(station, bs = "re") +  
                    offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_cat_2b)
plot(gam_cat_2b, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)


# 3 - space x time
gam_cat_3 <- bam(cat ~ t2(longitude, latitude, hour, d = c(2, 1), bs = c("tp", "cc"), k = c(80, 8), full = TRUE) +
                   s(station, bs = "re") +  
                   offset(log(survey_duration)), 
                 data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_cat_3)
plot(gam_cat_3, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)


gam_fox_3 <- bam(fox ~ t2(longitude, latitude, hour, d = c(2, 1), bs = c("tp", "cc"), k = c(80, 8), full = TRUE) +
                       s(foxbaits, bs = "tp", k = 4) + 
                       s(station, bs = "re") +  
                       offset(log(survey_duration)), 
                 data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)
summary(gam_fox_3)
plot(gam_fox_3, pages = 1, scheme = 2, seWithMean = TRUE, shade = TRUE, scale = 0)



# PLOT MODEL 1 ------------------------------------------------------------------
# get average sunrise / sunset times for plots
mean(times(records$sunrise))
mean(times(records$sunset))

# function to run and plot GAMs for each species
gam_plot <- function(gam_model, data){
  ## Plot global smooth (using the gratia package)
  plot_global <- draw(evaluate_smooth(gam_model, overall_uncertainty = TRUE, "s(hour)"), xlab = "Hour",  title = "", subtitle = "i.   Global function") + 
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") 
  ## Plot the difference smooths using the gratia package
  plot_difference <- draw(evaluate_smooth(gam_model, "s(hour,vegetation_group)"), xlab = "Hour",  title = "", subtitle = "ii.   Vegetation-specific deviations from the global function") + 
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") 
  ## Plot the actual hourly activity for each vegetation type   
  # make a new dataframe to predict into
  df <- expand.grid(hour = 0:23,
                    vegetation_group = levels(data$vegetation_group),
                    foxbaits = 0, 
                    survey_duration = mean(data$survey_duration))
  # get smooths to exlude
  x = sapply(gam_model$smooth, "[[",  "label")
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:6], "s(survey_duration)")))
  plot_veg <- ggplot(data=df, aes(x=hour, y=fit, group=vegetation_group)) +
    # ylim(-13, -2) +
    geom_line(aes(y=fit, x=hour), lwd = 0.7) +
    facet_wrap(~vegetation_group, nrow =1) +
    geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
    labs(subtitle = "iii.   Predicted activity", x = "Hour", y = "log(count)") 
  ## Assemble plots
  patch <- (plot_global | plot_difference) / plot_veg 
  patch
}

# plot / save
png("figs/fox_veg.png", width = 9, height = 5, res = 600, units = "in")
gam_plot(gam_fox, records) + plot_annotation(title = "A.   Fox")
dev.off()

png("figs/cat_veg.png", width = 9, height = 5, res = 600, units = "in")
gam_plot(gam_cat_1, records) + plot_annotation(title = "B.   Feral cat")
dev.off()

# to combine figures, type in the terminal (using imagemagick): 
# convert figs/fox_veg.png figs/cat_veg.png  -append figs/predator_veg.png


# PLOT MODEL 2B ------------------------------------------------------------------
df = expand.grid(hour = 0:23,
                 habitat_type = levels(records$habitat_type),
                 longitude = mean(records$longitude),
                 latitude = mean(records$latitude),
                 fox_count_adj = seq(min(records$fox_count_adj), max(records$fox_count_adj), length=150),
                 station = "T052",
                 survey_duration = mean(records$survey_duration))
# predict model results into dataframe
x = sapply(gam_cat_2b$smooth, "[[",  "label")
df_int <- cbind(df, predict(gam_cat_2b_st, newdata = df, se.fit = TRUE, type = "link", exclude = c(x[4:5], "s(survey_duration)")))
# plot
plot_range <- ggplot(aes(hour, fox_count_adj, fill = fit), data = df_int) +
  geom_tile() +
  facet_wrap(~habitat_type, nrow =1) +
  scale_fill_viridis("log(count)", option = "viridis") + 
  ggtitle("", subtitle = "Feral cat activity") +
  geom_vline(xintercept = c(6.16,18.34), colour = "grey75", size = 1) + 
  ylab("log(adjusted fox count)") + 
  xlab("Hour") 

# plot standard error 
plot_range_se <- ggplot(aes(hour, fox_count_adj, fill = se.fit), data = df_int) +
  geom_tile() +
  facet_wrap(~habitat_type, nrow =1) +
  scale_fill_viridis("Standard error", option = "viridis") + 
  ggtitle("", subtitle = "Uncertainty") +
  geom_vline(xintercept = c(6.16,18.34), colour = "grey75", size = 1) + 
  ylab("log(adjusted fox count)") + 
  xlab("Hour") 

# plot / save
png("figs/cat_fox_count.png", width = 9, height = 6.5, res = 600, units = "in")
plot_range / plot_range_se + plot_annotation(tag_levels = "A")
dev.off()

# 1B) min / max fox 
# min df 
df_int_min <- df_int %>%
  group_by(habitat_type) %>%
  slice_min(fox_count_adj, n = 1) %>%
  mutate(minmax = "min")
# max df 
df_int_max <- df_int %>%
  group_by(habitat_type) %>%
  slice_max(fox_count_adj, n = 1) %>%
  mutate(minmax = "max")
# join df
df_int_minmax <- rbind(df_int_min, df_int_max)
# plot
plot_wet_lowhigh <- ggplot(data=df_int_minmax, aes(x=hour, y=fit, group=minmax)) +
  geom_line(aes(y=fit, x=hour), lwd = 1.2) +
  facet_wrap(~habitat_type, nrow =1) +
  ggtitle("", subtitle = "Feral cat activity") +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  ylab("Effect") + 
  xlab("Hour")


# END