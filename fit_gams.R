## FIT GAMS TO MODEL DIEL ACTIVITY FOR FOXES AND FERAL CATS
# CODE FOR OIKOS PAPER 'DYNAMIC SHIFTS IN PREDATOR DIEL ACTIVITY PATTERNS ACROSS LANDSCAPES AND THREAT-LEVELS'
# MATTHEW REES


# SET-UP ------------------------------------------------------------------

# load r-packages (all except for mgcv and dplyr are for plotting)
library(dplyr)
library(mgcv)
library(gratia)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(sp)
library(sf)
library(viridis)
library(forcats)

# load plotting functions
source("plotting_functions/function_plot_spatiotemporal_interaction.R")
source("plotting_functions/function_plot_vegetation_type.R")
source("plotting_functions/function_plot_fox_human_footprint.R")
source("plotting_functions/function_plot_cat_v_fox.R")
source("plotting_functions/function_plot_fox_foxbaits.R")

# load data: response variable = fox / cat detection counts at a camera-trap site ('station') for each hour of the day (0-23)
records <- read.csv("data_predator_counts_hour.csv")

# set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())


# Transform variables -----------------------------------------------------
records <- transform(records,
                     region = factor(region, ordered = FALSE), 
                     station = factor(station, ordered = FALSE), 
                     x = as.numeric(x),
                     y = as.numeric(y),
                     hour = as.integer(hour),
                     cat = as.integer(cat), 
                     fox = as.integer(fox), 
                     fox_count_adj = as.numeric(fox_count_adj), 
                     foxbaits = as.integer(foxbaits), 
                     vegetation_group = factor(vegetation_group, ordered = FALSE), 
                     habitat_type = factor(habitat_type, ordered = FALSE), 
                     survey_duration = as.integer(survey_duration),
                     hfi = round(hfi, digits = 2))




# 1) Interaction between spatial and diel activity ('model 1') ---------------------

# specify penalties for 2D Duchon spline and 1D cc spline
m <- list(c(1,.5),rep(0,0)) 

# fox model 1 
gam_fox_sp <- bam(fox ~ s(hour, bs = "cc", k = 8) + 
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) + 
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +                     
                        s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# cat model 1 
gam_cat_sp <- bam(cat ~ s(hour, bs = "cc", k = 8) + 
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) + 
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +                     
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# plot and save
plot_model_1()
# to combine figures, type in the mac terminal (using imagemagick): 
# convert figs/spatial_interaction_fox.png figs/spatial_interaction_cat.png  -append figs/spatial_interaction.png




# 2) Spatiotemporal activity across vegetation types and human footprint index ('model 2')----------------------

# fox model 2
gam_fox_veg <- bam(fox ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) + 
                     s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                     s(hfi, bs = "ts", k = 4) + 
                     ti(hfi, hour, bs = c("ts", "cc"), k = c(4, 8)) +
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# cat model 2
gam_cat_veg <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                     s(hfi, bs = "ts", k = 4) + 
                     ti(hfi, hour, bs = c("ts", "cc"), k = c(4, 8)) +
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

## plot vegetation effects models (separately)
# fox
png("figs/fox_veg.png", width = 9, height = 5, res = 600, units = "in")
plot_model_2(gam_fox_veg, records) + plot_annotation(title = "a.   Fox")
dev.off()
# cat
png("figs/cat_veg.png", width = 9, height = 5, res = 600, units = "in")
plot_model_2(gam_cat_veg, records) + plot_annotation(title = "b.   Feral cat")
dev.off()
# to combine figures, type in the terminal (using imagemagick): 
# convert figs/fox_veg.png figs/cat_veg.png  -append figs/predator_veg.png
# convert figs/predator_veg.png figs/gams_legend_veg.png  -append figs/predator_veg_leg.png

## plot and save fox hfi effects
plot_fox_hfi()




# Cat activity in response to fox counts ('model 3')----------------------

# cat model 3
gam_cat_fox <- bam(cat ~ habitat_type + 
                         t2(hour, fox_count_adj, by = habitat_type, bs = c("cc", "ts"), k = c(8, 5), full = TRUE) +  
                         s(station, bs = "re") +  
                         offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)


# fox model 3 (fox activity in each habitat type - for figure)
gam_fox_ht <- bam(fox ~ habitat_type + s(hour, by = habitat_type, bs = "cc", k = 8) + 
                        s(foxbaits, bs = "tp", k = 4) + 
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)


#  plot and save
png("figs/cat_fox_count.png", width = 9, height = 10, res = 600, units = "in")
plot_model_3()
dev.off()


## Plot the effect of foxbaits on foxes on each model
plot_foxbaits()


# END