## FIT GAMS TO MODEL DIEL ACTIVITY FOR FOXES AND FERAL CATS
# MATTHEW REES

library(dplyr)
library(mgcv)

# load data
records <- read.csv("derived_data/predator_counts_hour.csv")


# Transform variables -----------------------------------------------------

# new variable for fox counts adjusted for survey effort 
records$fox_count_adj <- records$fox_count / records$survey_duration
# cats more likely care about fox activity on the log scale
records$fox_count_adj <- log(records$fox_count_adj + 1)

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

# rename habitat type variable (for plotting)
records$habitat_type <- recode_factor(records$habitat_type,  "dry_glenelg" = "Dry vegetation (Glenelg region)",
                                      "dry_otway" = "Dry vegetation (Otway Ranges)",
                                      "wet_otway" = "Wet forest (Otway Ranges)")



# Changes in spatial activity across the daily cycle ----------------------

# specify penalties for 2D Duchon spline and 1D cc spline
m <- list(c(1,.5),rep(0,0)) 

gam_fox_sp <- bam(fox ~ t2(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m, full = TRUE) +
                        s(foxbaits, bs = "tp", k = 4) + 
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)


gam_cat_sp <- bam(cat ~ t2(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m, full = TRUE) +
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)




# Spatiotemporal activity across vegetation types ----------------------
gam_fox_veg <- bam(fox ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) + 
                     s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                     s(region, bs = "re") +  
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)


gam_cat_veg <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                     s(region, bs = "re") +  
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)


# Cat activity in response to fox counts ----------------------
gam_cat_fox <- bam(cat ~ habitat_type + t2(hour, fox_count_adj, by = habitat_type, bs = c("cc", "ts"), k = c(8, 5), full = TRUE) +  
                         s(station, bs = "re") +  
                         offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)




# Fox activity in each habitat type ----------------------
gam_fox_ht <- bam(fox ~ habitat_type + s(hour, by = habitat_type, bs = "cc", k = 8) + 
                        s(foxbaits, bs = "tp", k = 4) + 
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)



# Model summary table -----------------------------------------------------
extract_fits <- function(model, model_name){
  x <- summary(model)
  df <- expand.grid(species = as.character(x$formula)[2],
                    model = model_name,
                    edf = sum(x$edf),
                    dev.expl = x$dev.expl,
                    r.sq = x$r.sq)
return(df)
}

summaries <- bind_rows(extract_fits(gam_fox_sp, "1_spatial"),
                       extract_fits(gam_cat_sp, "1_spatial"),
                       extract_fits(gam_fox_veg, "2_vegetation_type"),
                       extract_fits(gam_cat_veg, "2_vegetation_type"),
                       extract_fits(gam_cat_fox, "3_fox_by_habitat_type"),
                       extract_fits(gam_fox_ht, "3_habitat_type"))  %>%
              arrange(species)
summaries
write.csv(summaries, "derived_data/model_summaries.csv")



# Plot set-up -------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(ggpubr)
library(viridis)
library(gratia)
library(forcats)

# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

# now run plot_x scripts...
