## FIT GAMS TO MODEL DIEL ACTIVITY FOR FOXES AND FERAL CATS
# MATTHEW REES

library(dplyr)
library(mgcv)
# for plots
library(gratia)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(viridis)
library(forcats)
# Set the default theme for ggplot objects 
theme_set(theme_bw())
theme_update(panel.grid = element_blank())


# load data
records <- read.csv("derived_data/predator_counts_hour.csv")


# Transform variables -----------------------------------------------------
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
                      survey_duration = as.integer(survey_duration),
                      hfi = round(hfi, digits = 2))




# 1) Interaction between spatial and diel activity ('model 1') ---------------------
# specify penalties for 2D Duchon spline and 1D cc spline
m <- list(c(1,.5),rep(0,0)) 

gam_fox_sp <- bam(fox ~ s(hour, bs = "cc", k = 8) + 
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) + 
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +                     
                        s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

gam_cat_sp <- bam(cat ~ s(hour, bs = "cc", k = 8) + 
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) + 
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +                     
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# plot model 1 with "plot_spatiotemporal_interaction.R script. 


# 2) Spatiotemporal activity across vegetation types and human footprint index ('model 2')----------------------
gam_fox_veg <- bam(fox ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) + 
                     s(foxbaits, region, bs = "fs", xt = list(bs = "tp"), k = 4) + 
                     s(hfi, bs = "ts", k = 4) + 
                     ti(hfi, hour, bs = c("ts", "cc"), k = c(4, 8)) +
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)


gam_cat_veg <- bam(cat ~ s(hour, bs = "cc", k = 8) +   
                     s(hour, vegetation_group, bs = "fs", xt = list(bs = "cc"), k = 8) +  
                     s(hfi, bs = "ts", k = 4) + 
                     ti(hfi, hour, bs = c("ts", "cc"), k = c(4, 8)) +
                     s(station, bs = "re") +  
                     offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# plot model 2 with "plot_vegetation_type.R script (plots across vegetation types) and "plot_fox_baits_hfi.R" (plots of hfi)


# Cat activity in response to fox counts ('model 3')----------------------
gam_cat_fox <- bam(cat ~ habitat_type + 
                         t2(hour, fox_count_adj, by = habitat_type, bs = c("cc", "ts"), k = c(8, 5), full = TRUE) +  
                         s(station, bs = "re") +  
                         offset(log(survey_duration)), 
                   data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# Fox activity in each habitat type ----------------------
gam_fox_ht <- bam(fox ~ habitat_type + s(hour, by = habitat_type, bs = "cc", k = 8) + 
                        s(foxbaits, bs = "tp", k = 4) + 
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), discrete = TRUE)

# plot model 3 with "cat_v_fox.R script. 


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
write.csv(summaries, "derived_data/model_summaries.csv")


# now run plot_x scripts...