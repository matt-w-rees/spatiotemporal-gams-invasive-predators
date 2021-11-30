# PLOT SPATIO-TEMPORAL MODELS
# 1) DIEL ACTIVITY STRENGTH ACROSS SPACE
# 2) ANIMATED GIF FOR SPATIOTEMPORAL ACTIVITY


gam_fox_sp <- bam(fox ~ t2(x, y, hour, d = c(2, 1), bs = c("tp", "cc"), k = c(80, 8), full = TRUE) +
                    s(foxbaits, bs = "tp", k = 4) + 
                    s(station, bs = "re") +  
                    offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)

gam_cat_sp <- bam(cat ~ t2(x, y, hour, d = c(2, 1), bs = c("tp", "cc"), k = c(80, 8),full = TRUE) +
                    s(station, bs = "re") +  
                    offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)




gam_fox_sp <- bam(fox ~ ti(x, y, bs = "tp", k = 80) + 
                        ti(hour, bs = "cc", k = 8) +
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "tp"), k = c(80, 8)) +
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = FALSE)



# additional packages
library(sp)
library(sf)
library(viridis)

# split records
records_g <- filter(records, region == "glenelg")
records_o <- filter(records, region == "otways")

# cam locations 
camdata <- records %>%
  distinct(., station, .keep_all = TRUE) %>%
  select(region, station, x, y) 
# add new variable to match variable we will facet by - region x species
camdata_fox <- camdata
camdata_fox$species <- "Red fox"
camdata_cat <- camdata
camdata_cat$species <- "Feral cat"
camdata <- rbind(camdata_fox, camdata_cat)
camdata$facet <- paste0(camdata$species, "_", camdata$region)
unique(camdata$facet)

# MAKE SPATIAL DATAFRAMES TO PREDICT INTO ----------------------------------------------

## GLENELG:
# take just the cams
records_g_cams <- distinct(records_g, station, .keep_all = TRUE)
# change it to sf class
records_g_cams <- st_as_sf(records_g_cams, coords = c("x", "y"), crs = 32754) 
# make a 4km buffer around each camera
glenelg_cams_buffer = st_buffer(records_g_cams, 4000)
# dissolve the buffer
glenelg_cams_buffer = st_union(glenelg_cams_buffer)
# convert back to dataframe
glenelg_buffer_df <- fortify(as_Spatial(glenelg_cams_buffer))%>%
  transmute(x = long, y = lat, order = order, group = group)
# split by each group (seperate polygon for each grid)
buffer_df1 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.1"),]
buffer_df2 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.2"),]
buffer_df3 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.3"),]
buffer_df4 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.4"),]
# make a dataframe from this to predict into
data_g_plot = expand.grid(
  x = seq(min(glenelg_buffer_df$x), 
          max(glenelg_buffer_df$x),
          length=50),
  y = seq(min(glenelg_buffer_df$y),
          max(glenelg_buffer_df$y),
          length=50),
  survey_duration = 60,
  foxbaits = 0, 
  hour = c(0:23),
  station = "A001")
# subset this data to just locations within buffers
data_g_plot1 = data_g_plot[with(data_g_plot, inSide(buffer_df1, x, y)),]
data_g_plot2 = data_g_plot[with(data_g_plot, inSide(buffer_df2, x, y)),]
data_g_plot3 = data_g_plot[with(data_g_plot, inSide(buffer_df3, x, y)),]
data_g_plot4 = data_g_plot[with(data_g_plot, inSide(buffer_df4, x, y)),]
data_g_plot <- rbind(data_g_plot1, data_g_plot2, data_g_plot3, data_g_plot4)


## OTWAYS
#take just the cams
records_o_cams <- distinct(records_o, station, .keep_all = TRUE)
# change it to sf class
records_o_cams <- st_as_sf(records_o_cams, coords = c("x", "y"), crs = 32754) 
# make a 4km buffer around each camera
otways_cams_buffer = st_buffer(records_o_cams, 6500)
# dissolve the buffer
otways_cams_buffer = st_union(otways_cams_buffer)
# convert back to dataframe
otways_buffer_df <- fortify(as_Spatial(otways_cams_buffer))%>%
  transmute(x = long, y = lat, order = order, group = group)
## make dataframe to predict intro from this
data_o_plot = expand.grid(
  x = seq(min(otways_buffer_df$x), 
          max(otways_buffer_df$x),
          length=50),
  y = seq(min(otways_buffer_df$y),
          max(otways_buffer_df$y),
          length=50),
  station = "T053", 
  year = 2018,
  survey_duration = 60,
  hour = 0:23,
  foxbaits = 0)
# subset df just locations within buffer zone
data_o_plot = data_o_plot[with(data_o_plot, inSide(otways_buffer_df, x, y)),]



# PREDICT MODELS ----------------------------------------------------------
# predict model estimates into dataframes and combine

## GLENELG
# fox
x = sapply(gam_fox_sp$smooth, "[[",  "label")
data_g_plot_fox <- cbind(data_g_plot, predict.gam(gam_fox_sp, newdata = data_g_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)", "s(survey_duration)")))
data_g_plot_fox$species <- "fox"
# cat
data_g_plot_cat <- cbind(data_g_plot, predict.gam(gam_cat_sp, newdata = data_g_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)", "s(survey_duration)")))
data_g_plot_cat$species <- "cat"
# bind dataframes
pred_glenelg <- bind_rows(data_g_plot_fox, data_g_plot_cat)
pred_glenelg$year <- NA
head(pred_glenelg)

## OTWAYS
# fox
data_o_plot_fox <- cbind(data_o_plot, predict.gam(gam_fox_sp, newdata = data_o_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2],"s(station)", "s(survey_duration)")))
data_o_plot_fox$species <- "fox"
# cat
data_o_plot_cat <- cbind(data_o_plot, predict.gam(gam_cat_sp, newdata = data_o_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2],"s(station)", "s(survey_duration)")))
data_o_plot_cat$species <- "cat"
# bind dataframes
pred_otways <- bind_rows(data_o_plot_fox, data_o_plot_cat)
head(pred_otways)

# now bind regions
pred_glenelg$region <- "glenelg"
pred_otways$region <- "otways"
model_predictions <- bind_rows(pred_glenelg, pred_otways)
# new variable for faceting
model_predictions$facet <- paste0(model_predictions$species, "_", model_predictions$region)
head(model_predictions)




# ADJUST DATAFRAME --------------------------------------------------------

# rename species
model_predictions$species <- as.factor(model_predictions$species)
model_predictions$species <- recode_factor(model_predictions$species,  
                                           "cat" = "Feral cat",
                                           "fox" = "Red fox")


# filter everything to region - currently there is a bug with geom_tile which fucks up faceting by species x region UGH
# https://github.com/tidyverse/ggplot2/issues/849
model_predictions_g <- filter(model_predictions, region == "glenelg")
model_predictions_o <- filter(model_predictions, region == "otways")
camdata_g <- filter(camdata, region == "glenelg")
camdata_o <- filter(camdata, region == "otways")


# FACET - HOUR ------------------------------------------------------------
# non-gif version for pdf's

## GLENELG
plot_facet_g_fox <- ggplot(aes(x, y, fill = fit), data = filter(model_predictions_g, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("Activity", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Red fox") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  



plot_facet_g_cat <- ggplot(aes(x, y, fill = fit), data = filter(model_predictions_g, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("Activity", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Feral cat") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  



## OTWAYS
plot_facet_o_fox <- ggplot(aes(x, y, fill = fit), data = filter(model_predictions_o, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("Activity", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Red fox") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  



plot_facet_o_cat <- ggplot(aes(x, y, fill = fit), data = filter(model_predictions_o, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("Activity", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Feral cat") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  


# save 
png("figs/spte_facet_g_fox.png", width = 7.5, height = 10, res = 600, units = "in")
plot_facet_g_fox
dev.off()

png("figs/spte_facet_g_cat.png", width = 7.5, height = 10, res = 600, units = "in")
plot_facet_g_cat
dev.off()

png("figs/spte_facet_o_fox.png", width = 8, height = 10, res = 600, units = "in")
plot_facet_o_fox
dev.off()

png("figs/spte_facet_o_cat.png", width = 8, height = 10, res = 600, units = "in")
plot_facet_o_cat
dev.off()


# END