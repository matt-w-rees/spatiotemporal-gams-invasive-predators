# PLOT SEPERATED MARGINAL AND INTERACTION EFFECTS FOR SPATIAL MODELS (1)


# RUN MODEL WITH TI TENSOR PRODUCT ----------------------
# so we can separate marginal effects from interaction

# specify penalties for 2D Duchon spline and 1D cc spline
m <- list(c(1,.5),rep(0,0)) 

gam_fox_sp <- bam(fox ~ s(hour, bs = "cc", k = 8) +
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) +
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +
                        s(foxbaits, bs = "tp", k = 4) + 
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)


gam_cat_sp <- bam(cat ~ s(hour, bs = "cc", k = 8) +
                        s(x, y, bs = "ds", k = 80, m = c(1,.5)) +
                        ti(x, y, hour, d = c(2, 1), bs = c("ds", "cc"), k = c(80, 8), m = m) +
                        s(station, bs = "re") +  
                        offset(log(survey_duration)), 
                  data = records, family = nb, knots = list(hour = c(0, 23)), nthreads = 3, discrete = TRUE)





# PLOT AVERAGE HOUR -------------------------------------------------------

# new df 
df <- expand.grid(x = 625380,
                  y = 5760730,
                  survey_duration = 60,
                  foxbaits = 0, 
                  hour = 0:23,
                  station = "A001")


# predict fox
x = sapply(gam_fox_sp$smooth, "[[",  "label")
df_hr_fox <- cbind(df, predict.gam(gam_fox_sp, newdata = df, se.fit = TRUE, type = "link", exclude = c(x[2:4], x[5])))
df_hr_fox$species <- "fox"

# predict cat
x = sapply(gam_cat_sp$smooth, "[[",  "label")
df_hr_cat <- cbind(df, predict.gam(gam_cat_sp, newdata = df, se.fit = TRUE, type = "link", exclude = c(x[2:4])))
df_hr_cat$species <- "cat"

# merge dataframes
df_hr <- bind_rows(df_hr_fox, df_hr_cat)

# rename species
df_hr$species <- as.factor(df_hr$species)
df_hr$species <- recode_factor(df_hr$species,  
                                           "fox" = "Red fox",
                                           "cat" = "Feral cat",
)

# plot
plot_hr <- ggplot(data=df_hr, aes(x=hour, y=fit)) +
  geom_line(aes(y=fit, x=hour), lwd = 1) +
  facet_wrap(~species, nrow =1, scales = "free") +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  labs(title = "Average diel pattern", subtitle = "Marginal effect (model 1)", x = "Hour", y = "log(count)") 
plot_hr


png("figs/avg_diel_predator.png", width = 10, height = 5, res = 600, units = "in")
plot_hr
dev.off()



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


# MAKE SPATIAL DATAFRAMES AND PREDICT INTO ----------------------------------------------

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



# predict model estimates into dataframes and combine

## GLENELG
# fox
x = sapply(gam_fox_sp$smooth, "[[",  "label")
# predict interaction
data_g_plot_fox <- cbind(data_g_plot, predict.gam(gam_fox_sp, newdata = data_g_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)")))
# predict marginal - space
data_g_plot_fox <- cbind(data_g_plot_fox, predict.gam(gam_fox_sp, newdata = data_g_plot_fox, se.fit = TRUE, type = "link", exclude = c(x[1], x[3], "s(station)")))
# rename fits
names(data_g_plot_fox) <- c(names(data_g_plot_fox)[1:6], "fit_int", "se.fit_int", "fit_mar", "se.fit_mar")
# add species col
data_g_plot_fox$species <- "fox"

# cat
x = sapply(gam_cat_sp$smooth, "[[",  "label")
# predict interaction
data_g_plot_cat <- cbind(data_g_plot, predict.gam(gam_cat_sp, newdata = data_g_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)")))
# predict marginal - space
data_g_plot_cat <- cbind(data_g_plot_cat, predict.gam(gam_cat_sp, newdata = data_g_plot_cat, se.fit = TRUE, type = "link", exclude = c(x[1], x[3], "s(station)")))
# rename fits
names(data_g_plot_cat) <- c(names(data_g_plot_cat)[1:6], "fit_int", "se.fit_int", "fit_mar", "se.fit_mar")
# add species col
data_g_plot_cat$species <- "cat"

# bind dataframes
pred_glenelg <- bind_rows(data_g_plot_fox, data_g_plot_cat)
pred_glenelg$year <- NA
head(pred_glenelg)

## OTWAYS
# fox
x = sapply(gam_fox_sp$smooth, "[[",  "label")
# predict interaction
data_o_plot_fox <- cbind(data_o_plot, predict.gam(gam_fox_sp, newdata = data_o_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)")))
# predict marginal - space
data_o_plot_fox <- cbind(data_o_plot_fox, predict.gam(gam_fox_sp, newdata = data_o_plot_fox, se.fit = TRUE, type = "link", exclude = c(x[1], x[3], "s(station)")))
# rename fits
names(data_o_plot_fox) <- c(names(data_o_plot_fox)[1:7], "fit_int", "se.fit_int", "fit_mar", "se.fit_mar")
# add species col
data_o_plot_fox$species <- "fox"

# cat
x = sapply(gam_cat_sp$smooth, "[[",  "label")
# predict interaction
data_o_plot_cat <- cbind(data_o_plot, predict.gam(gam_cat_sp, newdata = data_o_plot, se.fit = TRUE, type = "link", exclude = c(x[1:2], "s(station)")))
# predict marginal - space
data_o_plot_cat <- cbind(data_o_plot_cat, predict.gam(gam_cat_sp, newdata = data_o_plot_cat, se.fit = TRUE, type = "link", exclude = c(x[1], x[3], "s(station)")))
# rename fits
names(data_o_plot_cat) <- c(names(data_o_plot_cat)[1:7], "fit_int", "se.fit_int", "fit_mar", "se.fit_mar")
# add species col
data_o_plot_cat$species <- "cat"

# bind sp dataframes
pred_otways <- bind_rows(data_o_plot_fox, data_o_plot_cat)
head(pred_otways)

# now bind regions
pred_glenelg$region <- "glenelg"
pred_otways$region <- "otways"
model_predictions <- bind_rows(pred_glenelg, pred_otways)
# new variable for faceting
model_predictions$facet <- paste0(model_predictions$species, "_", model_predictions$region)
head(model_predictions)

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


# PLOT INTERACTION --------------------------------------------------------

## GLENELG
plot_facet_g_fox <- ggplot(aes(x, y, fill = fit_int), data = filter(model_predictions_g, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Red fox", subtitle = "Space-time interaction effect (excluding marginal effects)") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  

plot_facet_g_cat <- ggplot(aes(x, y, fill = fit_int), data = filter(model_predictions_g, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Feral cat", subtitle = "Space-time interaction effect (excluding marginal effects)") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  

## OTWAYS
plot_facet_o_fox <- ggplot(aes(x, y, fill = fit_int), data = filter(model_predictions_o, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Red fox", subtitle = "Space-time interaction effect (excluding marginal effects)") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  

plot_facet_o_cat <- ggplot(aes(x, y, fill = fit_int), data = filter(model_predictions_o, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.05, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Feral cat", subtitle = "Space-time interaction effect (excluding marginal effects)") + 
  facet_wrap(~hour) + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  


# save 
png("figs/spte_diff_avg_g_fox.png", width = 7.5, height = 10, res = 600, units = "in")
plot_facet_g_fox
dev.off()

png("figs/spte_diff_avg_g_cat.png", width = 7.5, height = 10, res = 600, units = "in")
plot_facet_g_cat
dev.off()

png("figs/spte_diff_avg_o_fox.png", width = 8, height = 10, res = 600, units = "in")
plot_facet_o_fox
dev.off()

png("figs/spte_diff_avg_o_cat.png", width = 8, height = 10, res = 600, units = "in")
plot_facet_o_cat
dev.off()


# PLOT MARGINAL SPACE -----------------------------------------------------

## GLENELG
plot_g_fox <- ggplot(aes(x, y, fill = fit_mar), data = filter(model_predictions_g, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.5, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Red fox", subtitle = "Space marginal effect") + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  

plot_g_cat <- ggplot(aes(x, y, fill = fit_mar), data = filter(model_predictions_g, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_g, fill = NA, col = "white", size = 0.02, alpha = 0.5, shape = 3) +
  theme_bw(10) + 
  ggtitle("Glenelg region: Feral cat", subtitle = "Space marginal effect") + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  


## OTWAYS
plot_o_fox <- ggplot(aes(x, y, fill = fit_mar), data = filter(model_predictions_o, species == "Red fox")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.5, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Red fox", subtitle = "Space marginal effect") + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  

plot_o_cat <- ggplot(aes(x, y, fill = fit_mar), data = filter(model_predictions_o, species == "Feral cat")) +
  geom_tile() +
  scale_fill_viridis("log(count)", option = "viridis") +
  geom_point(data = camdata_o, fill = NA, col = "white", size = 0.02, alpha = 0.5, shape = 3) +
  theme_bw(10) + 
  ggtitle("Otway Ranges: Feral cat", subtitle = "Space marginal effect") + 
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom")  


# save 
png("figs/sp_marginal_o.png", width = 10, height = 5,res = 600, units = "in")
plot_o_fox + plot_o_cat 
dev.off()

png("figs/sp_marginal_g.png", width = 10, height = 5, res = 600, units = "in")
plot_g_fox + plot_g_cat
dev.off()


# END