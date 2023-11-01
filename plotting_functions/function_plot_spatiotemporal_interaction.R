plot_model_1 <- function(){
  
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
  
  
  # MAKE SPATIAL DATAFRAMES TO PREDICT INTO ----------------------------------------------
  
  ## GLENELG:
  # take just the cams
  records_g_cams <- distinct(records_g, station, .keep_all = TRUE)
  # change it to sf class
  records_g_cams <- st_as_sf(records_g_cams, coords = c("x", "y"), crs = 3111) 
  # make a 4km buffer around each camera
  glenelg_cams_buffer = st_buffer(records_g_cams, 4000)
  # dissolve the buffer
  glenelg_cams_buffer = st_union(glenelg_cams_buffer)
  # but bring in the edges - don't want to predict outside of surveyed area
  glenelg_cams_buffer <- st_buffer(glenelg_cams_buffer, dist = -3500, joinStyle  = "ROUND")
  # convert back to dataframe
  glenelg_buffer_df <- fortify(as_Spatial(glenelg_cams_buffer))%>%
    transmute(x = long, y = lat, order = order, group = group)
  # split by each group (seperate polygon for each grid)
  buffer_df1 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.1"),]
  buffer_df2 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.2"),]
  buffer_df3 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.3"),]
  buffer_df4 <- glenelg_buffer_df[which(glenelg_buffer_df$group == "ID1.4"),]
  # make this a list for x and y cols
  buffer_df1 <- list(buffer_df1$x, buffer_df1$y)
  buffer_df2 <- list(buffer_df2$x, buffer_df2$y)
  buffer_df3 <- list(buffer_df3$x, buffer_df3$y)
  buffer_df4 <- list(buffer_df4$x, buffer_df4$y)
  # rename lists
  names(buffer_df1) <- c("x", "y")
  names(buffer_df2) <- c("x", "y")
  names(buffer_df3) <- c("x", "y")
  names(buffer_df4) <- c("x", "y")
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
    hour = c(0:23))
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
  records_o_cams <- st_as_sf(records_o_cams, coords = c("x", "y"), crs = 3111) 
  # make a 4km buffer around each camera
  otways_cams_buffer = st_buffer(records_o_cams, 6500)
  # dissolve the buffer
  otways_cams_buffer = st_union(otways_cams_buffer)
  # but bring in the edges - don't want to predict outside of surveyed area
  otways_cams_buffer <- st_buffer(otways_cams_buffer, dist = -6000, joinStyle  = "ROUND")
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
    year = 2018,
    survey_duration = 60,
    hour = 0:23,
    foxbaits = 0)
  # make this a list for x and y cols
  otways_buffer_df <- list(otways_buffer_df$x, otways_buffer_df$y)
  # rename lists
  names(otways_buffer_df) <- c("x", "y")
  # subset df just locations within buffer zone
  data_o_plot = data_o_plot[with(data_o_plot, inSide(otways_buffer_df, x, y)),]
  
  
  
  # PREDICT MODELS ----------------------------------------------------------
  # predict model estimates into dataframes and combine
  
  ## GLENELG
  # fox
  data_g_plot_fox <- cbind(data_g_plot, predict.gam(gam_fox_sp, newdata = data_g_plot, se.fit = FALSE, newdata.guaranteed = TRUE, type = "terms", terms = "ti(hour,x,y)"))
  data_g_plot_fox$species <- "fox"
  # cat
  data_g_plot_cat <- cbind(data_g_plot, predict.gam(gam_cat_sp, newdata = data_g_plot, se.fit = FALSE, newdata.guaranteed = TRUE, type = "terms", terms = "ti(hour,x,y)"))
  data_g_plot_cat$species <- "cat"
  # bind dataframes
  pred_glenelg <- bind_rows(data_g_plot_fox, data_g_plot_cat)
  
  ## OTWAYS
  # fox
  data_o_plot_fox <- cbind(data_o_plot, predict.gam(gam_fox_sp, newdata = data_o_plot, se.fit = FALSE, newdata.guaranteed = TRUE, type = "terms", terms = "ti(hour,x,y)"))
  data_o_plot_fox$species <- "fox"
  # cat
  data_o_plot_cat <- cbind(data_o_plot, predict.gam(gam_cat_sp, newdata = data_o_plot, se.fit = FALSE, newdata.guaranteed = TRUE, type = "terms", terms = "ti(hour,x,y)"))
  data_o_plot_cat$species <- "cat"
  # bind dataframes
  pred_otways <- bind_rows(data_o_plot_fox, data_o_plot_cat)
  
  # now bind regions
  pred_glenelg$region <- "glenelg"
  pred_otways$region <- "otways"
  model_predictions <- bind_rows(pred_glenelg, pred_otways)
  
  
  # ADJUST DATAFRAME --------------------------------------------------------
  # rename species
  model_predictions$species <- as.factor(model_predictions$species)
  model_predictions$species <- recode_factor(model_predictions$species,  
                                             "fox" = "Red fox",
                                             "cat" = "Feral cat")
  model_predictions <- rename(model_predictions, 'fit' = 'ti(hour,x,y)')
  
  # filter everything to region - currently there is a bug with geom_tile which fucks up faceting by species x region UGH
  # https://github.com/tidyverse/ggplot2/issues/849
  model_predictions_g <- filter(model_predictions, region == "glenelg")
  model_predictions_o <- filter(model_predictions, region == "otways")
  camdata_g <- filter(camdata, region == "glenelg")
  camdata_o <- filter(camdata, region == "otways")
  
  
  
  # SELECT FACET ------------------------------------------------------------
  # easier comparison of species at a select few hours: 0,6,12,18
  cat_g_select <- filter(model_predictions_g, species == "Feral cat" & (hour == 0 | hour == 6 | hour == 12 | hour == 18))
  fox_g_select <- filter(model_predictions_g, species == "Red fox" & (hour == 0 | hour == 6 | hour == 12 | hour == 18))
  cat_o_select <- filter(model_predictions_o, species == "Feral cat" & (hour == 0 | hour == 6 | hour == 12 | hour == 18))
  fox_o_select <- filter(model_predictions_o, species == "Red fox" & (hour == 0 | hour == 6 | hour == 12 | hour == 18))
  
  # get min max values so plots will be on the same scale (used in viridis function)
  fox_minval <- min(c(fox_g_select$fit, fox_o_select$fit))
  fox_maxval <- max(c(fox_g_select$fit, fox_o_select$fit))
  cat_minval <- min(c(cat_g_select$fit, cat_o_select$fit))
  cat_maxval <- max(c(cat_g_select$fit, cat_o_select$fit))
  
  # build plots
  plot_facet_select_g_fox <- ggplot(aes(x, y, fill = fit), data = fox_g_select) +
    geom_tile() +
    scale_fill_distiller("Partial effect", palette = "RdBu", direction = -1, limits = c(fox_minval, fox_maxval)) + 
   # geom_contour(aes(z = fit, fill = NULL), colour = "black") +
    geom_point(data = camdata_g, fill = NA, col = "black", size = 0.02, alpha = 0.05, shape = 3) +
    ggtitle("", subtitle = "ii)   Spatial deviation from the marginal diel activity: Glenelg region") + 
    facet_wrap(~hour, nrow = 1) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill='grey85'),
          legend.position = "none") 
  
  plot_facet_select_o_fox <- ggplot(aes(x, y, fill = fit), data = fox_o_select) +
    geom_tile() +
    scale_fill_distiller("Partial effect", palette = "RdBu", direction = -1, limits = c(fox_minval, fox_maxval)) + 
   # geom_contour(aes(z = fit, fill = NULL), colour = "black") +
    geom_point(data = camdata_o, fill = NA, col = "black", size = 0.02, alpha = 0.05, shape = 3) +
    ggtitle("", subtitle = "iii)   Spatial deviation from the marginal diel activity: Otway Ranges") + 
    facet_wrap(~hour, nrow = 1) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill='grey85'),
          legend.position = "bottom")  
  
  plot_facet_select_g_cat <- ggplot(aes(x, y, fill = fit), data = cat_g_select) +
    geom_tile() +
    scale_fill_distiller("Partial effect", palette = "RdBu", direction = -1, limits = c(cat_minval, cat_maxval)) + 
    #geom_contour(aes(z = fit, fill = NULL), colour = "black") +
    geom_point(data = camdata_g, fill = NA, col = "black", size = 0.02, alpha = 0.05, shape = 3) +
    ggtitle("", subtitle = "ii)   Spatial deviation from the marginal diel activity: Glenelg region") + 
    facet_wrap(~hour, nrow = 1) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill='grey85'),
          legend.position = "none")  
  
  plot_facet_select_o_cat <- ggplot(aes(x, y, fill = fit), data = cat_o_select) +
    geom_tile() +
    scale_fill_distiller("Partial effect", palette = "RdBu", direction = -1, limits = c(cat_minval, cat_maxval)) + 
    #geom_contour(aes(z = fit, fill = NULL), colour = "black") +
    geom_point(data = camdata_o, fill = NA, col = "black", size = 0.02, alpha = 0.05, shape = 3) +
    ggtitle("", subtitle = "iii)   Spatial deviation from the marginal diel activity: Otway Ranges") + 
    facet_wrap(~hour, nrow = 1) + 
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill='grey85'),
          legend.position = "bottom")  
  
  
  
  ## Plot global smooth 
  fox_hour <- smooth_estimates(gam_fox_sp, overall_uncertainty = TRUE, "s(hour)") %>%
    mutate(lower = est + (1.96 * se),
           upper = est - (1.96 * se)) %>% 
    ggplot(aes(hour, est)) + 
     geom_line() + 
     geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
     geom_vline(xintercept = c(6.16,18.34), colour = "black", linewidth = 0.6, linetype="dotted") + 
     labs(subtitle = "i)   Marginal diel activity", x = "Hour", y = "Partial effect")
    
  cat_hour <- smooth_estimates(gam_cat_sp, overall_uncertainty = TRUE, "s(hour)") %>%
    mutate(lower = est + (1.96 * se),
           upper = est - (1.96 * se)) %>% 
    ggplot(aes(hour, est)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
    geom_vline(xintercept = c(6.16,18.34), colour = "black", linewidth = 0.6, linetype="dotted") + 
    labs(subtitle = "i)   Marginal diel activity", x = "Hour", y = "Partial effect")
  
  
  
  png("figs/spatial_interaction_fox.png", width = 8, height = 7, res = 600, units = "in")
  print(fox_hour / plot_facet_select_g_fox / plot_facet_select_o_fox + plot_annotation(title = "(a)   Fox"))
  dev.off()
  
  png("figs/spatial_interaction_cat.png", width = 8, height = 7, res = 600, units = "in")
  print(cat_hour / plot_facet_select_g_cat / plot_facet_select_o_cat + plot_annotation(title = "(b)   Feral cat"))
  dev.off()

} 