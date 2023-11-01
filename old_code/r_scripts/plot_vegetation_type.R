

# function to run and plot GAMs for each species
gam_plot <- function(gam_model, data){
  ## Plot global smooth (using the gratia package)
  plot_global <- draw(evaluate_smooth(gam_model, overall_uncertainty = TRUE, "s(hour)"), xlab = "Hour",  title = "", subtitle = "i.   Average") + 
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") 
  ## Plot the difference smooths using the gratia package
  plot_difference <- draw(evaluate_smooth(gam_model, "s(hour,vegetation_group)"), xlab = "Hour",  title = "", subtitle = "ii.   Vegetation-specific deviations from the average") + 
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") 
  ## Plot the actual hourly activity for each vegetation type   
  # make a new dataframe to predict into
  df <- expand.grid(hour = 0:23,
                    vegetation_group = levels(data$vegetation_group),
                    foxbaits = 0, 
                    hfi = mean(data$hfi),
                    region = c("glenelg", "otways"),
                    survey_duration = 60)
  # get smooths to exlude
  x = sapply(gam_model$smooth, "[[",  "label")
  # predict model estimates (and uncertainty) into this dataframe, excluding the impact of variables other than hour:
  df <- cbind(df, predict.gam(gam_model, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[3:7])))
  plot_veg <- ggplot(data=df, aes(x=hour, y=fit, group=vegetation_group)) +
    # ylim(-13, -2) +
    geom_line(aes(y=fit, x=hour), lwd = 0.7) +
    facet_wrap(~vegetation_group, nrow =1) +
    geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
    geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
    labs(subtitle = "iii.   Predicted activity in each vegetation type", x = "Hour", y = "log(count)") 
  ## Assemble plots
  patch <- (plot_global | plot_difference) / plot_veg 
  patch
}

# save
png("figs/fox_veg.png", width = 9, height = 5, res = 600, units = "in")
gam_plot(gam_fox_veg, records) + plot_annotation(title = "a.   Fox")
dev.off()

png("figs/cat_veg.png", width = 9, height = 5, res = 600, units = "in")
gam_plot(gam_cat_veg, records) + plot_annotation(title = "b.   Feral cat")
dev.off()

# to combine figures, type in the terminal (using imagemagick): 
# convert figs/fox_veg.png figs/cat_veg.png  -append figs/predator_veg.png
# convert figs/predator_veg.png figs/gams_legend_veg.png  -append figs/predator_veg_leg.png

