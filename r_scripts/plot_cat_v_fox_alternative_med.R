# PLOT FIGURE

## CAT PLOT - 3D 
df = expand.grid(hour = 0:23,
                 habitat_type = levels(records$habitat_type),
                 longitude = mean(records$longitude),
                 latitude = mean(records$latitude),
                 fox_count_adj = seq(min(records$fox_count_adj), max(records$fox_count_adj), length=150),
                 station = "T052",
                 survey_duration = 60)
# predict model results into dataframe
x = sapply(gam_cat_fox$smooth, "[[",  "label")
df_int <- cbind(df, predict(gam_cat_fox, newdata = df, se.fit = TRUE, type = "link", exclude = c(x[4:5], "s(survey_duration)")))

## CAT PLOT - 2D (same model at min / max fox)
# min df 
df_int_min <- df_int %>%
  group_by(habitat_type) %>%
  slice_min(fox_count_adj, n = 1) %>%
  mutate(minmax = "0min")
# median df
df_int_med <- df_int %>% 
  filter(abs(fox_count_adj - median(fox_count_adj)) == min(abs(fox_count_adj - median(fox_count_adj)))) %>% 
  mutate(minmax = "1mid")
# max df 
df_int_max <- df_int %>%
  group_by(habitat_type) %>%
  slice_max(fox_count_adj, n = 1) %>%
  mutate(minmax = "2max")
# join df
df_int_minmax <- bind_rows(df_int_min, df_int_med, df_int_max)
# rename facet
df_int_minmax$minmax <- as.factor(df_int_minmax$minmax)
df_int_minmax$minmax <- recode_factor(df_int_minmax$minmax,  
                                      "0min" = "Minimum",
                                      "1mid" = "Median",
                                      "2max" = "Maximum")
# plot
plot_wet_lowhigh <- ggplot(data=df_int_minmax, aes(x=hour, y=fit, group=minmax, col = minmax)) +
  geom_line(aes(y=fit, x=hour), lwd = 1.2) +
  facet_wrap(~minmax ~habitat_type, nrow = 3, scales = "free") +
  ggtitle("", subtitle = "Feral cat") +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit), fill = minmax), alpha=0.2) +
  scale_color_manual(values = c("dodgerblue", "purple", "tomato"), labels = c("Minimum", "Median", "Maximum"), name = "Fox activity") +   
  scale_fill_manual(values = c("dodgerblue", "purple", "tomato"), labels = c("Minimum", "Median", "Maximum"), name = "Fox activity") +    
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  ylab("log(count)") + 
  xlab("Hour") + 
  theme(legend.position = "bottom")

png("figs/cat_fox_min_med_max.png", width = 7, height = 10, res = 600, units = "in")
plot_fox / plot_wet_lowhigh + plot_annotation(tag_levels = "a") +
  plot_layout(heights = unit(c(1.5, 5), c('in', 'in')))
dev.off()

