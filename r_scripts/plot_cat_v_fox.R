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
# plot
plot_range <- ggplot(aes(hour, fox_count_adj, fill = fit), data = df_int) +
  geom_tile() +
  facet_wrap(~habitat_type, nrow =1) +
  scale_fill_viridis("Cat activity", option = "viridis", breaks = c(-5, -3), labels = c("Low", "High")) + 
  ggtitle("", subtitle = "Feral cat") +
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  scale_y_continuous(breaks=c(min(records$fox_count_adj), max(records$fox_count_adj)), labels = c("Low", "High")) +
  ylab("Fox activity") + 
  xlab("Hour") 

## CAT PLOT - 2D (same model at min / max fox)
# min df 
df_int_min <- df_int %>%
  group_by(habitat_type) %>%
  slice_min(fox_count_adj, n = 1) %>%
  mutate(minmax = "0min")
# max df 
df_int_max <- df_int %>%
  group_by(habitat_type) %>%
  slice_max(fox_count_adj, n = 1) %>%
  mutate(minmax = "1max")
# join df
df_int_minmax <- rbind(df_int_min, df_int_max)
# plot
plot_wet_lowhigh <- ggplot(data=df_int_minmax, aes(x=hour, y=fit, group=minmax, col = minmax)) +
  geom_line(aes(y=fit, x=hour), lwd = 1.2) +
  facet_wrap(~habitat_type, nrow =1) +
  ggtitle("", subtitle = "Feral cat") +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit), fill = minmax), alpha=0.2) +
  scale_color_manual(values = c("dodgerblue", "tomato"), labels = c("Low", "High"), name = "Fox activity") +   
  scale_fill_manual(values = c("dodgerblue", "tomato"), labels = c("Low", "High"), name = "Fox activity") +    
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  ylab("log(count)") + 
  xlab("Hour")


## FOX PLOT - DIEL ACTIVITY X HABITAT TYPE
# 1B) fox diel x habitat type 
df = expand.grid(hour = 0:23,
                 habitat_type = levels(records$habitat_type),
                 station = "T052",
                 foxbaits = 0,
                 survey_duration = 60)
# predict model results into dataframe
x = sapply(gam_fox_ht$smooth, "[[",  "label")
df_int <- cbind(df, predict(gam_fox_ht, newdata = df, se.fit = TRUE, type = "link", exclude = c(x[4:8], "s(survey_duration)")))
# plot
plot_fox <- ggplot(aes(hour, fit), data = df_int) +
  geom_line(aes(y=fit, x=hour), lwd = 1.2) +
  facet_wrap(~habitat_type, nrow =1) +
  ggtitle("", subtitle = "Red Fox") +
  geom_ribbon(aes(ymin=(fit-2*se.fit), ymax=(fit+2*se.fit)), alpha=0.2) +
  geom_vline(xintercept = c(6.16,18.34), colour = "black", size = 0.6, linetype="dotted") + 
  ylab("log(count)") + 
  xlab("Hour")

#  save
png("figs/cat_fox_count.png", width = 9, height = 10, res = 600, units = "in")
plot_fox / plot_range / plot_wet_lowhigh + plot_annotation(tag_levels = "a")
dev.off()


# for slides
png("figs/cat_fox_count_pres1.png", width = 7.2, height = 5, res = 600, units = "in")
plot_fox / plot_spacer() + plot_annotation(tag_levels = "a")
dev.off()
  
png("figs/cat_fox_count_pres2.png", width = 7.2, height = 5.5, res = 600, units = "in")
plot_fox / plot_range + plot_annotation(tag_levels = "a")
dev.off()

png("figs/cat_fox_count_pres3.png", width = 7.2, height = 5.5, res = 600, units = "in")
plot_fox / plot_wet_lowhigh + plot_annotation(tag_levels = "a")
dev.off()

# END