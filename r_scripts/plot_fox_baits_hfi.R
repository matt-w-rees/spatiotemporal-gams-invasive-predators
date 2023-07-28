# plot spatial effects of human footprint and fox control 


# HUMAN FOOTPRINT INDEX PLOT  ---------------------------------------------
# plot human footprint index  (no effect on diel activity - plot marginal effect only)
fox_hfi <- smooth_estimates(gam_fox_veg, overall_uncertainty = TRUE, "s(hfi)") %>%
  mutate(lower = est + (1.96 * se),
         upper = est - (1.96 * se)) %>% 
  ggplot(aes(hfi, est)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
  labs(x = "Human footprint index", y = "Partial effect")

# plot together
png("figs/fox_hfi.png", width = 6, height = 4, res = 600, units = "in")
fox_hfi
dev.off()


# ESTIMATES OF FOX SUPPRESSION ----------------------------------------------------------
# MODEL 1 
x = sapply(gam_fox_sp$smooth, "[[",  "label")
df <- expand.grid(hour = 0,
                  x = mean(records$x),
                  y = mean(records$y),
                  foxbaits = min(records$foxbaits):max(records$foxbaits), 
                  hfi = mean(records$hfi),
                  region = c("glenelg", "otways"),
                  survey_duration = 60)
newdf_1 <- cbind(df, predict.gam(gam_fox_sp, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[c(1,2,3)], "s(survey_duration)"))) %>%
  mutate(model = "1",
         fit_resp = exp(fit), 
         lower = exp(fit - (1.96 * se.fit)),
         upper = exp(fit + (1.96 * se.fit))) %>%
  filter(foxbaits == 0 | foxbaits == 19) %>%
  select(model, region, fit_resp, lower, upper) %>% 
  mutate(across(c(fit_resp, lower, upper), round, 3))

# MODEL 2
x = sapply(gam_fox_veg$smooth, "[[",  "label")
df <- expand.grid(hour = 0,
                     vegetation_group = "Dry Forest",
                     foxbaits = min(records$foxbaits):max(records$foxbaits), 
                     hfi = mean(records$hfi),
                     region = c("glenelg", "otways"),
                     survey_duration = 60)

newdf_2 <- cbind(df, predict.gam(gam_fox_veg, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[c(1,2,4,5,6, 7)], "s(survey_duration)"))) %>%
  mutate(model = "2",
         fit_resp = exp(fit), 
         lower = exp(fit - (1.96 * se.fit)),
         upper = exp(fit + (1.96 * se.fit))) %>%
  filter(foxbaits == 0 | foxbaits == 19) %>%
  select(model, region, fit_resp, lower, upper) %>% 
  mutate(across(c(fit_resp, lower, upper), round, 3))
         

# MODEL 3
x = sapply(gam_fox_ht$smooth, "[[",  "label")
df <- expand.grid(hour = 0,
                  habitat_type = "Dry vegetation (Otway Ranges)",
                  region = c("glenelg", "otways"),
                  foxbaits = min(records$foxbaits):max(records$foxbaits), 
                  survey_duration = 60)
newdf_3 <- cbind(df, predict.gam(gam_fox_ht, newdata = df, se.fit = TRUE, type = "link", newdata.guaranteed = TRUE, exclude = c(x[c(1,2,3,5)], "s(survey_duration)"))) %>%
  mutate(model = "3",
         fit_resp = exp(fit), 
         lower = exp(fit - (1.96 * se.fit)),
         upper = exp(fit + (1.96 * se.fit))) %>%
  filter(foxbaits == 0 | foxbaits == 19) %>%
  select(model, region, fit_resp, lower, upper) %>% 
  mutate(across(c(fit_resp, lower, upper), round, 3))

# bind dataframes together and calcuate factor of suppression
dfs <- bind_rows(newdf_1, newdf_2, newdf_3) %>%
  group_by(model, region) %>%
  mutate(suppression = max(fit_resp) / min(fit_resp)) %>%
  arrange(region)




# PARTIAL EFFECT PLOTS OF FOX SUPPRESSION ----------------------------------------------------------

# fox control plot - model 1 
fox_baits_a <- smooth_estimates(gam_fox_sp, overall_uncertainty = TRUE, "s(foxbaits,region)") %>%
  mutate(lower = est + (1.96 * se),
         upper = est - (1.96 * se),
         region = if_else(region == "glenelg", "Glenelg region", "Otway Ranges")) %>% 
  ggplot(aes(foxbaits, est, group = region)) + 
  facet_wrap(~region) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
  labs(subtitle = "Model 1", x = "Poison fox-bait density within 2.3 km radius", y = "Partial effect")

# fox control plot - model 2
fox_baits_b <- smooth_estimates(gam_fox_veg, overall_uncertainty = TRUE, "s(foxbaits,region)") %>%
  mutate(lower = est + (1.96 * se),
         upper = est - (1.96 * se),
         region = if_else(region == "glenelg", "Glenelg region", "Otway Ranges")) %>% 
  ggplot(aes(foxbaits, est, group = region)) + 
  facet_wrap(~region) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
  labs(subtitle = "Model 2", x = "Poison fox-bait density within 2.3 km radius", y = "Partial effect")


# fox control plot - model 3
fox_baits_c <- smooth_estimates(gam_fox_ht, overall_uncertainty = TRUE, "s(foxbaits,region)") %>%
  mutate(lower = est + (1.96 * se),
         upper = est - (1.96 * se),
         region = if_else(region == "glenelg", "Glenelg region", "Otway Ranges")) %>% 
  ggplot(aes(foxbaits, est, group = region)) + 
  facet_wrap(~region) +
  geom_line() + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
  labs(subtitle = "Model 3", x = "Poison fox-bait density within 2.3 km radius", y = "Partial effect")

fox_baits_a / fox_baits_b / fox_baits_c + plot_annotation(tag_levels = "a")
 

# plot together
png("figs/foxbaiting_effects.png", width = 7, height = 10, res = 600, units = "in")
fox_baits_a / fox_baits_b / fox_baits_c + plot_annotation(tag_levels = "a")
dev.off()

# END
