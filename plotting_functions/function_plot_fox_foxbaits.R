plot_foxbaits <- function(){

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
   # filter(foxbaits == 0 | foxbaits == 19) %>%
    select(model, region, foxbaits, fit_resp, lower, upper) %>% 
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
    #filter(foxbaits == 0 | foxbaits == 19) %>%
    select(model, region, foxbaits, fit_resp, lower, upper) %>% 
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
    select(model, region, foxbaits, fit_resp, lower, upper) %>% 
    mutate(across(c(fit_resp, lower, upper), round, 3))
  
  # bind dataframes together and calculate percentage of suppression
  dfs <- bind_rows(newdf_1, newdf_2, newdf_3) %>%
    group_by(model, region) %>%
    mutate(suppression = 100 - min(fit_resp) / max(fit_resp) * 100) %>%
    arrange(region)
  
  # plots
  mod1 <- ggplot(filter(dfs, model == 1), aes(foxbaits, fit_resp, group = region)) + 
    facet_wrap(~region) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
    labs(subtitle = "Model 1", x = "Number of poison-baits within a 2.3 km radius", y = "Predicted counts")
  
  mod2 <- ggplot(filter(dfs, model == 2), aes(foxbaits, fit_resp, group = region)) + 
    facet_wrap(~region) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
    labs(subtitle = "Model 2", x = "Number of poison-baits within a 2.3 km radius", y = "Predicted counts")
  
  mod3 <- ggplot(filter(dfs, model == 3), aes(foxbaits, fit_resp, group = region)) + 
    facet_wrap(~region) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
    labs(subtitle = "Model 3", x = "Number of poison-baits within a 2.3 km radius", y = "Predicted counts")
  
  # plot together
  png("figs/foxbaiting_effects.png", width = 7, height = 10, res = 600, units = "in")
  print(mod1 / mod2 / mod3 + plot_annotation(tag_levels = "a"))
  dev.off()
  
}
