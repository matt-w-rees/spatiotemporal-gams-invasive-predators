plot_fox_hfi <- function(){

  # plot human footprint index  (no effect on diel activity - plot marginal effect only)
  plot <- smooth_estimates(gam_fox_veg, overall_uncertainty = TRUE, "s(hfi)") %>%
    mutate(lower = est + (1.96 * se),
           upper = est - (1.96 * se)) %>% 
    ggplot(aes(hfi, est)) + 
    geom_line() + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2)  + 
    labs(x = "Human footprint index", y = "Partial effect")
  
  # plot together
  png("figs/fox_hfi.png", width = 6, height = 4, res = 600, units = "in")
  print(plot)
  dev.off()
  
}
