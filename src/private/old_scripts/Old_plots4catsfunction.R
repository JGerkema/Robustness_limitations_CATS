plot_cats4_function <- function(y_var, x_var, plot_title,
                                x_lab, y_lab, figure_table){
  cats <- raw_plot %>%
    filter(traits == y_var) %>%
    filter(Model_type %notin% c("trait_minus_bias", "raw_meta_trait_uncor")) %>%
    select(-c(FEve, FRic, qual.FRic, RaoQ, traits,
              obs.stat_u)) %>%
    distinct() %>%
    select(fit, pval.meta, study_group, plot_study, Model_type) 
  
  signif <- raw_plot %>%
    filter(traits %in% c(y_var)) %>%
    filter(pval.meta < 0.05) 
  
  cats_both <- raw_plot %>%
    filter(traits == x_var) %>%
    filter(Model_type %notin% c("trait_minus_bias", "raw_meta_trait_uncor")) %>%
    select(-c(FEve, FRic, FDiv, qual.FRic, RaoQ, traits, 
              obs.stat_u)) %>%
    select(fit, pval.meta, study_group, plot_study, Model_type) %>%
    rename(fit_xvar = fit, pval.meta_xvar = pval.meta) %>%
    distinct() %>%
    full_join(cats) %>%
    mutate(Model_type = recode(Model_type, model.bias = "Estimate of model bias (u)", raw_meta = "Metacommunity (m)",
                               raw_meta_trait = "Metacommunity and traits (m,t)", raw_trait = "Traits (u,t)",
                               trait_minus_bias = "Explanatory power of traits beyond model bias"))  %>%
    mutate(Model_type = factor(Model_type, levels = c("Estimate of model bias (u)", 
                                                      "Metacommunity (m)",
                                                      "Traits (u,t)",
                                                      "Metacommunity and traits (m,t)")))
  
  
  cats_both_signif <- cats_both  %>%
    filter(plot_study %in% signif$plot_study)
  
  
  model_types <- c("Estimate of model bias (u)", "Metacommunity (m)",
                   "Traits (u,t)",  "Metacommunity and traits (m,t)")
  
  cats_both_models <- model_summary_all <- data.frame()
  
  for(i in 1:length(model_types) ){
    
    model_type_select <- model_types[i]
    
    cats_one_model_type <- cats_both %>%
      filter(Model_type %in% model_type_select) %>%
      drop_na(fit)
    
    model <- lmerTest::lmer(fit ~ fit_xvar + (1 | study_group), data = cats_one_model_type)
    
    model_summary <- summary(model)$coefficients
    
    model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                                   t_value_slope = round(model_summary[2,4], digits = 3),
                                   p_value_slope = round(model_summary[2,5], digits = 3),
                                   estimate_intercept = round(model_summary[1], digits = 3),
                                   t_value_intercept = round(model_summary[1,4], digits = 3),
                                   p_value_intercept = round(model_summary[1,5], digits = 3),
                                   model_type = model_type_select)
    
    
    cats_one_model_type$F0 <- predictSE(model, cats_one_model_type, level = 0)$fit
    cats_one_model_type$SE <- predictSE(model, cats_one_model_type, level = 0)$se.fit
    
    cats_both_models <- bind_rows(cats_both_models, cats_one_model_type)
    model_summary_all <- bind_rows(model_summary_all, model_summary_df)
  }
  
  plot_count <- cats_both_models %>% 
    group_by(study_group) %>% 
    tally() %>%
    mutate(n = n/4)
  
  plot <-  ggplot(cats_both_models, aes(x = fit_xvar, y = fit)) +
    geom_point(alpha = 0.2, aes(colour = study_group), size = 4) +
    theme_classic() +
    labs(title = NULL,  
         y = y_lab, x = x_lab) +
    labs(col = "Dataset")+
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = 2) +
    facet_wrap(~Model_type) +
    theme(legend.position = "bottom") +
    geom_ribbon(aes( ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
                alpha = 0.2, linetype= 1, linewidth=0.5, data = cats_both_models, fill = "navy",
                colour = "navy") +
    geom_line(aes(y=(F0)), size=1, colour = "navy", data = cats_both_models) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
    guides(color = guide_legend(nrow = 4)) + 
    theme(axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          plot.title = element_text(size=15),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    scale_colour_discrete(labels = paste0(plot_count$study_group, " (", plot_count$n, ")"))
  
  
  if(figure_table == "figure"){
    return(plot)}else{
      return(model_summary_all)}
}