plot_function_raw <- function(model_type, plot_title){

signif_both <- signif_both %>%
  filter(Model_type == model_type) 

plot_count2 <- signif_both %>%
  group_by(signif_cat) %>% 
  tally() %>%
  mutate(signif_cat = factor(signif_cat, levels = c("Only with original data", 
                                                    "Only with TRY data",
                                                    "With both"))) %>%
  mutate(plot_count_sign = paste0(signif_cat,  " (", n, ")"))

signif_both <- signif_both %>%
  full_join(plot_count2)

signif_both_orig <- signif_both %>%
  filter(signif_cat %in% c("Only with original data", "With both")) 

model_orig <- lmerTest::lmer(Try ~ Orig + (1 | study), data = signif_both_orig)
summary(model_orig)
# Calculate marginal and conditional R-squared
#r_squared <- r.squaredGLMM(model)

signif_both_orig$F0 <- predictSE(model_orig, signif_both_orig, level = 0)$fit
signif_both_orig$SE <- predictSE(model_orig, signif_both_orig, level = 0)$se.fit

signif_both_try <- signif_both %>%
  filter(signif_cat %in% c("Only with TRY data", "With both")) 

model_try <- lmerTest::lmer(Try ~ Orig + (1 | study), data = signif_both_try)
summary(model_try)
# Calculate marginal and conditional R-squared
#r_squared <- r.squaredGLMM(model)

signif_both_try$F0 <- predictSE(model_try, signif_both_try, level = 0)$fit
signif_both_try$SE <- predictSE(model_try, signif_both_try, level = 0)$se.fit



plot_count <- signif_both %>%
  group_by(study_group) %>% 
  tally()

orig <- pull(plot_count2[2,3])
try <- pull(plot_count2[1,3])
both <- pull(plot_count2[3,3])


shape_colors <- c(orig = "navy", 
                  try = "darkgreen",
                  both = "gold3")


plot <- ggplot(signif_both, aes(x = Orig, y = Try)) +
  geom_point(alpha = 0.3, size = 5, aes(colour = study_group, shape = plot_count_sign)) +
  scale_shape_manual(values = c("Only with original data (120)" = 16, 
                                "Only with TRY data (62)" = 17,  
                                "With both (108)" = 15), 
                     labels = c("Only with original data (120)", 
                                "Only with TRY data (62)",  
                                "With both (108)"),
                     breaks = c("Only with original data (120)", 
                                "Only with TRY data (62)",
                                "With both (108)"),
                     guide = guide_legend(nrow = 3, override.aes = list(color = shape_colors, alpha = 1))) +
  theme_classic() +
  #geom_smooth(method = "lm", se = T) + 
  labs(title = NULL,  
       y = "KLR2 with TRY data", x = "KLR2 with original data") +
  labs(col = "Dataset", shape = "Significance")+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = 2) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 6)) + 
  geom_ribbon(aes( ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
              alpha = 0.2, linetype= 1, size=0.5, data = signif_both_try, fill = "darkseagreen3",
              colour = "darkseagreen3") +
  geom_line(aes(y=(F0)), size=1, colour = "darkgreen", data = signif_both_try) +
  geom_ribbon(aes( ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
              alpha = 0.2, linetype= 1, size=0.5, data = signif_both_orig,
              fill = "slategray3", colour = "slategray3") +
  geom_line(aes(y=(F0)), size=1, colour = "navy", data = signif_both_orig) +
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size=15),
        legend.title = element_text(size=14))  +
  coord_cartesian(ylim = c(0, 0.8), xlim = c(0, 0.8)) +
    scale_colour_discrete(labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) # +
  #ggplot2::annotate("text", x = 0, y = 1, 
 #                   label = paste("Estimated slope (plots significant with original trait values - blue):", 
  #                                round(fixef(model_orig)["Orig"], 3)), color = "black", hjust = 0) +
 # ggplot2::annotate("text", x = 0, y = 0.97, 
 #                   label = paste("Standard Error:", round(sqrt(diag(vcov(model_orig))["Orig"]), 3)), color = "black", hjust = 0) +
 # ggplot2::annotate("text", x = 0, y = 0.94, 
 #                   label = paste("t-value:", round(summary(model_orig)$coefficients["Orig", "t value"], 3)), color = "black", hjust = 0) +
 # ggplot2::annotate("text", x = 0, y = 0.88, 
 #                   label = paste("Estimated slope (plots significant with TRY trait values - green):", 
 #                                 round(fixef(model_try)["Orig"], 3)), color = "black", hjust = 0) +
 # ggplot2::annotate("text", x = 0, y = 0.85, 
 #                   label = paste("Standard Error:", round(sqrt(diag(vcov(model_try))["Orig"]), 3)), color = "black", hjust = 0) +
 # ggplot2::annotate("text", x = 0, y = 0.82, 
  #                  label = paste("t-value:", round(summary(model_try)$coefficients["Orig", "t value"], 3)), color = "black", hjust = 0)   

return(plot)

}


plot_cats4_function <- function(y_var, x_var, plot_title,
                                x_lab, y_lab, figure_table){
  cats <- raw_plot %>%
    filter(traits == y_var) %>%
    filter(Model_type %notin% c("trait_minus_bias", "raw_meta_trait_uncor")) %>%
    select(-c(FEve, FRic, qual.FRic, RaoQ, traits,
              obs.stat_u)) %>%
    distinct() %>%
    select(fit, pval.meta, FDiv, study_group, plot_study, Model_type) 
  
  signif <- raw_plot %>%
    filter(traits %in% c(y_var)) %>%
    filter(pval.meta < 0.05) 
  
  cats_both <- raw_plot %>%
    filter(traits == x_var) %>%
    filter(Model_type %notin% c("trait_minus_bias", "raw_meta_trait_uncor")) %>%
    select(-c(FEve, FRic, qual.FRic, RaoQ, traits, 
              obs.stat_u)) %>%
    select(fit, pval.meta, FDiv, study_group, plot_study, Model_type) %>%
    rename(fit_try = fit, pval.meta_try = pval.meta, FDiv_try = FDiv) %>%
    distinct() %>%
    full_join(cats) %>%
   # drop_na() %>%
    mutate(dif_fit = fit - fit_try, 
           dif_FDiv = FDiv - FDiv_try) %>%
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
                   "Metacommunity and traits (m,t)","Traits (u,t)")
  
  cats_both_models <- model_summary_all <- data.frame()

    for(i in 1:length(model_types) ){
    
    model_type_select <- model_types[i]
    
    cats_one_model_type <- cats_both %>%
      filter(Model_type %in% model_type_select) %>%
      drop_na(fit)

    model <- lmerTest::lmer(fit ~ fit_try + (1 | study_group), data = cats_one_model_type)
    
    model_summary <- summary(model)$coefficients
    
    model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                                   t_value_slope = round(model_summary[2,4], digits = 3),
                                   p_value_slope = model_summary[2,5],
                                   estimate_intercept = round(model_summary[1], digits = 3),
                                   t_value_intercept = round(model_summary[1,4], digits = 3),
                                   p_value_intercept = model_summary[1,5],
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

plot <-  ggplot(cats_both_models, aes(x = fit_try, y = fit)) +
    geom_point(alpha = 0.2, aes(colour = study_group), size = 4) +
    theme_classic() +
    #geom_smooth(method = "lm", se = F) + 
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




plot_dif_function_zoom <- function(variable, y_lab, x_lab, df) {
  
  
  div_df <- df %>%
    rename("var_int" = {{variable}}) %>%
    drop_na(var_int)
  
  model <- lmerTest::lmer(Fit_div ~ var_int + (1|study_group), data = div_df)
  summary(model)

  div_df$F0 <- predictSE(model, div_df, level = 0)$fit
  div_df$SE <- predictSE(model, div_df, level = 0)$se.fit
  
  plot_count <- div_df %>%
    group_by(study_group) %>% 
    tally()
  
  plot <- ggplot(div_df, aes(y = Fit_div, x = var_int)) +
    geom_point(alpha = 0.3, size = 5, aes(colour = study_group)) +
    theme_classic() +
    labs(title = NULL,  
         y = y_lab, x = x_lab) +
    labs(col = "Dataset")+
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3)) + 
    coord_cartesian(ylim = c(-0.25, 0.75)) +
    geom_ribbon(aes(ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
                alpha = 0.2, linetype = 1, size = 0.5, data = div_df,
                fill = "slategray3", colour = "slategray3") +
    geom_line(aes(y = F0), size = 1, colour = "navy", data = div_df) +
    theme(axis.text = element_text(size = 15),
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 15),
          legend.title = element_text(size = 16))  +
    scale_colour_discrete(labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) 
  
  return(plot)
}

plot_dif_function <- function(variable, plot_title, x_lab, df) {
  
  
  div_df <- df %>%
    rename("var_int" = {{variable}}) %>%
    drop_na(var_int)
  
  model <- lmer(Fit_div ~ var_int + (1|study_group), data = div_df)
  summary(model)
  
  div_df$F0 <- predictSE(model, div_df, level = 0)$fit
  div_df$SE <- predictSE(model, div_df, level = 0)$se.fit
  
  plot_count <- div_df %>%
    group_by(study_group) %>% 
    tally()
  
  plot <- ggplot(div_df, aes(y = Fit_div, x = var_int)) +
    geom_point(alpha = 0.3, size = 3.5, aes(colour = study_group)) +
    theme_classic() +
    labs(title = NULL,  
         y = "Difference in KLR2 (Original - TRY)", x = x_lab) +
    labs(col = "Dataset")+
    theme(legend.position = "bottom") +
    guides(color = guide_legend(nrow = 3)) + 
    geom_ribbon(aes(ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
                alpha = 0.2, linetype = 1, size = 0.5, data = div_df,
                fill = "slategray3", colour = "slategray3") +
    geom_line(aes(y = F0), size = 1, colour = "navy", data = div_df) +
    theme(axis.text = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 13),
          plot.title = element_text(size = 15))  +
    scale_colour_discrete(labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) +
    ggplot2::annotate("text", x = -0.4, y = 0.8, 
                      label = paste("Estimated slope:", 
                                    round(fixef(model)[["var_int"]], 3)), color = "black", hjust = 0) +
    ggplot2::annotate("text", x = -0.4, y = 0.77, 
                      label = paste("Standard Error:", round(sqrt(diag(vcov(model)))[["var_int"]], 3)), color = "black", hjust = 0) +
    ggplot2::annotate("text", x = -0.4, y = 0.74, 
                      label = paste("t-value:", round(summary(model)$coefficients[["var_int", "t value"]], 3)), color = "black", hjust = 0) 
  
  return(plot)
}

plot_boxplot <- function(model_type, plot_title){

cats_compare <- raw_plot %>%
  filter(traits %in% c("Try",  "Orig")) %>%
  filter(plot_study %in% signif$plot_study) %>%
  select(c(study, Plot_code, Model_type, fit, traits, study_group)) %>%
  filter(Model_type == model_type) %>%
  mutate(traits = recode(traits, Orig = "Original", Try = "TRY database")) %>%
  mutate(traits = recode(traits, Orig_genus = "Original", Try_genus = "TRY database")) 


plot_count <- cats_compare %>% 
  group_by(study_group) %>% 
  tally() %>%
  mutate(n = n/2)

stat.test <- cats_compare %>%
  full_join(plot_count) %>%
  mutate(study_group_n = paste0(study_group, " (", n, ")")) %>%
  group_by(study_group_n) %>%
  filter(n() >2) %>%
  pairwise_t_test(
    fit ~ traits, paired = TRUE, detailed = T) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  select(-df,  -method, -alternative, ) 

extra <- data.frame(study_group_n = "Osuri (1)",
                    estimate = (0.410408486 - 0.466798108))


levels_vector <- stat.test %>%
  bind_rows(extra) %>%
  arrange(-estimate) %>%
  pull(study_group_n) %>%
  unique()

cats_compare <- cats_compare %>%
  full_join(plot_count) %>%
  mutate(study_group_n = paste0(study_group, " (", n, ")")) %>%
  mutate(study_group_n = factor(study_group_n, levels = levels_vector, ordered = TRUE)) #%>%
  #filter(study_group_n %in% stat.test$study_group_n) 


plot <- ggplot(cats_compare, aes(x = study_group_n, y = fit)) +
  geom_boxplot(aes(fill = traits)) +
  labs(title = NULL,
       x = NULL,
       y = "KLR2",
       fill = "Source trait values") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "bottom") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", x = "study_group_n",
                     step.increase = 0.08, hide.ns = TRUE, tip.length = 2,
                     y.position = 1.1) +
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 13))

return(plot)

}

