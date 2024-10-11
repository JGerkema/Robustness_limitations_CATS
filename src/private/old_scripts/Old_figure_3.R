signif <- decomp_plot %>%
  filter(traits %in% c("Local")) %>%
  filter(pval.meta < 0.05) 

signif_try <- decomp_plot %>%
  filter(traits %in% c("Try")) %>%
  filter(pval.meta < 0.05) 

signif_temp <- decomp_plot %>%
  filter(traits %in% c("Local", "Try")) %>%
  filter(plot_study %in% signif$plot_study) %>%
  mutate(signif_local = "Yes") %>%
  select(KLR2, plot_study, signif_local, traits, study, Deviance, study_group) 

signif_both <- decomp_plot %>%
  filter(traits %in% c("Local", "Try")) %>%
  filter(plot_study %in% signif_try$plot_study) %>%
  mutate(signif_try = "Yes") %>%
  select(KLR2, plot_study, species_richness, signif_try, traits, study, Deviance, study_group) %>%
  full_join(signif_temp) %>%
  mutate_all(~ ifelse(is.na(.), "No", .)) %>%
  pivot_wider(names_from = "traits", values_from = "KLR2") %>%
  mutate(signif_cat = case_when(signif_local == "Yes" & signif_try == "No" ~ "Only with local data",
                                signif_local == "No" & signif_try == "Yes" ~ "Only with TRY data",
                                signif_local == "Yes" & signif_try == "Yes" ~ "With both")) %>%
  mutate(Deviance = factor(Deviance, levels = c("Dispersal mass effect",
                                                "Trait-based filtering",
                                                "Joint trait & dispersal mass effect",
                                                "Unexplained")))

plot_count_all <- decomp_plot %>%
  filter(traits %in% c("Local")) %>%
  select(all_of(c("Plot_code", "study_group"))) %>%
  distinct() %>%
  count(study_group) %>%
  rename(Dataset = study_group)


plots_signif_table <- signif_both %>%
  select(all_of(c("plot_study", "study_group", "signif_cat"))) %>%
  distinct() %>%
  group_by(study_group) %>%
  count(signif_cat) %>%
  ungroup() %>%
  pivot_wider(names_from = signif_cat, values_from = n) %>%
  mutate_all(~ ifelse(is.na(.), 0, .)) %>%
  relocate(study_group, "Only with local data", "Only with TRY data", "With both") %>%
  rename(Dataset = study_group) %>%
  full_join(plot_count_all) %>%
  select(all_of(c("Dataset", "Only with local data", "Only with TRY data", 
                  "With both", "n"))) %>%
  rename(`Analyzed plots (n)` = `n`)


write_xlsx(plots_signif_table, path = "docs/tableD1_plots_signif_table.xlsx")

model_types <- c("Dispersal mass effect", "Trait-based filtering",
                 "Joint trait & dispersal mass effect", "Unexplained")

cats_all_models_local <- cats_all_models_try <- model_summary_all <- data.frame()
for(i in 1:length(model_types) ){
  
  model_type_select <- model_types[i]
  
  cats_one_model_type <- signif_both %>%
    filter(Deviance %in% model_type_select) 
  
  
  plot_count2 <- cats_one_model_type %>%
    group_by(signif_cat) %>% 
    tally() %>%
    mutate(signif_cat = factor(signif_cat, levels = c("Only with local data", 
                                                      "Only with TRY data",
                                                      "With both"))) %>%
    mutate(plot_count_sign = paste0(signif_cat,  " (", n, ")"))
  
  cats_one_model_type <- cats_one_model_type %>%
    full_join(plot_count2)
  
  cats_one_model_type_local <- cats_one_model_type %>%
    filter(signif_cat %in% c("Only with local data", "With both")) 
  
  
  
  model_local <- lmerTest::lmer(Local ~ Try + (1 | study_group), 
                                data = cats_one_model_type_local)
  summary(model_local)
  
  model_summary <- summary(model_local)$coefficients
  model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                                 t_value_slope = round(model_summary[2,4], digits = 3),
                                 p_value_slope = round(model_summary[2,5], digits = 3),
                                 estimate_intercept = round(model_summary[1], digits = 3),
                                 t_value_intercept = round(model_summary[1,4], digits = 3),
                                 p_value_intercept = round(model_summary[1,5], digits = 3),
                                 Deviance = model_type_select)
  
  model_summary_all <- bind_rows(model_summary_all, model_summary_df) %>%
    relocate(Deviance)
  
  
  cats_one_model_type_local$F0 <- predictSE(model_local, cats_one_model_type_local, level = 0)$fit
  cats_one_model_type_local$SE <- predictSE(model_local, cats_one_model_type_local, level = 0)$se.fit
  
  cats_one_model_type_try <- signif_both %>%
    filter(signif_cat %in% c("Only with TRY data", "With both")) 
  
  model_try <- lmerTest::lmer(Local ~ Try + (1 | study_group), data = cats_one_model_type_try)
  summary(model_try)
  
  
  
  cats_one_model_type_try$F0 <- predictSE(model_try, cats_one_model_type_try, level = 0)$fit
  cats_one_model_type_try$SE <- predictSE(model_try, cats_one_model_type_try, level = 0)$se.fit
  
  cats_all_models_local <- bind_rows(cats_all_models_local, cats_one_model_type_local)
  cats_all_models_try <- bind_rows(cats_all_models_try, cats_one_model_type_try)
  
  plot_count <- signif_both %>%
    group_by(study_group) %>% 
    tally() %>%
    mutate(n = n/4)
  
}

local <- pull(plot_count2[2,3])
try <- pull(plot_count2[1,3])
both <- pull(plot_count2[3,3])


shape_colors <- c(local = "navy", 
                  try = "darkgreen",
                  both = "gold3")

signif_both_final <- signif_both %>%
  full_join(plot_count2)

cats_decomposed_local_try <-   ggplot(signif_both_final, aes(x = Try, y = Local)) +
  geom_point(alpha = 0.3, size = 3.5, aes(colour = study_group, shape = plot_count_sign)) +
  facet_wrap(~Deviance, scales = "free") +
  theme_classic() +
  #geom_smooth(method = "lm", se = T) + 
  labs(title = NULL,  
       y = "KLR2 with local data", x = "KLR2 with TRY data") +
  labs(col = "Dataset", shape = "Significance")+
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = 2) +
  guides(color = guide_legend(nrow = 6, order = 1), 
         shape = guide_legend(nrow = 3, order = 2, override.aes = list(color = shape_colors, alpha = 1))) + 
  # guides(shape = guide_legend(nrow = 3), order = 2) +
  scale_shape_manual(values = c("Only with local data (110)" = 16, 
                                "Only with TRY data (34)" = 17,  
                                "With both (47)" = 15), 
                     labels = c("Only with local data (110)", 
                                "Only with TRY data (34)",  
                                "With both (47)"),
                     breaks = c("Only with local data (110)", 
                                "Only with TRY data (34)",
                                "With both (47)")) +
  # geom_ribbon(aes( ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
  #        alpha = 0.2, linetype= 1, size=0.5, data = cats_all_models_try, fill = "darkseagreen3",
  #         colour = "darkseagreen3") +
  # geom_line(aes(y=(F0)), linewidth=1, colour = "darkgreen", data = cats_all_models_try) +
  geom_ribbon(aes( ymin = (F0 - 1.96 * SE), ymax = (F0 + 1.96 * SE)), 
              alpha = 0.2, linetype= 1, size=0.5, data = cats_all_models_local,
              fill = "slategray3", colour = "slategray3") +
  geom_line(aes(y=(F0)), linewidth=1, colour = "navy", data = cats_all_models_local) +
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size=15),
        legend.title = element_text(size=14),
        strip.text = element_text(size = 14),
        legend.position = "bottom")  +
  #coord_cartesian(ylim = c(0, axis_limit), xlim = c(0, axis_limit)) +
  scale_colour_discrete(labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) 
