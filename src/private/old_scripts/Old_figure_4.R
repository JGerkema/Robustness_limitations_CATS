cats_decomposed_specieslvl <- decomp_plot %>%
  filter(traits == "Local") %>%
  mutate(Deviance = factor(Deviance, levels = c("Trait-based filtering", 
                                                "Dispersal mass effect",
                                                "Joint trait & dispersal mass effect",
                                                "Unexplained"))) %>%
  filter(plot_study %in% signif$plot_study) %>%
  ungroup() %>%
  mutate(species_richness = as.numeric(species_richness))

plot_count <- cats_decomposed_specieslvl %>%
  group_by(study_group) %>% 
  tally() %>%
  mutate(n = n/4)

model_types <- unique(cats_decomposed_specieslvl$Deviance)

cats_all_models <- model_summary_all <- data.frame()

## Model 1
model_type_select <- model_types[1]

cats_one_model_type <- cats_decomposed_specieslvl %>%
  filter(Deviance %in% model_type_select) %>%
  mutate(log_KLR2 = log(KLR2)) %>%
  filter(!is.infinite(log_KLR2))


model1 <- lmerTest::lmer(log_KLR2 ~ species_richness + (1 | study_group), 
                         data = cats_one_model_type)

model2 <- lmerTest::lmer(log_KLR2 ~ species_richness + (species_richness | study_group),
                         data = cats_one_model_type)

anova(model1, model2) #AIC model1: 190.24 and model2: 192.18  

model_summary <- summary(model1)$coefficients
model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                               t_value_slope = round(model_summary[2,4], digits = 3),
                               p_value_slope = round(model_summary[2,5], digits = 3),
                               estimate_intercept = round(model_summary[1], digits = 3),
                               t_value_intercept = round(model_summary[1,4], digits = 3),
                               p_value_intercept = round(model_summary[1,5], digits = 3),
                               Deviance = model_type_select)

model_summary_all <- bind_rows(model_summary_all, model_summary_df) %>%
  relocate(Deviance)

cats_one_model_type$F0 <- (predictSE(model1, cats_one_model_type, level = 0)$fit)
cats_one_model_type$SE <- (predictSE(model1, cats_one_model_type, level = 0)$se.fit)
cats_one_model_type$F1 <- fitted(model1, level = 1) 


cats_one_model_type <- cats_one_model_type %>%
  mutate(ymin = exp(F0 - 1.96 * SE),
         ymax = exp(F0 + 1.96 * SE),
         y=exp(F0),
         F1=exp(F1))

cats_all_models <- bind_rows(cats_all_models, cats_one_model_type)


# Model 2
model_type_select <- model_types[2]

cats_one_model_type <- cats_decomposed_specieslvl %>%
  filter(Deviance %in% model_type_select)

model1 <- lmerTest::lmer(KLR2 ~ species_richness + (1 | study_group), 
                         data = cats_one_model_type)

model2 <- lmerTest::lmer(KLR2 ~ species_richness + (species_richness | study_group), 
                         data = cats_one_model_type)

anova(model1, model2) #AIC model1: -85.579 and model2: -82.055

model_summary <- summary(model1)$coefficients
model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                               t_value_slope = round(model_summary[2,4], digits = 3),
                               p_value_slope = round(model_summary[2,5], digits = 3),
                               estimate_intercept = round(model_summary[1], digits = 3),
                               t_value_intercept = round(model_summary[1,4], digits = 3),
                               p_value_intercept = round(model_summary[1,5], digits = 3),
                               Deviance = model_type_select)

model_summary_all <- bind_rows(model_summary_all, model_summary_df) %>%
  relocate(Deviance)

cats_one_model_type$F0 <- predictSE(model1, cats_one_model_type, level = 0)$fit
cats_one_model_type$SE <- predictSE(model1, cats_one_model_type, level = 0)$se.fit
cats_one_model_type$F1 <- fitted(model1, level = 1) 

cats_one_model_type <- cats_one_model_type %>%
  mutate(ymin = (F0 - 1.96 * SE),
         ymax = (F0 + 1.96 * SE),
         y=(F0))

cats_all_models <- bind_rows(cats_all_models, cats_one_model_type)


# Model 3
model_type_select <- model_types[3]

cats_one_model_type <- cats_decomposed_specieslvl %>%
  filter(Deviance %in% model_type_select)


model1 <- lmerTest::lmer(KLR2 ~ species_richness + (1 | study_group), 
                         data = cats_one_model_type)

model2 <- lmerTest::lmer(KLR2 ~ species_richness + (species_richness | study_group), 
                         data = cats_one_model_type)

anova(model1, model2) #AIC model1: -61.269 and model2: -60.772 

model_summary <- summary(model1)$coefficients
model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                               t_value_slope = round(model_summary[2,4], digits = 3),
                               p_value_slope = round(model_summary[2,5], digits = 3),
                               estimate_intercept = round(model_summary[1], digits = 3),
                               t_value_intercept = round(model_summary[1,4], digits = 3),
                               p_value_intercept = round(model_summary[1,5], digits = 3),
                               Deviance = model_type_select)

model_summary_all <- bind_rows(model_summary_all, model_summary_df) %>%
  relocate(Deviance)

cats_one_model_type$F0 <- predictSE(model1, cats_one_model_type, level = 0)$fit
cats_one_model_type$SE <- predictSE(model1, cats_one_model_type, level = 0)$se.fit
cats_one_model_type$F1 <- fitted(model1, level = 1) 

cats_one_model_type <- cats_one_model_type %>%
  mutate(ymin = (F0 - 1.96 * SE),
         ymax = (F0 + 1.96 * SE),
         y=(F0))

cats_all_models <- bind_rows(cats_all_models, cats_one_model_type)

## Model 4
model_type_select <- model_types[4]

cats_one_model_type <- cats_decomposed_specieslvl %>%
  filter(Deviance %in% model_type_select) 


model1 <- lme4::lmer(exp(KLR2) ~ species_richness + (1 | study_group), 
                     data = cats_one_model_type)

model2 <- lmerTest::lmer(exp(KLR2) ~ species_richness + 
                           (species_richness | study_group), 
                         data = cats_one_model_type)
fixef(model1)
ranef(model1)

anova(model1, model2) #AIC model1: 78.334 and model2: 63.521

model_summary <- summary(model2)$coefficients
model_summary_df <- data.frame(estimate_slope = round(model_summary[2], digits = 3),
                               t_value_slope = round(model_summary[2,4], digits = 3),
                               p_value_slope = round(model_summary[2,5], digits = 3),
                               estimate_intercept = round(model_summary[1], digits = 3),
                               t_value_intercept = round(model_summary[1,4], digits = 3),
                               p_value_intercept = round(model_summary[1,5], digits = 3),
                               Deviance = model_type_select)

model_summary_all <- bind_rows(model_summary_all, model_summary_df) %>%
  relocate(Deviance)

cats_one_model_type$F0 <- (predictSE(model2, cats_one_model_type, level = 0)$fit)
cats_one_model_type$SE <- (predictSE(model2, cats_one_model_type, level = 0)$se.fit)
cats_one_model_type$F1 <- fitted(model2, level = 1) 

cats_one_model_type <- cats_one_model_type %>%
  mutate(ymin = log(F0 - 1.96 * SE),
         ymax = log(F0 + 1.96 * SE),
         y=log(F0),
         F1=log(F1))

cats_all_models <- bind_rows(cats_all_models, cats_one_model_type)

cats_all_models <- cats_all_models %>%
  mutate(Deviance = factor(Deviance, levels = c("Dispersal mass effect",
                                                "Trait-based filtering", 
                                                "Joint trait & dispersal mass effect",
                                                "Unexplained")))


cats_decomposed_speciesrichness <- ggplot(cats_all_models, aes(y = KLR2, x = species_richness)) +
  geom_point(alpha = 0.1, size = 3.5, aes(colour = study_group)) +
  facet_wrap(~Deviance, scales = "free") +
  theme_classic() +
  labs(title = NULL,  
       y = "KLR2 with local data", x = "Species richness") +
  labs(col = "Dataset", shape = "Significance") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 3)) + 
  geom_line(aes(y = F1, colour = study_group), alpha = 0.8, size = 1)+
  geom_ribbon(aes(ymin = ymin, ymax = ymax), 
              alpha = 0.2, linetype = 1, linewidth = 0.5, data = cats_all_models,
              fill = "gray80", colour = "gray80") +
  geom_line(aes(y = y), linewidth = 1, colour = "gray20", data = cats_all_models) +
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 15),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  scale_colour_viridis_d(option = "H", 
                         labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) 


write_xlsx(model_summary_all, path = "docs/tableE3_cats_decomposed_speciesrichness_table.xlsx")  

jpeg("Figures/Figure4_cats_decomposed_speciesrichness.jpeg",  width = 3150, height = 2300, units = "px", res = 300)
plot(cats_decomposed_speciesrichness) #Figure 4.
dev.off()