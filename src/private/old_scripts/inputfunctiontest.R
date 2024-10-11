test_function3 <- function(input_veg, input_veg_fd, input_trait, input_trait_fd, trait_information, 
                           cutoff_value, trait_detail, study_name){
  
  extract <- new_names_synonyms %>%
    filter(study == study_name)
  
  not_included_traits <- c("Leaf_angle", "Photo_path", "Root_density", "Lateral_spread",
                           "Resprout", "Bark_thickness_rel", "Rhizome", "Stem_emerg",
                           "ShootDMC", "Leaf_shape", "Root_diam", "twFWC", "FWC")
  
  traits_orig <- trait_information[["traitlist"]]
  
  traits_filtered <- traits_orig[traits_orig %notin% not_included_traits]  
  
  traitlist <- sort(c(traits_filtered))
  cwm_names <- paste0(traitlist, "_cwm")
  cwv_names <- paste0(traitlist, "_cwv")
  
  trait_information <- lst(traitlist, cwm_names, cwv_names)
  
  if(nrow(extract) > 0){
    
    input_veg <- input_veg %>%
      rename(!!!setNames(extract$Orig_name, extract$Name)) %>%
      select(order(colnames(.)))
    
    input_trait <- input_trait %>%
      mutate(Name = case_when(
        Name %in% extract$Orig_name ~ extract$Name[match(Name, extract$Orig_name)],
        TRUE ~ Name)) %>%
      arrange(Name)
  } 
  
  
  if(trait_detail == "Study"){
    
    input_trait <- input_trait %>%
      select(Name, all_of(traits_filtered))

    input_trait_fd <- input_trait %>%
      filter(Name %in% trait_imp_specieslvl$Name) %>%
      column_to_rownames(var = "Name")
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
    
  } else{
    input_trait <- input_trait %>%
      select(Name, Plot_code, all_of(traits_filtered))
    
      input_trait_fd <- input_trait %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
  }
  
  
  
  prep_orig <- cats_preparation(input_veg = input_veg, 
                                input_trait = input_trait,
                                trait_information = trait_information, 
                                cutoff_value = 0.6, 
                                trait_detail = trait_detail)
  
 
  diversity_info_species <- calculate_diversity_info(veg_data = input_veg, 
                                                     input_veg_fd = input_veg_fd,
                                                     input_processed = prep_orig, 
                                                     trait_fd = input_trait_fd,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)

  
  input_trait_try <- trait_imp_specieslvl %>% 
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]])) 
  

  
  
  
  prep_try <- cats_preparation(input_veg = input_veg, 
                               input_trait = input_trait_try,
                               trait_information = trait_information, 
                               cutoff_value = 0.6, 
                               trait_detail = "Study")
  

  
  
  input_trait_try_fd <- input_trait_try %>%
    column_to_rownames(var = "Name") %>%
    filter(row.names(.) %in% colnames(prep_try[["relative_abundances"]][["matrix4traits_prop"]]))
  

  
  diversity_info_try_species <- calculate_diversity_info(veg_data = input_veg, 
                                                         input_veg_fd = input_veg_fd,
                                                         input_processed = prep_try, 
                                                         trait_fd = input_trait_try_fd,  
                                                         trait_detail =  "Study", 
                                                         study_name = study_name)
  
  

  if(!all.equal(colnames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]]), 
                colnames(prep_try[["relative_abundances"]][["matrix4traits_prop"]]))) {
    warning("Not the same species in original vs. try dataset")
  }
  
  if(!all.equal(rownames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]]), 
                rownames(prep_try[["relative_abundances"]][["matrix4traits_prop"]]))) {
    warning("Not the same plots in original vs. try dataset")
  }
  
  if(!all.equal(prep_orig[["relative_abundances"]][["matrix4traits_prop"]], 
                prep_try[["relative_abundances"]][["matrix4traits_prop"]])) {
    warning("Vegetation matrix is not the same")
  }
  
  if(foreach::getDoParRegistered() != TRUE){
    n_cores <- parallel::detectCores() - 5 
    
    my_cluster <- makeSOCKcluster(n_cores)
    
    registerDoSNOW(cl = my_cluster)}   
  
  cats_orig <- cats_calculations(prep_orig, trait_information, trait_detail, nperm = 30) %>%
    mutate(study = study_name)
  
  cats_try <- cats_calculations(prep_try, trait_information, trait_detail = "Study", nperm = 30) %>%
    mutate(study = study_name)
 
  
  all_ouput <- lst(diversity_info_species,  
                   diversity_info_try_species,  
                   prep_orig, prep_try,  
                   cats_orig,  cats_try) 
  
  return(all_ouput)
  
}







temp_fig_plot <- temp_raw_plot <- traits_n_plot <- data.frame()

for(i in 1:length(output_list)){
  
  
  cats_orig <- output_list[[i]][["cats_orig"]] %>%
    mutate(model.prior.plus.traits_uncorrected = model.prior.plus.traits,
           model.prior.plus.traits = 
             case_when(model.prior.plus.traits < model.mean.null.given.prior ~ model.mean.null.given.prior,
                       TRUE ~ model.prior.plus.traits))
  
  cats_try <- output_list[[i]][["cats_try"]] %>%
    mutate(model.prior.plus.traits_uncorrected = model.prior.plus.traits,
           model.prior.plus.traits = 
             case_when(model.prior.plus.traits < model.mean.null.given.prior ~ model.mean.null.given.prior,
                       TRUE ~ model.prior.plus.traits))
  
    diversity_info <- output_list[[i]][["diversity_info_species"]]
  study_name <- unique(diversity_info$study)
  
  cats_combined_orig <- cats_visualisation(cats_orig, "Orig")
  
  raw_orig <- cats_orig %>%
   select(study, Plot_code, model.bias, pval.uniform,
           model.uniform.prior.plus.traits,
           model.mean.null.given.prior,
           model.prior.plus.traits,
           model.prior.plus.traits_uncorrected,
           obs.stat_u) %>%
    rename(raw_meta = model.mean.null.given.prior,
           raw_trait = model.uniform.prior.plus.traits,
           raw_meta_trait = model.prior.plus.traits,
           raw_meta_trait_uncor = model.prior.plus.traits_uncorrected) %>%
    mutate(trait_minus_bias = raw_trait - model.bias) %>%
    relocate(pval.uniform,   obs.stat_u) %>%
    pivot_longer(cols = 5:10, values_to = "fit", names_to = "Model_type") %>%
    left_join(diversity_info) %>%
    mutate(traits = "Orig")

diversity_info <- output_list[[i]][["diversity_info_try_species"]]  

  raw_try <- cats_try %>%
   select(study, Plot_code, model.bias, pval.uniform,
           model.uniform.prior.plus.traits,
           model.mean.null.given.prior,
           model.prior.plus.traits,
           model.prior.plus.traits_uncorrected,
           obs.stat_u) %>%
    rename(raw_meta = model.mean.null.given.prior,
           raw_trait = model.uniform.prior.plus.traits,
           raw_meta_trait = model.prior.plus.traits,
           raw_meta_trait_uncor = model.prior.plus.traits_uncorrected) %>%
    mutate(trait_minus_bias = raw_trait - model.bias) %>%
    relocate(pval.uniform,   obs.stat_u) %>%
    pivot_longer(cols = 5:10, values_to = "fit", names_to = "Model_type") %>%
    left_join(diversity_info) %>%
    mutate(traits = "Try")
    

  cats_combined_try <- cats_visualisation(cats_try, "Try")
  
  
  temp_fig_plot <- bind_rows(temp_fig_plot, cats_combined_orig,
                              cats_combined_try
                        )
  
  temp_raw_plot <- bind_rows(temp_raw_plot, raw_orig, raw_try )
  
  dataset <- data_input_list[[i]]
  
  not_included_traits <- c("Leaf_angle", "Photo_path", "Root_density", "Lateral_spread",
                           "Resprout", "Bark_thickness_rel", "Rhizome", "Stem_emerg",
                           "ShootDMC", "Leaf_shape", "Root_diam", "twFWC", "FWC")
  
  traits_orig <- dataset[["trait_information"]][["traitlist"]]

  traits_filtered <- traits_orig[traits_orig %notin% not_included_traits]  
  
  traits <- data.frame(traits_n = length(traits_filtered), study = dataset[["study_name"]])
  
  traits_n_plot <- bind_rows(traits_n_plot, traits)
  
}

temp_raw_plot <- temp_raw_plot %>%
  mutate(plot_study = paste(Plot_code, study)) %>%
  filter(study %notin% c("Aiello1992", "Aiello2011" , "Asefa", "Osuri", "Frenette2009")) %>%
  mutate(study_group = case_when(study == "Liane_finland_1"  ~ "Liane",
                                 study == "Liane_finland_2"  ~ "Liane",
                                 study == "Liane_russia"  ~ "Liane",
                                 study == "Liane_sweden"  ~ "Liane",
                                 study == "Frenette2010" ~ "Frenette",
                                 study == "kober" ~ "Kober",
                                 TRUE ~ study))
library(lme4)
cats_species <- plot_cats4_function("Orig", "Try", "CATS at species level with original and TRY trait values",
                                    "KLR2 with original data", "KLR2 with TRY data",
                                    difference = "No")
set.seed(4)
cats_species

perm_values2 <- lst(perm_values)
all_KLR2 <- cbind(perm_values2,
  model.bias, model.uniform.prior.plus.traits,
  model.mean.null.given.prior, model.prior.plus.traits, model.mean.null.given.uniform,
  lambda_fit1, lambda_fit2, pval.uniform, pval.meta, 
  ci.pval_lower_u, ci.pval_upper_u, ci.pval_lower_m, ci.pval_upper_m,
  # mean_statdif_u, sd_statdif_u, obs.stat_u, ci.pval_lower_u, ci.pval_upper_u, 
  #  mean_statdif_m, sd_statdif_m, obs.stat_m, ci.pval_lower_m, ci.pval_upper_m
  obs.stat_u, obs.stat_m
) %>%
  as.data.frame() %>%
  mutate(Plot_code = plot)
