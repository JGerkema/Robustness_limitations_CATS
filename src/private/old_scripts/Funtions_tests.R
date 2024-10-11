test_function <- function(input_veg, input_veg_fd, input_trait, input_trait_fd, trait_information, 
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
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup()
    
    input_trait_genus_fd <- input_trait_genus %>%
      column_to_rownames(var = "Name")
    
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F)))) 
    
  } else{
    input_trait <- input_trait %>%
      select(Name, Plot_code, all_of(traits_filtered))
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name, Plot_code) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup()
    
    input_trait_genus_fd <- input_trait_genus
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F)))) 
  }
  
  input_veg_genus <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    column_to_rownames("Name") %>%
    t()
  
  input_veg_genus_fd <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    filter(Name %in% input_trait_genus$Name) %>%
    column_to_rownames("Name") %>%
    t()
  
  if(trait_detail == "Study"){
    input_trait_genus_fd <- input_trait_genus_fd %>%
      filter(row.names(.) %in%  colnames(input_veg_genus_fd))
  } else{
    input_trait_genus_fd <- input_trait_genus_fd %>%
      filter(Name %in% colnames(input_veg_genus_fd))
  }
  
  prep_orig <- cats_preparation(input_veg = input_veg, 
                                input_trait = input_trait,
                                trait_information = trait_information, 
                                cutoff_value = cutoff_value, 
                                trait_detail = trait_detail)
  
  prep_orig_genus <- cats_preparation(input_veg = input_veg_genus, 
                                      input_trait = input_trait_genus,
                                      trait_information = trait_information, 
                                      cutoff_value = cutoff_value, 
                                      trait_detail = trait_detail)
  
  diversity_info_species <- calculate_diversity_info(veg_data = input_veg, 
                                                     input_veg_fd = input_veg_fd,
                                                     input_processed = prep_orig, 
                                                     trait_fd = input_trait_fd,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)
  
  diversity_info_genus <- calculate_diversity_info(veg_data = input_veg_genus,
                                                   input_veg_fd = input_veg_genus_fd,
                                                   input_processed = prep_orig_genus,
                                                   trait_fd = input_trait_genus_fd,
                                                   trait_detail = trait_detail,
                                                   study_name = study_name)
  
  input_trait_try <- trait_imp_specieslvl %>% 
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]])) 
  
  input_trait_try_genus <- trait_imp_genuslvl %>%
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]]))
  
  prep_try <- cats_preparation(input_veg = input_veg, 
                               input_trait = input_trait_try,
                               trait_information = trait_information, 
                               cutoff_value = 0.6, 
                               trait_detail = "Study")
  
  prep_try_genus <- cats_preparation(input_veg = input_veg_genus, 
                                     input_trait = input_trait_try_genus,
                                     trait_information = trait_information, 
                                     cutoff_value = 0.6, 
                                     trait_detail = "Study")
  
  if(foreach::getDoParRegistered() != TRUE){
    n_cores <- parallel::detectCores() - 5 
    
    my_cluster <- makeSOCKcluster(n_cores)
    
    registerDoSNOW(cl = my_cluster)}   
  
  
  cats_orig <- cats_calculations(prep_orig, trait_information, 
                                 trait_detail, nperm = 400) %>%
    mutate(study = study_name)
  
  cats_try <- cats_calculations(prep_try, trait_information, 
                                trait_detail = "Study", nperm = 400) %>%
    mutate(study = study_name)
  
  cats_orig_genus <- cats_calculations(prep_orig_genus, trait_information, 
                                       trait_detail, nperm = 400) %>%
    mutate(study = study_name)
  
  cats_try_genus <- cats_calculations(prep_try_genus, trait_information, 
                                      trait_detail = "Study", nperm = 400) %>%
    mutate(study = study_name)
  
  
  all_ouput <- lst(diversity_info_species, diversity_info_genus, prep_orig, 
                   prep_try, cats_orig,  cats_try,
                   prep_orig_genus, prep_try_genus, cats_orig_genus, cats_try_genus) 
  
  return(all_ouput)
  
}

test_function2 <- function(input_veg, input_veg_fd, input_trait, input_trait_fd, trait_information, 
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
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      filter(Name %in% trait_imp_genuslvl$Name)
      
    
    input_trait_genus_fd <- input_trait_genus %>%
      column_to_rownames(var = "Name")
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
    
    
  } else{
    input_trait <- input_trait %>%
      select(Name, Plot_code, all_of(traits_filtered))
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name, Plot_code) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      filter(Name %in% trait_imp_genuslvl$Name)
    
    input_trait_genus_fd <- input_trait_genus
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)
  }
  
  input_veg_genus <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    column_to_rownames("Name") %>%
    t()
  
  input_veg_genus_fd <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    filter(Name %in% input_trait_genus$Name) %>%
    column_to_rownames("Name") %>%
    t() %>%
    as.data.frame()
  
  if(trait_detail == "Study"){
    input_trait_genus_fd <- input_trait_genus_fd %>%
      filter(row.names(.) %in%  colnames(input_veg_genus_fd))
  } else{
    input_trait_genus_fd <- input_trait_genus_fd %>%
      filter(Name %in% colnames(input_veg_genus_fd))
  }
  
  prep_orig <- cats_preparation(input_veg = input_veg, 
                                input_trait = input_trait,
                                trait_information = trait_information, 
                                cutoff_value = 0.6, 
                                trait_detail = trait_detail)
  
  prep_orig_genus <- cats_preparation(input_veg = input_veg_genus, 
                                      input_trait = input_trait_genus,
                                      trait_information = trait_information, 
                                      cutoff_value = 0.6, 
                                      trait_detail = trait_detail)
  
  input_veg_genus_fd <- input_veg_genus_fd %>%
    filter(row.names(.) %in% rownames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]]))
  
  input_veg_genus_fd <- input_veg_genus_fd %>%
    select(all_of(colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]])))
  
  if(trait_detail == "Study"){
  input_trait_genus_fd <- input_trait_genus_fd %>%
    filter(row.names(.) %in% colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]]))
  } else{
    input_trait_genus_fd <- input_trait_genus_fd %>%
      filter(Name %in% colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]]))
  }

  diversity_info_species <- calculate_diversity_info(veg_data = input_veg, 
                                                     input_veg_fd = input_veg_fd,
                                                     input_processed = prep_orig, 
                                                     trait_fd = input_trait_fd,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)
  

  diversity_info_genus <- calculate_diversity_info(veg_data = input_veg_genus,
                                                   input_veg_fd = input_veg_genus_fd,
                                                   input_processed = prep_orig_genus,
                                                   trait_fd = input_trait_genus_fd,
                                                   trait_detail = trait_detail,
                                                   study_name = study_name)

  
  input_trait_try <- trait_imp_specieslvl %>% 
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]])) 
  
  diversity_info_try_species <- calculate_diversity_info(veg_data = input_veg, 
                                                     input_veg_fd = input_veg_fd,
                                                     input_processed = prep_try, 
                                                     trait_fd = input_trait_try,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)
  
  input_trait_try_genus <- trait_imp_genuslvl %>%
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]]))
  
  
  diversity_info_try_genus <- calculate_diversity_info(veg_data = input_veg_genus, 
                                                         input_veg_fd = input_veg_genus_fd,
                                                         input_processed = prep_try_genus, 
                                                         trait_fd = input_trait_try_genus,  
                                                         trait_detail = trait_detail, 
                                                         study_name = study_name)
  
  prep_try <- cats_preparation(input_veg = input_veg, 
                               input_trait = input_trait_try,
                               trait_information = trait_information, 
                               cutoff_value = 0.6, 
                               trait_detail = "Study")
  
  prep_try_genus <- cats_preparation(input_veg = input_veg_genus, 
                                     input_trait = input_trait_try_genus,
                                     trait_information = trait_information, 
                                     cutoff_value = 0.6, 
                                     trait_detail = "Study")
  
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
  
  if(!all.equal(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]], 
                prep_try_genus[["relative_abundances"]][["matrix4traits_prop"]])) {
    warning("Vegetation matrix is not the same")
  }

  if(foreach::getDoParRegistered() != TRUE){
    n_cores <- parallel::detectCores() - 5 
    
    my_cluster <- makeSOCKcluster(n_cores)
    
    registerDoSNOW(cl = my_cluster)}   
  
  cats_orig <- cats_calculations(prep_orig, trait_information, trait_detail) %>%
    mutate(study = study_name)
  
  cats_try <- cats_calculations(prep_try, trait_information, trait_detail = "Study") %>%
    mutate(study = study_name)
  
  cats_orig_genus <- cats_calculations(prep_orig_genus, trait_information, 
                                       trait_detail) %>%
    mutate(study = study_name)
  
  cats_try_genus <- cats_calculations(prep_try_genus, trait_information, 
                                      trait_detail = "Study") %>%
    mutate(study = study_name)
  
  
  all_ouput <- lst(diversity_info_species, diversity_info_genus, 
                   diversity_info_try_species, diversity_info_try_genus, 
                   prep_orig, prep_try, cats_orig,  cats_try,
                   prep_orig_genus, prep_try_genus, cats_orig_genus, cats_try_genus) 
  
  return(all_ouput)
  
}



test_function3 <- function(input_veg, input_veg_fd, input_trait, input_trait_fd, trait_information, 
                           cutoff_value, trait_detail, study_name){
  
  extract <- new_names_synonyms %>%
    filter(study == study_name)
  
  not_included_traits <- c("Leaf_angle", "Photo_path", "Root_density", "Lateral_spread",
                           "Resprout", "Bark_thickness_rel", "Rhizome", "Stem_emerg",
                           "ShootDMC", "Leaf_shape", "twFWC", "FWC")
  
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
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      filter(Name %in% trait_imp_genuslvl$Name)
    
    input_trait_genus_fd <- input_trait_genus %>%
      column_to_rownames(var = "Name")
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait_fd <- input_trait %>%
      filter(Name %in% trait_imp_specieslvl$Name) %>%
      column_to_rownames(var = "Name")
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)

    
  } else{
    input_trait <- input_trait %>%
      select(Name, Plot_code, all_of(traits_filtered))
    
    input_trait_genus <- input_trait %>%
      separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
      group_by(Name, Plot_code) %>%
      summarise(across(where(is.numeric), mean)) %>%
      ungroup() %>%
      filter(Name %in% trait_imp_genuslvl$Name)
    
    input_trait_genus_fd <- input_trait_genus
    
    input_trait_genus <- input_trait_genus %>%
      mutate(across(where(is.numeric), ~c(scale(., center = F))))
    
    input_trait_fd <- input_trait %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
    input_trait <- input_trait %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = F))))  %>%
      filter(Name %in% trait_imp_specieslvl$Name)
    
  }
  
  input_veg_genus <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    column_to_rownames("Name") %>%
    t()
  
  input_veg_genus_fd <- input_veg %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Species") %>%
    separate(Species, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), sum)) %>%
    ungroup() %>%
    filter(Name %in% input_trait_genus$Name) %>%
    column_to_rownames("Name") %>%
    t() %>%
    as.data.frame()
  
  prep_orig <- cats_preparation(input_veg = input_veg, 
                                input_trait = input_trait,
                                trait_information = trait_information, 
                                cutoff_value = 0.6, 
                                trait_detail = trait_detail)
  
  prep_orig_genus <- cats_preparation(input_veg = input_veg_genus, 
                                      input_trait = input_trait_genus,
                                      trait_information = trait_information, 
                                      cutoff_value = 0.6, 
                                      trait_detail = trait_detail)

  diversity_info_species <- calculate_diversity_info(veg_data = input_veg, 
                                                     input_veg_fd = input_veg_fd,
                                                     input_processed = prep_orig, 
                                                     trait_fd = input_trait_fd,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)
  
  diversity_info_genus <- calculate_diversity_info(veg_data = input_veg_genus,
                                                   input_veg_fd = input_veg_genus_fd,
                                                   input_processed = prep_orig_genus,
                                                   trait_fd = input_trait_genus_fd,
                                                   trait_detail = trait_detail,
                                                   study_name = study_name)
  
  
  input_trait_try <- trait_imp_specieslvl %>% 
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig[["relative_abundances"]][["matrix4traits_prop"]])) 
  
  
  input_trait_try_genus <- trait_imp_genuslvl %>% 
    select(c(Name, traits_filtered)) %>%
    filter(Name %in% 
             colnames(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]])) 

  
  
  prep_try <- cats_preparation(input_veg = input_veg, 
                               input_trait = input_trait_try,
                               trait_information = trait_information, 
                               cutoff_value = 0.6, 
                               trait_detail = "Study")
  
  prep_try_genus <- cats_preparation(input_veg = input_veg_genus, 
                                     input_trait = input_trait_try_genus,
                                     trait_information = trait_information, 
                                     cutoff_value = 0.6, 
                                     trait_detail = "Study")
  

    input_trait_try_fd <- input_trait_try %>%
      column_to_rownames(var = "Name") %>%
      filter(row.names(.) %in% colnames(prep_try[["relative_abundances"]][["matrix4traits_prop"]]))
    
    input_trait_try_genus_fd <- input_trait_try_genus %>%
      column_to_rownames(var = "Name") %>%
      filter(row.names(.) %in% colnames(prep_try_genus[["relative_abundances"]][["matrix4traits_prop"]]))

  
  diversity_info_try_species <- calculate_diversity_info(veg_data = input_veg, 
                                                         input_veg_fd = input_veg_fd,
                                                         input_processed = prep_try, 
                                                         trait_fd = input_trait_try_fd,  
                                                         trait_detail =  "Study", 
                                                         study_name = study_name)
  

  diversity_info_try_genus <- calculate_diversity_info(veg_data = input_veg_genus, 
                                                       input_veg_fd = input_veg_genus_fd,
                                                       input_processed = prep_try_genus, 
                                                       trait_fd = input_trait_try_genus_fd,  
                                                       trait_detail = "Study", 
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
  
  if(!all.equal(prep_orig_genus[["relative_abundances"]][["matrix4traits_prop"]], 
                prep_try_genus[["relative_abundances"]][["matrix4traits_prop"]])) {
    warning("Vegetation matrix is not the same")
  }
  
  if(foreach::getDoParRegistered() != TRUE){
    n_cores <- parallel::detectCores() - 5 
    
    my_cluster <- makeSOCKcluster(n_cores)
    
    registerDoSNOW(cl = my_cluster)}   
  
  cats_orig <- cats_calculations(prep_orig, trait_information, trait_detail, 
                                 nperm = 300) %>%
    mutate(study = study_name)
  
  cats_try <- cats_calculations(prep_try, trait_information, 
                                trait_detail = "Study", nperm = 300) %>%
    mutate(study = study_name)
  
  cats_orig_genus <- cats_calculations(prep_orig_genus, trait_information, 
                                       trait_detail, nperm = 300) %>%
    mutate(study = study_name)
  
  cats_try_genus <- cats_calculations(prep_try_genus, trait_information, 
                                     trait_detail = "Study", nperm = 300) %>%
    mutate(study = study_name)
  

  all_ouput <- lst(diversity_info_species, diversity_info_genus, 
                   diversity_info_try_species, diversity_info_try_genus, 
                   prep_orig, prep_try, prep_orig_genus, prep_try_genus,
                   cats_orig,  cats_try,
                   cats_orig_genus, cats_try_genus) 
  
  return(all_ouput)
  
}
  
  
  
  
  

