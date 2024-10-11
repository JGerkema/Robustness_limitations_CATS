data_summary <- data.frame()

for(i in 1:length(data_input_list)){
  
  dataset <- data_input_list[[i]]
  
  extract <- new_names_synonyms %>%
  filter(study == dataset[["study_name"]])

not_included_traits <- c("Leaf_angle", "Photo_path", "Root_density", "Lateral_spread",
                         "Resprout", "Bark_thickness_rel", "Rhizome", "Stem_emerg",
                         "ShootDMC", "Leaf_shape", "Root_diam", "twFWC", "FWC")

traits_orig <- dataset[["trait_information"]][["traitlist"]]

traits_filtered <- traits_orig[traits_orig %notin% not_included_traits]  

traitlist <- sort(c(traits_filtered))
cwm_names <- paste0(traitlist, "_cwm")
cwv_names <- paste0(traitlist, "_cwv")

trait_information <- lst(traitlist, cwm_names, cwv_names)

if(nrow(extract) > 0){
  
  input_veg <- dataset[["input_veg"]] %>%
    rename(!!!setNames(extract$Orig_name, extract$Name)) %>%
    select(order(colnames(.)))
  
  input_trait <- dataset[["input_trait"]]  %>%
    mutate(Name = case_when(
      Name %in% extract$Orig_name ~ extract$Name[match(Name, extract$Orig_name)],
      TRUE ~ Name)) %>%
    arrange(Name)
}else{
  input_veg <- dataset[["input_veg"]] 
  input_trait <- dataset[["input_trait"]]
}


if(dataset[["trait_detail"]] == "Study"){
  
  input_trait <- input_trait %>%
    select(Name, all_of(traits_filtered))
  
  input_trait_genus <- input_trait %>%
    separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup() %>%
    filter(Name %in% trait_imp_genuslvl$Name)

    mutate(across(where(is.numeric), ~c(scale(., center = F))))

  
  
} else{
  input_trait <- input_trait %>%
    select(Name, Plot_code, all_of(traits_filtered))
  
  input_trait_genus <- input_trait %>%
    separate(Name, into = c("Name"), extra = "drop", sep = " ") %>%
    group_by(Name, Plot_code) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup() %>%
    filter(Name %in% trait_imp_genuslvl$Name)
  
  input_trait <- input_trait %>%
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



prep_orig <- cats_preparation(input_veg = input_veg, 
                              input_trait = input_trait,
                              trait_information = trait_information, 
                              cutoff_value = 0.6, 
                              trait_detail = dataset[["trait_detail"]])


n_plot_before <- ncol(prep_orig[["relative_abundances"]][["matrix"]])

n_plot_after <- ncol(prep_orig[["relative_abundances"]][["matrix4traits_prop"]])


regional_trait_coverage <- prep_orig[["relative_abundances"]][["matrix"]] %>% 
  as.data.frame() %>%
  summarise(across(where(is.numeric), ~ sum(.))) %>%
  t() %>%
  as.data.frame() %>%
  mutate(Total_cover = rowSums(.)) %>%
  select(Total_cover) %>%
  filter(Total_cover != 0) %>%
  mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
  rownames_to_column(var = "Name") %>%
  filter(Name %in% input_trait$Name) %>%
  select(Name, Relative_meta_cover) %>%
  summarise(across(where(is.numeric), ~sum(.))) %>%
  mutate(study = dataset[["study_name"]])

prep_orig_genus <- cats_preparation(input_veg = input_veg_genus, 
                                    input_trait = input_trait_genus,
                                    trait_information = trait_information, 
                                    cutoff_value = 0.6, 
                                    trait_detail = dataset[["trait_detail"]])

plot_count <- bind_cols(n_plot_before = n_plot_before, 
                        n_plot_after =  n_plot_after,
                        study = dataset[["study_name"]]) %>%
  full_join(regional_trait_coverage) 

data_summary <- bind_rows(data_summary, plot_count)

} 



 