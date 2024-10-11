

veg_data <- calculate_abundance(chalmandrier_cover, chalmandrier_trait, 0.8)
trait_data <- chalmandrier_trait

cwm_cwv_all_traits <- as.data.frame(plot_code_vector) %>%
  rename(Plot_code = plot_code_vector)



for(t in 1:length(traitlist)) {
  
  trait <- traitlist[t]
  cwm_cwv_one_trait <- data.frame()
  trait_cwm <- cwm_names[t]
  trait_cwv <- cwv_names[t]
  
  for(p in plot_code_vector){
    # Calculate 'community weighted means' (CWMs) by multiplying species proportional abundances by matched trait values and summing by site
    
    # Select one specific plot
    veg_perplot <- as.data.frame(veg_data$matrix4traits_prop) %>%
      filter(rownames(.) %in% p) 
    
    # Select the trait data of the same plot
    trait_perplot <- trait_data %>%
      filter(Plot_code == p)
    
    # The trait data contains only the data of present species while the veg data
    # also contains empty states, ie. species with a relative abundance of 0.
    # By matching trait data with the species list of the respective plot, 
    # empty states are also created in the trait data. 
    trait_perplot_temp <- trait_perplot[[trait]][match(
      colnames(veg_perplot),
      trait_perplot$Name)]
    
    # To ensure a 0*0 multiplication, the NA's created in the previous step
    # are substituted with 0.
    trait_perplot_temp[is.na(trait_perplot_temp)] <- 0
    
    # The community weighted mean is calculated by multiplying the relative
    # abundance of each trait with it respective value.
    cwm <- as.data.frame((as.matrix(veg_perplot) %*% trait_perplot_temp)) %>%
      rename({{trait_cwm}} := V1) 
    

    
    veg_for_cwv <- veg_perplot %>%
      pivot_longer(cols = 1:ncol(.), 
                   values_to = "abundance", names_to = "Name") %>%
      filter(abundance > 0) 
    
    #
    trait_for_cwv <- trait_perplot %>%
      filter(Name %in% veg_for_cwv$Name) %>%
      pull({{trait}})
    
    veg_for_cwv <- veg_for_cwv %>%
      pull(abundance) %>% as.matrix()
      
    trait_variance <- (trait_for_cwv - cwm[1,1])^2 
  
    #The brackets around veg_for_cwv * trait_variance are essential. 
    # Otherwise it starts multiplying and colSumming simultaneously instead or 
    # sequentially. Not good. 
    cwv <- (veg_for_cwv * trait_variance) %>%
      colSums() %>%
      as.data.frame() %>%
      rename({{trait_cwv}} := .)

    cwm_cwv <- bind_cols(cwm, cwv)
    
    cwm_cwv_one_trait <- bind_rows(cwm_cwv_one_trait, cwm_cwv) 
    
    
    
  }
  
  cwm_cwv_one_trait <- cwm_cwv_one_trait %>%
    rownames_to_column(var = "Plot_code")  
  
  cwm_cwv_all_traits <- full_join(cwm_cwv_all_traits, cwm_cwv_one_trait)
  
} 
