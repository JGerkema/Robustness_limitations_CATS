complete_analysis <- function(input_veg, input_veg_fd, input_trait, input_trait_fd, trait_information, 
                         cutoff_value, trait_detail, study_name){
  #' @title Necessary preparations and calculations to create input-ready data for the CATS model and the 
  #' @description 
  #' @param input_veg The vegetation input should be a matrix where each row represents a plot and each column a species. The species columns should be alphabetically ordered. 
  #' @param input_trait The trait input should be a dataframe where each row represents a species. The species names should be stored in a column called "Name". The corresponding trait values should be stored in the rest of the columns, in alphabetical order. 
  #' @param trait_information A list with the traits used in the study
  #' @param cutoff_value In some occurrences, not all species will have trait data. These species are removed. The total relative abundance of species with trait data is then calculated. If it exceeds the specified cutoff value, the plot remains in the data. Otherwise it is removed. 
  #' @param trait_detail 
  #' @returns A list with the processed input, diversity information and output of the CATS model.  
  #'
  #'

  input_processed_orig <- cats_preparation(input_veg = input_veg, 
                                      input_trait = input_trait,
                                      trait_information = trait_information, 
                                      cutoff_value = cutoff_value, 
                                      trait_detail = trait_detail)

input_trait_try <- TRY_DATA %>% 
    select(c(colnames(input_trait))) %>%
    filter(Name %in% colnames(temp))

  
  input_processed_try <- cats_preparation(input_veg = input_veg, 
                                          input_trait = input_trait_try,
                                          trait_information = trait_information, 
                                          cutoff_value = cutoff_value, 
                                          trait_detail = trait_detail)

  diversity_info <- calculate_diversity_info(veg_data = input_veg, 
                                             veg_data_fd = input_veg_fd,
                                             input_processed = input_processed, 
                                             trait_fd = input_trait_fd,  
                                             trait_detail = trait_detail, 
                                             study_name = study_name)
  
  
  if(foreach::getDoParRegistered() != TRUE){
    n_cores <- parallel::detectCores() - 5 
    
    my_cluster <- makeSOCKcluster(n_cores)
    
    registerDoSNOW(cl = my_cluster)}   

  
  cats_output <- cats_calculations(input_processed, trait_information, trait_detail)

  cats_output_clean <- cats_visualisation(cats_output) %>%
                      mutate(study = study_name)

  all_ouput <- lst(diversity_info, input_processed, cats_output, cats_output_clean)
  
  return(all_ouput)
  
}





