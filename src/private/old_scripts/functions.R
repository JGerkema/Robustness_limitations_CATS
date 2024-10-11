#--------------------------------------------------------#
#   Useful functions Chapter I
#
#   Jelyn Gerkema (May 2023)
#--------------------------------------------------------#

#To do: put in warning messages

######## Set Not in 
`%notin%` <- Negate(`%in%`)



######## Calculate relative abundances
calculate_abundance <- function(input_veg_data, input_trait_data, cutoff_value) {
  #' @title Calculate the relative species abundance of each plot.
  #' @description The relative species abundance is calculated for each plot. Only sites were the species for which trait data is available 
  #' cover more than 80% are retained in the dataset. Species that don't occur when the dataset is trimmed to 80% coverage are removed. 
  #' The abundance data of the remaining sites is reproportioned so the abundances sum to 1 for each site. 
  #' Code adapted from: Guerin et al. 'Environmental correlates of key functional traits in Australian plant communities based on relative abundance in plots'.
  #' @param input_veg_data The vegetation input data. Should be formatted as a matrix where each row is a single site and each column corresponds with a single species.
  #' @param input_trait_data The trait input data. Should be formatted as a dataframe where each row contains the trait values of a single
  #' entities (species, genus, family etc). Entity names should be stored in column called 'Name'.
  #' @param cutoff_value Determines the coverage threshold value below which plots get removed from the dataset. 
  #' @returns A list containing the output of the steps explained above. The final matrix with the proportional abundances of the remaining sites is 
  #' returned as 'matrix4traits_prop'.
  
  
  veg_data <- list()
  
  veg_data$matrix <- input_veg_data
  
  trait_data <- input_trait_data
  
  # add a matrix of proportional abundances (abundances divided by total site abundances)
  veg_data$matrix_prop <- veg_data$matrix / rowSums(veg_data$matrix)
  
  # add a matrix of proportional abundances that is cut down to species that are matched in trait the data
  veg_data$matrix_prop_match <- veg_data$matrix_prop[, colnames(veg_data$matrix_prop) %in% 
                                trait_data$Name]
  
  # cut down to high coverage sites by abundance (x%) now that species with no trait data have been removed
  veg_data$matrix_prop_match <- veg_data$matrix_prop_match[rowSums(veg_data$matrix_prop_match) > cutoff_value, ]
  
  # remove also species that don't occur in the dataset when trimmed to x% coverage
  veg_data$matrix_prop_match <- veg_data$matrix_prop_match[, which(colSums(veg_data$matrix_prop_match) > 0)]
  
  if(ncol(veg_data$matrix_prop_match) == 0){
  warning("No plots above the cutoff value")
  return(NULL)  # Break off the function execution
}
  
  veg_data$matrix4traits <- veg_data$matrix[rownames(veg_data$matrix) %in% 
                            rownames(veg_data$matrix_prop_match), colnames(veg_data$matrix) %in% 
                            colnames(veg_data$matrix_prop_match)]
  
  # convert the above into proportional abundance based on only species included in the analysis so that abundances sum to 1 for each community
  veg_data$matrix4traits_prop <- veg_data$matrix4traits / rowSums(veg_data$matrix4traits)
  
  return(veg_data)
} 



######## Calculate community weighted means and community weighted variance
calculate_cwm_cwv <- function(input_veg_data, input_trait_data, trait, trait_cwm, trait_cwv) {
  #' @title Calculate community weighted mean (cwm) and variance (cwv) of each plot
  #' @description The proportional abundance of each species in each site is 
  #' multiplied by the respective trait value. For each site, the multiplied 
  #' products are summed together, resulting in the community weighted mean.
  #' @param input_veg_data The dataset for which the cwm and cwv should be calculated. 
  #' The input data should be the output of the 'calculate_abundance' function.
  #' @param input_trait_data Dataframe containing the trait data.
  #' @param trait The name of the trait
  #' @param trait_cwm The name under which the cwm's should be stored.
  #' @param trait_cwv The name under which the cwv's should be stored.
  #' @returns A new dataset which the cwm and cwv of each site.
  
  # Calculate 'community weighted means' (CWMs) by multiplying species proportional abundances by matched trait values and summing by site
  veg_data <- input_veg_data
  trait_data <- input_trait_data
  

  tmp <- trait_data[[trait]][match(
    colnames(veg_data$matrix4traits_prop),
    trait_data$Name)]
  cwm <- as.data.frame((as.matrix(veg_data$matrix4traits_prop) %*%
                          trait_data[[trait]][match(
                            colnames(veg_data$matrix4traits_prop),
                            trait_data$Name
                          )])) %>% 
    # The function 'match' operates based on the order of the first object 
    #of the function, in this case the order of the column names. 
    rename({{ trait_cwm }} := V1) %>%
    rownames_to_column(var = "Plot_code") 
  
  veg_data$TEMP <- veg_data$matrix4traits_prop # copy the proportional abundance matrix to store variable components
  
  for (site_row in 1:nrow(veg_data$TEMP)) { # for each site
    for (species_col in 1:ncol(veg_data$TEMP)) { # for each species
      prop <- veg_data$TEMP[site_row, species_col] # store the abundance of species in site 
      if (prop > 0) { # if the species occurs
        CWM <- cwm[site_row, 2] # store the relevant CWM for that site from above
        value <- trait_data[[trait]][match(colnames(veg_data$TEMP)[species_col], 
                                                    trait_data$Name)]
        # extract the trait value that matches the current species
        veg_data$TEMP[site_row, species_col] <- prop * (value - CWM)^2
      }
    } # end species_col
  } # end site_row
  
  cwv <- as.data.frame(as.matrix(rowSums(veg_data$TEMP))) %>%
    rename({{ trait_cwv }} := V1) %>%
    rownames_to_column(var = "Plot_code") 
  
  output <- inner_join(cwm, cwv, by = "Plot_code")
  
  return(output)
}


cats_preparation <- function(input_veg, input_trait, trait_information, cutoff_value, trait_detail){
  #' @title Necessary preparations and calculations to create input-ready data for the CATS model. 
  #' @description First, the raw abundance data is converted to relative abundances. 
  #' Only the plots with enough trait coverage are retained in the dataset using the 'calculate_abundance' function. 
  #' Then, the community weighted mean (cwm) and variance (cwv) of each plot for multiple traits is calculated using the 'calculate_cwm_cwv' function.
  #' @param input_veg The vegetation input should be a matrix where each row represents a plot and each column a species. The species columns should be alphabetically ordered. 
  #' @param input_trait The trait input should be a dataframe where each row represents a species. The species names should be stored in a column called "Name". The corresponding trait values should be stored in the rest of the columns, in alphabetical order. 
  #' @param trait_information 
  #' @param cutoff_value In some occurrences, not all species will have trait data. These species are removed. The total relative abundance of species with trait data is then calculated. If it exceeds the specified cutoff value, the plot remains in the data. Otherwise it is removed. 
  #' @param trait_detail
  #' @returns A list with the relative abundances, trait data with only the species actually occurring in the dataset and the calculated CWM and CWV values for each plot. 
  #'
  #'
  if (!all(colnames(input_veg) == sort(colnames(input_veg)))) {
    stop("Column names input_veg (species names) are not in alphabetical order")
  } # end If loop
  
  if(identical(trait_information[["traitlist"]], sort(trait_information[["traitlist"]])) == FALSE){
    stop("List with trait names is not in alphabetical order")
  }
  
  if(identical(trait_information[["cwm_names"]], sort(trait_information[["cwm_names"]])) == FALSE){
    stop("List with trait_cwm names is not in alphabetical order")
  }
  
  if(identical(trait_information[["cwv_names"]], sort(trait_information[["cwv_names"]])) == FALSE){
    stop("List with trait_cwv names is not in alphabetical order")
  }
  
  cutoff_value <- cutoff_value
  
  if(trait_detail == "Study"){ 
    
    check_trait <- input_trait %>%
      column_to_rownames(var = "Name") 
    
    
    if (!all(colnames(check_trait) == sort(colnames(check_trait)))) {
      stop("Column names input_trait (trait names) are not in alphabetical order")
    } 
    
    if (!all(rownames(check_trait) == sort(rownames(check_trait)))) {
      stop("Row names input_trait (species names) are not in alphabetical order")
    } 
    
    
    relative_abundances <- calculate_abundance(input_veg, input_trait, cutoff_value)
    
    if(is.null(relative_abundances) == TRUE){
      warning("Relative abundances is empty")
      return(NULL)
    }
    # Filter the trait species so that they match up with the species actually occurring in the dataset
    trait_data <- input_trait %>%
      subset(Name %in% colnames(relative_abundances$matrix_prop_match)) %>%
      arrange(Name)
    
    calculations_cwm_cwv <- as.data.frame(relative_abundances[["matrix4traits_prop"]]) %>%
      rownames_to_column(var = "Plot_code") %>%
      select(Plot_code)
    
    # Using the function 'calculate_cwm_cwv' for each trait 
    for (t in 1:length(trait_information[["traitlist"]])) {
      
      trait <- trait_information[["traitlist"]][t]
      trait_cwm <- trait_information[["cwm_names"]][t]
      trait_cwv <- trait_information[["cwv_names"]][t]
      
      output <- calculate_cwm_cwv(relative_abundances, trait_data, {{ trait }}, {{ trait_cwm }}, 
                                  {{ trait_cwv }})
      
      calculations_cwm_cwv <- calculations_cwm_cwv %>% 
        inner_join(output, by = "Plot_code")
    } 
    
    
    
  }  else {
    
    
    relative_abundances <- calculate_abundance(input_veg, input_trait, cutoff_value )
    
    if(is.null(relative_abundances) == TRUE){
      warning("Relative abundances is empty")
      return(NULL)
    }
    
    # Select only the trait data for species that are present in the dataset after it has been filtered
    trait_data <- input_trait %>%
      filter(Name %in% colnames(relative_abundances$matrix4traits_prop)) %>%
      arrange(Name)
    
    plot_code_vector <- rownames(relative_abundances$matrix4traits_prop)
    
    calculations_cwm_cwv <- as.data.frame(plot_code_vector) %>%
      rename(Plot_code = plot_code_vector)
    
    for(t in 1:length(trait_information[["traitlist"]])) {
      
      cwm_cwv_one_trait <- data.frame()
      
      trait <- trait_information[["traitlist"]][t]
      trait_cwm <- trait_information[["cwm_names"]][t]
      trait_cwv <- trait_information[["cwv_names"]][t]
      
      for(p in plot_code_vector){
        # Calculate 'community weighted means' (CWMs) by multiplying species proportional abundances by matched trait values and summing by site
        
        # Select one specific plot
        veg_perplot <- as.data.frame(relative_abundances$matrix4traits_prop) %>%
          filter(rownames(.) %in% p) 
        
        # Select the trait data of the same plot
        trait_perplot <- trait_data %>%
          filter(Plot_code == p)
        
        trait_perplot_check <- trait_perplot %>%
          select(-c(Plot_code)) %>%
          column_to_rownames(var = "Name")
        
        if (!all(colnames(trait_perplot_check) == sort(colnames(trait_perplot_check)))) {
          stop("Column names input_trait (trait names) are not in alphabetical order")
        } 
        
        if (!all(rownames(trait_perplot_check) == sort(rownames(trait_perplot_check)))) {
          stop("Row names input_trait (species names) are not in alphabetical order")
        }         
        
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
        # abundance of each trait with it respective value. With matrix multiplications, 
        # 
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
        # Otherwise it starts multiplying and colSumming simultaneously instead of
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
      
      calculations_cwm_cwv <- full_join(calculations_cwm_cwv, cwm_cwv_one_trait)
      
    }
    
    
  } # end else loop
  
  output <- lst(calculations_cwm_cwv, relative_abundances, trait_data)
  
  return(output)
  
  
}



######## Calculate the average Bray-Curtis dissimilarity index for each plot

calculate_mean_dist <- function(input_veg, input_processed){

all_dist <- as.matrix(vegdist(input_processed[["relative_abundances"]][["matrix4traits"]], method="bray"))

plot_code_vector <- rownames(input_processed[["relative_abundances"]][["matrix4traits_prop"]])

all_dist <- all_dist %>%
  as.data.frame() %>%
  select(all_of(plot_code_vector)) %>%
  filter(rownames(.) %in% plot_code_vector)

mean_dist <- data.frame()

for(p in plot_code_vector){
  
  dist_oneplot <- all_dist %>%
    as.data.frame() %>%
    select(all_of(p)) %>%
    filter(rownames(.) %notin% c(p)) %>%
    summarise(across(where(is.numeric), ~mean(., na.rm = T))) %>%
    t() %>%
    as.data.frame() %>%
    rename(mean_dist = V1)
  
  mean_dist <- bind_rows(mean_dist, dist_oneplot)
  
    }

mean_dist <- mean_dist %>%
  rownames_to_column(var = "Plot_code")

return(mean_dist)

}

######## Calculate the functional diversity indices of each plot


calculate_FD_plotlvl <- function(input_processed, input_veg, input_trait){
  
  plot_code_vector <- rownames(input_processed[["relative_abundances"]][["matrix4traits_prop"]])
  
  FD_df <- data.frame(matrix(nrow = length(plot_code_vector), ncol = 6))
  
  assigned_colnames <- c("FEve", "FDiv", "FRic", "RaoQ", "qual.FRic",  "Plot_code")
  colnames(FD_df) <- assigned_colnames
  m <- 0
  
  for(p in plot_code_vector){
    
    m <- m + 1
    veg_temp <- input_veg %>%
      subset(rownames(.) %in% p) %>%
      as.data.frame() %>% 
      select(where(~ any(. != 0))) # removes the columns (species) with abundance 0
    
    trait_temp <- input_trait %>%
      filter(Plot_code %in% p) %>%
      filter(Name %in% colnames(veg_temp)) %>%
      column_to_rownames(var = "Name") %>%
      select(-c(Plot_code)) %>%
      mutate(across(where(is.numeric), ~ c(scale(., center = T)))) 
    
    if(ncol(veg_temp) < 3) {
     next  # If the amount of present species is less than 3, the algorithm is 
      # not able to compute the FD. Therefore, these plots are skipped. 
    }
    
    else {
    FD <- dbFD(trait_temp, veg_temp, stand.x = TRUE)
    
    
    FD_df[m, 1] <- FD$FEve
    FD_df[m, 2] <- FD$FDiv
    FD_df[m, 3] <- FD$FRic
    FD_df[m, 4] <- FD$RaoQ  
    FD_df[m, 5] <- FD$qual.FRic
    FD_df[m, 6] <- p
    
      }
  }
    
  return(FD_df)
  
}

calculate_diversity_info <- function(veg_data, input_veg_fd, input_processed, trait_fd,  
                                     trait_detail, study_name){
  
  mean_dist <- calculate_mean_dist(veg_data, input_processed)
  
  if(trait_detail == "Study"){
    
    FD <- dbFD(trait_fd, input_veg_fd, stand.x = TRUE)
    
    FD_output<- data.frame(FEve = FD$FEve,
                          FDiv = FD$FDiv,
                          FRic = FD$FRic,
                          RaoQ = FD$RaoQ,
                          qual.FRic = FD$qual.FRic) %>%
      rownames_to_column(var = "Plot_code")
    
    

  } else {
    
    FD_output <- calculate_FD_plotlvl(input_processed, input_veg_fd, trait_fd)

  }
  
  diversity_info <- input_processed[["relative_abundances"]][["matrix4traits_prop"]] %>% 
    as.data.frame() %>% 
    mutate(shannon = diversity(.[,1:ncol(.)], index = "shannon", MARGIN = 1, base = exp(1)),
           simpson = diversity(.[,1:ncol(.)], index = "simpson", MARGIN = 1, base =exp(1)),
           species_richness = specnumber(.[,1:ncol(.)], MARGIN = 1),
           empty_states = ncol(input_processed[["relative_abundances"]][["matrix4traits_prop"]]) - species_richness,
           meta_size = ncol(input_processed[["relative_abundances"]][["matrix4traits_prop"]])) %>%
    select(shannon, simpson, species_richness, empty_states, meta_size) %>% # base = exp(1) gives the shannon index a log e base
    rownames_to_column(var = "Plot_code")  %>% 
    full_join(mean_dist) %>%
    full_join(FD_output) %>%
    mutate(study = study_name)
  
  return(diversity_info)
  
}


cats_calculations <- function(prep_output, trait_information, trait_detail){
  
  
  input <- prep_output    
  
  if(trait_detail == "Study"){  
    
    plot_code_match <- input[["relative_abundances"]][["matrix4traits_prop"]] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "Plot_code") %>%
      pull(Plot_code) 
    
    meta_prior <- as.matrix(input[["relative_abundances"]][["matrix4traits"]]) %>%
      t() %>%
      as.data.frame() %>%
      #    mutate_if(is.numeric, ~1 * (. > 0)) %>%
      mutate(Total_cover = rowSums(.)) %>%
      select(Total_cover) %>%
      dplyr::filter(Total_cover != 0) %>%
      mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
      pull(Relative_meta_cover)
    
    states <- input[["trait_data"]] %>%
      column_to_rownames(var = "Name") %>%
      t() %>% as.data.frame()
    
    
    pb <- progress::progress_bar$new(
      format = "(Plot (n) = :Plot [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = length(plot_code_match),
      complete = "=", 
      incomplete = "-", 
      current = ">", 
      clear = FALSE, 
      width = 100
    ) 
    
    plot_n <- 1:length(plot_code_match)
    
    progress <- function(n) {
      pb$tick(tokens = list(Plot = plot_n[n]))
    }
    
    opts <- list(progress = progress)
    
    cats_output <- foreach(
      plot = plot_code_match, .combine = "rbind", .options.snow = opts,
      .errorhandling = "stop",
      .packages =
        list_of_packages
    ) %dopar% {
      
      source("src/maxent2.R")
      source("src/maxenttest2_annotated.R")
      set.seed(19970606)
      
      constraints_temp <- input[["calculations_cwm_cwv"]] %>%
        subset(Plot_code %in% plot) %>%
        select(all_of(trait_information[["cwm_names"]])) 
      
      if (!all(colnames(constraints_temp) == sort(colnames(constraints_temp)))) {
        stop("Column names (constraints) are not in alphabetical order")
      } 
      
      if (!all(rownames(states) == sort(rownames(states)))) {
        stop("Row names (states - traits) are not in alphabetical order")
      } 
      
      if (!all(colnames(states) == sort(colnames(states)))) {
        stop("Row names (states - species) are not in alphabetical order")
      } 
      
      relative_cover_temp <- input[["relative_abundances"]][["matrix4traits_prop"]] %>%
        as.data.frame() %>%
        subset(rownames(.) %in% plot) %>%
        select(all_of(colnames(states))) 
      
      # calculating model using the community-weighted trait plus a uniform prior...
      fit1 <- maxent2(constr = constraints_temp, states = states, lambda = TRUE)
      temp1 <- maxent.test2(model = fit1, obs = relative_cover_temp, nperm = 200, quick = FALSE, plot = FALSE)
      
      # Extract the lambda values
      lambda_fit1 <- list(fit1$lambda) %>%
        as.data.frame(col.names = c("Lambda_fit1")) %>%
        rownames_to_column(var = "Trait") %>%
        mutate(Trait = paste(Trait, "fit1")) %>%
        pivot_wider(names_from = "Trait", values_from = "Lambda_fit1")
      
      # an ESTIMATE (using 99 permutations) of model bias (model.mean.null.given.uniform)
      model.mean.null.given.uniform <- temp1$mean.KLR2.null
      model.bias <- model.mean.null.given.uniform
      pval.uniform <- temp1$pval
      
      # model fit using the CWM and only the uniform prior...
      model.uniform.prior.plus.traits <- temp1$KLR2.prior.plus.traits
      
      # if the model using the maximally uninformative prior plus traits is less than the
      # model bias, then correct...
      if (model.uniform.prior.plus.traits < model.bias) {
        model.uniform.prior.plus.traits <-
          model.bias
      }
      
      # calculating model using the CWM plus metacommunity prior
      fit2 <- maxent2(constr = constraints_temp, states = states, prior = meta_prior, lambda = TRUE)
      temp2 <- maxent.test2(model = fit2, obs = relative_cover_temp, nperm = 200, quick = FALSE, plot = FALSE)
      
      # Extract the lambda values
      lambda_fit2 <- list(fit2$lambda) %>%
        as.data.frame(col.names = c("Lambda_fit2")) %>%
        rownames_to_column(var = "Trait") %>%
        mutate(Trait = paste(Trait, "fit2")) %>%
        pivot_wider(names_from = "Trait", values_from = "Lambda_fit2")
      
      model.mean.null.given.prior <- temp2$mean.KLR2.null
      
      # if this null, given the metacommunity prior, is less than the mean given the maximally
      # uniformative prior, then correct...
      if (model.mean.null.given.prior < model.bias) {
        model.mean.null.given.prior <-
          model.bias
      }
      
      # fit using the CWM and the metacommunity prior...
      model.prior.plus.traits <- temp2$KLR2.prior.plus.traits
      
      # if this model, given the metacommunity prior and the traits, is less than the fit of the model using
      # the maximally uniformative prior and the traits, then correct...
      if (model.prior.plus.traits < model.uniform.prior.plus.traits) model.prior.plus.traits <- model.uniform.prior.plus.traits
      
      # OWN ADDITION. If the model given the metacommunity prior and the traits, is less than the fit of the model using the 
      #   metacommunity prior and permuted traits, then correct...
      if(model.prior.plus.traits < model.mean.null.given.prior) model.prior.plus.traits <- model.mean.null.given.prior
      
      stat_fit1 <- temp1$obs.stat
      stat_fit2 <- temp2$obs.stat
      
      pval.meta <- temp2$pval
      
      all_KLR2 <- cbind(
        model.bias, model.uniform.prior.plus.traits,
        model.mean.null.given.prior, model.prior.plus.traits, model.mean.null.given.uniform,
        lambda_fit1, lambda_fit2, pval.uniform, pval.meta, stat_fit1, stat_fit2
      ) %>%
        as.data.frame() %>%
        mutate(Plot_code = plot)
      
      gc()
      
      return(all_KLR2)
      
      
    } 
    
    
  } else {
    
    plot_code_match <- input[["relative_abundances"]][["matrix4traits_prop"]] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "Plot_code") %>%
      pull(Plot_code) 
    
    meta_prior <- as.matrix(input[["relative_abundances"]][["matrix4traits"]]) %>%
      t() %>%
      as.data.frame() %>%
      mutate(Total_cover = rowSums(.)) %>%
      select(Total_cover) %>%
      dplyr::filter(Total_cover != 0) %>%
      mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
      pull(Relative_meta_cover)
    
    pb <- progress::progress_bar$new(
      format = "(Plot (n) = :Plot [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = length(plot_code_match),
      complete = "=", 
      incomplete = "-", 
      current = ">", 
      clear = FALSE, 
      width = 100
    ) 
    
    plot_n <- 1:length(plot_code_match)
    
    progress <- function(n) {
      pb$tick(tokens = list(Plot = plot_n[n]))
    }
    
    opts <- list(progress = progress)
    
    cats_output <- foreach(
      plot = plot_code_match, .combine = "rbind", .options.snow = opts,
      .errorhandling = "stop",
      .packages =
        list_of_packages
    ) %dopar% {
      
      source("src/maxent2.R")
      source("src/maxenttest2_annotated.R")  
      set.seed(19970606)
      
      states <- input[["trait_data"]] %>%
        dplyr::filter(Plot_code == plot) %>%
        select(-Plot_code) %>%
        column_to_rownames(var = "Name") %>%
        t() %>% as.data.frame()
      
      constraints_temp <- input[["calculations_cwm_cwv"]] %>%
        dplyr::filter(Plot_code == plot) %>%
        select(all_of(trait_information[["cwm_names"]])) 
      
      relative_cover_temp <- input[["relative_abundances"]][["matrix4traits_prop"]] %>%
        as.data.frame() %>%
        subset(rownames(.) %in% plot) %>%
        select(all_of(colnames(states))) 
      
      # calculating model using the community-weighted trait plus a uniform prior...
      fit1 <- maxent2(constr = constraints_temp, states = states, lambda = TRUE)
      temp1 <- maxent.test2(model = fit1, obs = relative_cover_temp, nperm = 100, quick = FALSE, plot = FALSE)
      
      # Extract the lambda values
      lambda_fit1 <- list(fit1$lambda) %>%
        as.data.frame(col.names = c("Lambda_fit1")) %>%
        rownames_to_column(var = "Trait") %>%
        mutate(Trait = paste(Trait, "fit1")) %>%
        pivot_wider(names_from = "Trait", values_from = "Lambda_fit1")
      
      # an ESTIMATE (using 99 permutations) of model bias (model.mean.null.given.uniform)
      model.mean.null.given.uniform <- temp1$mean.KLR2.null
      model.bias <- model.mean.null.given.uniform
      pval.uniform <- temp1$pval
      
      # model fit using the CWM and only the uniform prior...
      model.uniform.prior.plus.traits <- temp1$KLR2.prior.plus.traits
      
      # if the model using the maximally uninformative prior plus traits is less than the
      # model bias, then correct...
      if (model.uniform.prior.plus.traits < model.bias) {
        model.uniform.prior.plus.traits <-
          model.bias
      }
      
      # calculating model using the CWM plus metacommunity prior
      fit2 <- maxent2(constr = constraints_temp, states = states, prior = meta_prior, lambda = TRUE)
      temp2 <- maxent.test2(model = fit2, obs = relative_cover_temp, nperm = 99, quick = FALSE, plot = FALSE)
      
      # Extract the lambda values
      lambda_fit2 <- list(fit2$lambda) %>%
        as.data.frame(col.names = c("Lambda_fit2")) %>%
        rownames_to_column(var = "Trait") %>%
        mutate(Trait = paste(Trait, "fit2")) %>%
        pivot_wider(names_from = "Trait", values_from = "Lambda_fit2")
      
      model.mean.null.given.prior <- temp2$mean.KLR2.null
      
      # if this null, given the metacommunity prior, is less than the mean given the maximally
      # uniformative prior, then correct...
      if (model.mean.null.given.prior < model.bias) {
        model.mean.null.given.prior <-
          model.bias
      }
      
      # fit using the CWM and the metacommunity prior...
      model.prior.plus.traits <- temp2$KLR2.prior.plus.traits
      
      # if this model, given the metacommunity prior and the traits, is less than the fit of the model using
      # the maximally uniformative prior and the traits, then correct...
      if (model.prior.plus.traits < model.uniform.prior.plus.traits) model.prior.plus.traits <- model.uniform.prior.plus.traits
      
      # OWN ADDITION. If the model given the metacommunity prior and the traits, is less than the fit of the model using the 
      #   metacommunity prior and permuted traits, then correct...
      if(model.prior.plus.traits < model.mean.null.given.prior) model.prior.plus.traits <- model.mean.null.given.prior
      
      
      pval.meta <- temp2$pval
      
      stat_fit1 <- temp1$obs.stat
      stat_fit2 <- temp2$obs.stat
      
      all_KLR2 <- cbind(
        model.bias, model.uniform.prior.plus.traits,
        model.mean.null.given.prior, model.prior.plus.traits, model.mean.null.given.uniform,
        lambda_fit1, lambda_fit2, pval.uniform, pval.meta, stat_fit1, stat_fit2
      ) %>%
        as.data.frame() %>%
        mutate(Plot_code = plot)
      
      gc()
      
      return(all_KLR2)
    }
    
    
    
    
  }
  
  
}




######## Unregister a cluster
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}


cats_visualisation <- function(cats_output, trait_source){
  
  cats_output_processed <- cats_output %>%
    mutate(
      just.neutral1 = model.prior.plus.traits - model.uniform.prior.plus.traits,
      just.neutral2 = model.mean.null.given.prior - model.mean.null.given.uniform,
      just.traits1 = model.uniform.prior.plus.traits - model.mean.null.given.uniform,
      just.traits2 = model.prior.plus.traits - model.mean.null.given.prior,
      information.unique.to.local.trait.constraints = just.traits2 / (1 - model.bias),
      information.unique.to.neutral.prior = just.neutral1 / (1 - model.bias),
      delta.R.trait = just.traits1,
      delta.R.neutral = just.neutral2,
      T1 = just.traits1 - just.traits2,
      T2 = just.neutral2 - just.neutral1,
      joint.information = T2 / (1 - model.bias),
      biologically.unexplained.information = (1 - model.prior.plus.traits ) / (1 - model.bias),
      check = information.unique.to.local.trait.constraints + information.unique.to.neutral.prior + joint.information +
        model.bias + biologically.unexplained.information
    ) %>%
    select(c(
      information.unique.to.local.trait.constraints, information.unique.to.neutral.prior,
      joint.information, biologically.unexplained.information, Plot_code, pval.meta, pval.uniform
    )) %>%
    pivot_longer(cols = 1:4, names_to = "Deviance", values_to = "KLR2") %>%
    mutate(Deviance = recode(Deviance,
                             biologically.unexplained.information = "Unexplained effect",
                             information.unique.to.neutral.prior = "Pure metacommunity effect",
                             joint.information = "Joint trait & metacommunity effect",
                             information.unique.to.local.trait.constraints = "Pure trait effect"
    ))
  
  cats_output_combined <- list(diversity_info, cats_output_processed) %>%
    reduce(inner_join) %>%
    mutate(traits = trait_source)
  
  
  return(cats_output_combined)
  
  
}




# Extract data from the same genus 
calculate_genusmean <- function(genus_name, traitnr){
  
  genus_mean <- try_filtered %>%
    filter(genus == genus_name) %>%
    filter(DataID %in% c(traitnr)) %>%
    select(AccSpeciesName, StdValue) %>%
    group_by(AccSpeciesName) %>%
    summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
    ungroup() %>%
    summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>%
    pull(StdValue)
}









