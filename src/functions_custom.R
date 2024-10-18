#--------------------------------------------------------#
#   Functions - Robustness and limitations of maximum 
#   entropy in plant community assembly
#
#   Jelyn Gerkema (October 2024)
#--------------------------------------------------------#

######## Set Not in 
`%notin%` <- Negate(`%in%`)


######## Calculate relative abundances
calculate_abundance <- function(input_veg_data, input_trait_data, cutoff_value) {
  #' @title Calculate the relative species abundance of each plot.
  #' @description The relative species abundance is calculated for each plot. Only sites were the species for which trait data is available 
  #' cover more than the cut-off value are retained in the dataset. Species that don't occur when the dataset is trimmed to the specified cut-off value coverage are removed. 
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
  #' @param input_veg_data The dataset for which the cwm and cwv should be calculated. Output of the calculate_abundance function.
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

######## Necessary preparations and calculations to create input-ready data for the CATS model
cats_preparation <- function(input_veg, input_trait, trait_information, cutoff_value, trait_detail){
  #' @title Necessary preparations and calculations to create input-ready data for the CATS model. 
  #' @description First, the raw abundance data is converted to relative abundances. 
  #' Only the plots with enough trait coverage are retained in the dataset using the 'calculate_abundance' function. 
  #' Then, the community weighted mean (cwm) and variance (cwv) of each plot for multiple traits is calculated using the 'calculate_cwm_cwv' function.
  #' @param input_veg The vegetation input should be a matrix where each row represents a plot and each column a species. The species columns should be alphabetically ordered. 
  #' @param input_trait The trait input should be a dataframe where each row represents a species. The species names should be stored in a column called "Name". The corresponding trait values should be stored in the rest of the columns, in alphabetical order. 
  #' @param trait_information Dataframe with the traitnames,  
  #' @param cutoff_value In some occurrences, not all species will have trait data. These species are removed. The total relative abundance of species with trait data is then calculated. If it exceeds the specified cutoff value, the plot remains in the data. Otherwise it is removed. 
  #' @param trait_detail Specifies if the traits are measured at study- or plotlevel. 
  #' @returns A list with the relative abundances, trait data with only the species actually occurring in the dataset and the calculated CWM and CWV values for each plot. 

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
        # sequentially. 
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


######## Calculate diversity information for eacht plot
calculate_diversity_info <- function(veg_data, input_veg_fd, prep_output, trait_fd,  
                                     trait_detail, study_name){
  #' @title Calculate diversity information for each plot.
  #' @description Calculates the (functional) diversity for each plot. 
  #' @param veg_data Dataframe with vegetation data. Each row should denote one plot. 
  #' Taxonomic unit (species or genus lvl) should be in columns
  #' @param input_veg_fd Dataframe with vegetation data, each row should denote one plot.
  #' Taxonomic unit (species or genus lvl) should be in columns. Should only 
  #' contain the taxonomic units with trait data. 
  #' @param prep_output Outcome of the "cats_preparation" function. 
  #' @param trait_fd Dataframe containing the trait data of each taxonomic unit.
  #' Each row should contain one taxonomic unit. Traits should be in columns.
  #' @param trait_detail Detail of the collected traits: either at plotlvl or studylvl.
  #' @param study_name Name of the study.
  #' @returns Returns a dataframe with all diversity information for each plot.
  
  input_processed <- prep_output

  all_dist <- as.matrix(vegdist(input_processed[["relative_abundances"]][["matrix4traits"]], method="bray"))
  
  plot_code_vector <- rownames(input_processed[["relative_abundances"]][["matrix4traits"]])
  
  all_dist <- all_dist %>%
    as.data.frame() 
  
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
  
 
  if(trait_detail == "Study"){
    
    trait_fd <- trait_fd %>%
      filter(rownames(.) %in% colnames(input_processed[["relative_abundances"]][["matrix4traits"]]))
    
    input_veg_fd <- input_processed[["relative_abundances"]][["matrix4traits"]]
    
    FD <- dbFD(trait_fd, input_veg_fd, stand.x = TRUE)
    
    FD_output<- data.frame(FEve = FD$FEve,
                          FDiv = FD$FDiv,
                          FRic = FD$FRic,
                          RaoQ = FD$RaoQ,
                          qual.FRic = FD$qual.FRic) %>%
      rownames_to_column(var = "Plot_code")
    
    
  
  } else {
    

    FD_output <- data.frame(matrix(nrow = length(plot_code_vector), ncol = 6))
    
    assigned_colnames <- c("FEve", "FDiv", "FRic", "RaoQ", "qual.FRic",  "Plot_code")
    colnames(FD_output) <- assigned_colnames
    m <- 0
    
    for(p in plot_code_vector){
      
      m <- m + 1
      veg_temp <- input_processed[["relative_abundances"]][["matrix4traits"]] %>%
        as.data.frame() %>%
        subset(rownames(.) %in% p) %>%
        select(where(~ any(. != 0))) # removes the columns (species) with abundance 0
      
      trait_temp <- trait_fd %>%
        filter(Plot_code %in% p) %>%
        filter(Name %in% colnames(veg_temp)) %>%
        column_to_rownames(var = "Name") %>%
        select(-c(Plot_code)) 
        
      if(ncol(veg_temp) < 3) {
        next  # If the amount of present species is less than 3, the algorithm is 
        # not able to compute the FD. Therefore, these plots are skipped. 
      }
      
      else {
        FD <- dbFD(trait_temp, veg_temp, stand.x = TRUE)
        
        
        FD_output[m, 1] <- FD$FEve
        FD_output[m, 2] <- FD$FDiv
        FD_output[m, 3] <- FD$FRic
        FD_output[m, 4] <- FD$RaoQ  
        FD_output[m, 5] <- FD$qual.FRic
        FD_output[m, 6] <- p
        
      }
    }
    

  }
  
  diversity_info <- input_processed[["relative_abundances"]][["matrix4traits_prop"]] %>% 
    as.data.frame() %>% 
    mutate(shannon = diversity(.[,1:ncol(.)], index = "shannon", MARGIN = 1, base = exp(1)),
           simpson = diversity(.[,1:ncol(.)], index = "simpson", MARGIN = 1, base =exp(1)),
           species_richness = specnumber(.[,1:ncol(.)], MARGIN = 1),
           evenness = shannon/log(species_richness),
           empty_states = ncol(input_processed[["relative_abundances"]][["matrix4traits_prop"]]) - species_richness,
           meta_size = ncol(input_processed[["relative_abundances"]][["matrix4traits_prop"]])) %>%
    select(shannon, simpson, species_richness, evenness, empty_states, meta_size) %>% # base = exp(1) gives the shannon index a log e base
    rownames_to_column(var = "Plot_code")  %>% 
    full_join(mean_dist) %>%
    full_join(FD_output) %>%
    mutate(study = study_name)
  
  return(diversity_info)
  
}


cats_calculations <- function(prep_output, trait_information, trait_detail, nperm){
  #' @title The calculation needed within the CATS framework.
  #' @description Performs the calculations of CATS.
  #' @param prep_output Outcome of the "cats_preparation" function.
  #' @param trait_information List with three vectors: the trait names, and their _cwm and _cmv extensions. 
  #' For example: Height, Heigh_cwm, Height_cwv. 
  #' @param trait_detail Detail of the collected traits: either at plotlvl or studylvl.
  #' @param nperm Number of permutations.
  #' @returns A dataframe with the outcomes of the calculations of CATS

  input <- prep_output    
 
  if(trait_detail == "Study"){  
    
    plot_code_match <- input[["relative_abundances"]][["matrix4traits_prop"]] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "Plot_code") %>%
      pull(Plot_code) 
    
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
        dplyr::filter(Plot_code %in% plot) %>%
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
        dplyr::filter(rownames(.) %in% plot) %>% 
        select(all_of(colnames(states))) 
      
      
      meta_prior <- as.matrix(input[["relative_abundances"]][["matrix4traits"]]) %>%
        t() %>%
        as.data.frame() %>%
        select(-any_of(plot)) %>%
        mutate(Total_cover = rowSums(.)) %>%
        select(Total_cover) %>%
        mutate(Total_cover = case_when(
          Total_cover == 0 ~ Total_cover + 0.001,
          TRUE ~ Total_cover)) %>%
        mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
        pull(Relative_meta_cover)
      
      
      # calculating model using the community-weighted trait plus a uniform prior...
      fit1 <- maxent2(constr = constraints_temp, states = states, lambda = TRUE)
      temp1 <- maxent.test2(model = fit1, obs = relative_cover_temp, nperm = nperm, quick = FALSE, plot = FALSE)
      
      ci.pval_lower_u <-  temp1$ci.pval[1]
      ci.pval_upper_u <- temp1$ci.pval[2]
      obs.stat_u <- temp1$obs.stat

 
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
      temp2 <- maxent.test2(model = fit2, obs = relative_cover_temp, nperm = nperm, quick = FALSE, plot = FALSE)
      
      ci.pval_lower_m <-  temp2$ci.pval[1]
      ci.pval_upper_m <- temp2$ci.pval[2]
      obs.stat_m <- temp2$obs.stat
      

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
      

      pval.meta <- temp2$pval
      
      all_KLR2 <- cbind(
        model.bias, model.uniform.prior.plus.traits,
        model.mean.null.given.prior, model.prior.plus.traits, model.mean.null.given.uniform,
        lambda_fit1, lambda_fit2, pval.uniform, pval.meta, 
        ci.pval_lower_u, ci.pval_upper_u, ci.pval_lower_m, ci.pval_upper_m,
     obs.stat_u, obs.stat_m
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
        dplyr::filter(Plot_code %in% plot) %>%
        select(-Plot_code) %>%
        column_to_rownames(var = "Name") %>%
        t() %>% as.data.frame()
      
      constraints_temp <- input[["calculations_cwm_cwv"]] %>%
        dplyr::filter(Plot_code == plot) %>%
        select(all_of(trait_information[["cwm_names"]])) 
      
      relative_cover_temp <- input[["relative_abundances"]][["matrix4traits_prop"]] %>%
        as.data.frame() %>%
        dplyr::filter(rownames(.) %in% plot) %>%
        select(all_of(colnames(states))) 
      
      meta_prior <- as.matrix(input[["relative_abundances"]][["matrix4traits"]]) %>%
        t() %>%
        as.data.frame() %>%
        select(-any_of(plot)) %>%
        mutate(Total_cover = rowSums(.)) %>%
        select(Total_cover) %>%
        mutate(Total_cover = case_when(
          Total_cover == 0 ~ Total_cover + 0.001,
          TRUE ~ Total_cover)) %>%
        mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
        pull(Relative_meta_cover)

      # calculating model using the community-weighted trait plus a uniform prior...
      fit1 <- maxent2(constr = constraints_temp, states = states, lambda = TRUE)
      temp1 <- maxent.test2(model = fit1, obs = relative_cover_temp, nperm = nperm, quick = FALSE, plot = FALSE)

    
      ci.pval_lower_u <-  temp1$ci.pval[1]
      ci.pval_upper_u <- temp1$ci.pval[2]
      obs.stat_u <- temp1$obs.stat

      
        # Extract the lambda values
      lambda_fit1 <- list(fit1$lambda) %>%
        as.data.frame(col.names = c("Lambda_fit1")) %>%
        rownames_to_column(var = "Trait") %>%
        mutate(Trait = paste(Trait, "fit1")) %>%
        pivot_wider(names_from = "Trait", values_from = "Lambda_fit1")
      
      # an ESTIMATE (using permutations) of model bias (model.mean.null.given.uniform)
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
      temp2 <- maxent.test2(model = fit2, obs = relative_cover_temp, nperm = nperm, quick = FALSE, plot = FALSE)
      
      ci.pval_lower_m <-  temp2$ci.pval[1]
      ci.pval_upper_m <- temp2$ci.pval[2]
      obs.stat_m <- temp2$obs.stat
      

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
      
      pval.meta <- temp2$pval
      
      all_KLR2 <- cbind(
        model.bias, model.uniform.prior.plus.traits,
        model.mean.null.given.prior, model.prior.plus.traits, model.mean.null.given.uniform,
        lambda_fit1, lambda_fit2, pval.uniform, pval.meta, 
        ci.pval_lower_u, ci.pval_upper_u, ci.pval_lower_m, ci.pval_upper_m,
       obs.stat_u, obs.stat_m
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

######## The complete CATS framework
complete_cats_framework <- function(input_veg, input_veg_fd, input_trait, trait_information, 
                                    trait_detail, study_name){
  #' @title The complete CATS framework, applied with local and TRY trait values, 
  #' at species and genus level. 
  #' @description Applies the CATS framework, including the necessary calculations and 
  #' preperation beforehand and extra calculations regarding (functional) diversity.
  #' @param input_veg Dataframe with vegetation data, each row should denote one plot.
  #' Species should be in columns.  
  #' @param input_veg_fd Dataframe with vegetation data, each row should denote one plot.
  #' Species are in columns. Should only contain the species with trait data. 
  #' @param input_trait Dataframe containing the trait data of each species. Each row
  #' should contain one species. Traits should be in columns.
  #' @param trait_information List with three vectors: the trait names, and their _cwm and _cmv extensions. 
  #' For example: Height, Heigh_cwm, Height_cwv. Used, among others, in the calculations of
  #' the cwv and cwm.
  #' @param trait_detail Detail of the collected traits: either at plotlvl or studylvl.
  #' @param study_name Name of the study.
  #' @returns A list with the output of the CATS framework and the preparation beforehand.
  
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
                                                     prep_output = prep_orig, 
                                                     trait_fd = input_trait_fd,  
                                                     trait_detail = trait_detail, 
                                                     study_name = study_name)
  
  diversity_info_genus <- calculate_diversity_info(veg_data = input_veg_genus,
                                                   input_veg_fd = input_veg_genus_fd,
                                                   prep_output = prep_orig_genus,
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
                                                         prep_output = prep_try, 
                                                         trait_fd = input_trait_try_fd,  
                                                         trait_detail =  "Study", 
                                                         study_name = study_name)
  
  
  diversity_info_try_genus <- calculate_diversity_info(veg_data = input_veg_genus, 
                                                       input_veg_fd = input_veg_genus_fd,
                                                       prep_output = prep_try_genus, 
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

######## Decomposition functions of CATS
cats_decomp <- function(cats_output, trait_source){
  #' @title Decomposition functions of CATS.
  #' @description Performs the decomposition functions of CATS. Final preparations 
  #' before the analyses of the outcomes and visualisation. 
  #' @param cats_output Output of the "complete_cats_framework" function. 
  #' @param trait_source Indicates the origin of the trait values, either "Local"
  #' or "TRY".
  #' @returns A dataframe containing the outcome of the decomposition functions. 
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
      biologically.unexplained.information = (1 - model.prior.plus.traits ) / (1 - model.bias)) %>%
    select(c(
      information.unique.to.local.trait.constraints, information.unique.to.neutral.prior,
      joint.information, biologically.unexplained.information, Plot_code, pval.meta, pval.uniform
    )) %>%
    pivot_longer(cols = 1:4, names_to = "Deviance", values_to = "KLR2") %>%
    mutate(Deviance = recode(Deviance,
                             biologically.unexplained.information = "Unexplained",
                             information.unique.to.neutral.prior = "Dispersal mass effect",
                             joint.information = "Joint trait & dispersal mass effect",
                             information.unique.to.local.trait.constraints = "Trait-based filtering"
    ))
  
  cats_output_combined <- list(diversity_info, cats_output_processed) %>%
    reduce(inner_join) %>%
    mutate(traits = trait_source)
  
  
  return(cats_output_combined)
  
  
}

######## Visualization of the raw model output from CATS
plot_cats4_function <- function(y_var, x_var, 
                                 x_lab = NULL, y_lab = NULL, 
                                 figure_table){
  #' @title Visualization of the raw model output from CATS
  #' @description Visualization of the raw model output from CATS (before the 
  #' decomposition functions). Data is analyzed using major axis regressions. 
  #' @param y_var y variable.
  #' @param x_var x variable.
  #' @param x_lab Title for the x-axis. Default is NULL.
  #' @param y_lab Title for the y-axis. Default is NULL.
  #' @param figure_table Specifies whether it returns the figure or the table. 
  #' @returns Either a figure or the corresponding table with the statistical information.
  cats <- raw_plot %>%
    filter(traits == y_var) %>%
    filter(Model_type %notin% c("trait_minus_bias", "raw_meta_trait_uncor")) %>%
    select(-c(FEve, FRic, qual.FRic, RaoQ, traits,
              obs.stat_u)) %>%
    distinct() %>%
    select(fit, pval.meta, study_group, plot_study, Model_type) 
  
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
  
  model_types <- c("Estimate of model bias (u)", "Metacommunity (m)",
                   "Traits (u,t)",  "Metacommunity and traits (m,t)")
  
  cats_both_models <- model_summary_all <- data.frame()
  
  for(i in 1:length(model_types) ){
    
    model_type_select <- model_types[i]
    
    cats_one_model_type <- cats_both %>%
      filter(Model_type %in% model_type_select) %>%
      drop_na(fit)
    
    
    rma_model <- lmodel2(fit ~ fit_xvar, data = cats_one_model_type,
                         "interval", "interval", 95)
    
    
    CI <- as.data.frame(rma_model$confidence.intervals) %>%
      filter(Method == "MA") %>%
      select(-`2.5%-Intercept`, -`97.5%-Intercept`)

    model_summary_df <- as.data.frame(rma_model$regression.results) %>%
      filter(Method == "MA") %>%
      mutate(Model_type = model_type_select) %>%
      relocate(Model_type) %>%
      full_join(CI) %>%
      select(-Method, -`Angle (degrees)`, -`P-perm (1-tailed)`) %>%
      mutate(`R²` = rma_model$rsquare)
    
    cats_both_models <- bind_rows(cats_both_models, cats_one_model_type)
    model_summary_all <- bind_rows(model_summary_all, model_summary_df)
  }
  
  plot_count <- cats_both_models %>% 
    group_by(study_group) %>% 
    tally() %>%
    mutate(n = n/4)
  
  plot <-  ggplot(cats_both, aes(fit_xvar, fit)) +
    geom_point(alpha = 0.15, aes(colour = study_group), size = 4) +
    theme_classic() +
    labs(title = NULL,  
         y = y_lab, x = x_lab) +
    labs(col = "Dataset")+
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = 2) +
    stat_ma_line(method = "MA",
                 range.y = "interval", range.x = "interval",  colour = "gray20",
                 nperm = 95) +
    facet_wrap(~Model_type) +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 1)) +
    guides(color = guide_legend(nrow = 4)) + 
    theme(axis.text = element_text(size = 13),
          legend.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          plot.title = element_text(size=15),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14)) +
    scale_colour_viridis_d(option = "H", 
                           labels = paste0(plot_count$study_group, " (", plot_count$n, ")")) 
  
  
  if(figure_table == "figure"){
    return(plot)}else{
      return(model_summary_all)}
}

