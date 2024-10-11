cats_calculations <- function(prep_output, trait_information, trait_detail, nperm){
  
  
  input <- prep_output    
  
  if(trait_detail == "Study"){  
    
    plot_code_match <- input[["relative_abundances"]][["matrix4traits_prop"]] %>% 
      as.data.frame() %>%
      rownames_to_column(var = "Plot_code") %>%
      pull(Plot_code) 
    
    # meta_prior <- as.matrix(input[["relative_abundances"]][["matrix4traits"]]) %>%
    #  t() %>%
    # as.data.frame() %>%
    #  mutate(Total_cover = rowSums(.)) %>%
    # select(Total_cover) %>%
    # dplyr::filter(Total_cover != 0) %>%
    #  mutate(Relative_meta_cover = Total_cover / sum(Total_cover)) %>%
    # pull(Relative_meta_cover)
    
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
      source("src/bootstrap.R")
      #conflict_prefer("filter", "dplyr")
      
      #  set.seed(19970606)
      
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
        dplyr::filter(rownames(.) %in% plot) %>% # is this last part unnecessary?
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
      values_u <-  temp1[["KLR2.null"]]
      
      # mean_statdif_u <- mean(temp1$values) - temp1$obs.stat 
      #  sd_statdif_u <- sd(temp1$values)
      ci.pval_lower_u <-  temp1$ci.pval[1]
      ci.pval_upper_u <- temp1$ci.pval[2]
      obs.stat_u <- temp1$obs.stat
      
      
      
      # This calculates the mean difference and sd between the permuted statistic and the original statistic
      
      
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
      
      values_m <-  temp2[["KLR2.null"]]
      # mean_statdif_m <- mean(temp2$values) - temp2$obs.stat 
      # sd_statdif_m <- sd(temp2$values)
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
      
      # OWN ADDITION. If the model given the metacommunity prior and the traits, is less than the fit of the model using the 
      #   metacommunity prior and permuted traits, then correct...
      # if(model.prior.plus.traits < model.mean.null.given.prior) model.prior.plus.traits <- model.mean.null.given.prior
      
      #stat_fit1 <- temp1$obs.stat
      #stat_fit2 <- temp2$obs.stat
      
      pval.meta <- temp2$pval
      
      all_KLR2 <- cbind(
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
      
      gc()
      
      all_values <- cbind(values_u, values_m) %>%
        as.data.frame() %>%
        mutate(Plot_code = plot)
      
      return(all_values)
      
      
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
      #source("src/bootstrap.R")
      # conflict_prefer("filter", "dplyr")
      
      # set.seed(19970605)
      
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
      # set.seed(1)
      # calculating model using the community-weighted trait plus a uniform prior...
      fit1 <- maxent2(constr = constraints_temp, states = states, lambda = TRUE)
      temp1 <- maxent.test2(model = fit1, obs = relative_cover_temp, nperm = nperm, quick = FALSE, plot = FALSE)
      values_u <-  temp1[["KLR2.null"]]
      
      
      # mean_statdif_u <- mean(temp1$values) - temp1$obs.stat 
      #  sd_statdif_u <- sd(temp1$values)
      ci.pval_lower_u <-  temp1$ci.pval[1]
      ci.pval_upper_u <- temp1$ci.pval[2]
      obs.stat_u <- temp1$obs.stat
      # This calculates the mean difference and sd between the permuted statistic and the original statistic
      
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
      
      values_m <-  temp2[["KLR2.null"]]
      # mean_statdif_m <- mean(temp2$values) - temp2$obs.stat 
      #  sd_statdif_m <- sd(temp2$values)
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
      
      # OWN ADDITION. If the model given the metacommunity prior and the traits, is less than the fit of the model using the 
      #   metacommunity prior and permuted traits, then correct...
      #if(model.prior.plus.traits < model.mean.null.given.prior) model.prior.plus.traits <- model.mean.null.given.prior
      
      
      pval.meta <- temp2$pval
      
      #stat_fit1 <- temp1$obs.stat
      #stat_fit2 <- temp2$obs.stat
      
      all_KLR2 <- cbind(
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
      
      gc()
      
      all_values <- cbind(values_u, values_m) %>%
        as.data.frame() %>%
        mutate(Plot_code = plot)
      
      return(all_values)
    }
    
    
    
    
  }
  
  
}

values_df <- data.frame()

for(i in 1:3){

temp <- output_list[[i]][["cats_orig"]] %>%
  as.data.frame()

values_df <- bind_rows(values_df, temp)

}
