plot_specificic <- as.matrix(input[["relative_abundances"]][["matrix4traits_prop"]]) %>%
  t() %>%
  as.data.frame() %>%
  select("TR06b") %>%
  rownames_to_column(var = "Name") %>%
  filter(Name %in%colnames(states)) %>%
  mutate(pred1 = fit1$prob,
         prior = 1/ncol(states),
         metaprior = meta_prior,
         pred2 = prob.temp[1,]) %>%
  rename(ob = "TR06b")

pred_matrix <- plot_specificic %>%
  select(ob, pred1) %>% as.matrix() %>% t()

kl_pred <- KL(pred_matrix)

1- (kl_pred/kl_prior)
tmp3 <- cats_output_processed %>% filter(Plot_code == 259888)



plot_specificic %>% 
  arrange(ob) %>%
  filter(ob > 0 ) %>%
  #mutate(pred1 = pred1/sum(pred1)) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=ob) ) +
  geom_segment( aes(x=Name ,xend=Name, y=0, yend=ob), color="grey") +
  geom_segment( aes(x=Name ,xend=Name, y=0, yend=pred1), color="grey") +
  geom_segment( aes(x=Name ,xend=Name, y=0, yend=pred2), color="grey") +
  geom_segment( aes(x=Name ,xend=Name, y=0, yend=prior), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  geom_point(aes(y=pred1), color = "#440154FF", size=3) +
  geom_point(aes(y=pred2), color = "#AA1304", size=3) +
  geom_point(aes(y=prior), color = "#334CFF", size=3) +
  coord_flip() +
  labs(x = "Species", y = "Probabilities", 
       title = "Observed (green) vs. predicted probabilities from (m,t) (red)") + 
  theme_hc() + 
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none") 




stat <- function(o, p, q) sum(o * log(p/q))

  obs.stat <- stat(obs, prob.temp, prior)
  
  sum(obs * log(prob/prior))
  
  prob_test <- prob %>% as.data.frame() 
  rownames(prob_test) <- c("prob_test")
  
  prob_test <- prob_test %>%
    as.matrix
  
  test <- rbind(obs2, prob_test, prob.temp, prior) %>% t() %>%
    as.data.frame() %>%
    rename(prior = set1, prob_perm = V3) %>%
    filter(TR06b != 0)
  
  sum(test$prob_perm)
  
  log(0.053902599 /0.0233)
  log(0.98/0.99)


test <- bergholz_veg %>%
  filter(plot == "TR07b") %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Name", values_to = "Abundance") %>%
  filter(Abundance != 0) %>%
  mutate(Rel_abund = Abundance / sum(Abundance)) 



meta_temp <- meta_prior %>%
  rownames_to_column(var = "Name")

traits_temp <- bergholz_trait %>%
  inner_join(meta_temp) %>%
  mutate(SLA_cwm = SLA * Relative_meta_cover,
         Height_cwm = Height * Relative_meta_cover,
         LDMC_cwm = LDMC * Relative_meta_cover) %>%
  mutate(SLA_cwv = Relative_meta_cover * (SLA - sum(SLA_cwm))^2,
         Height_cwv = Relative_meta_cover * (Height - sum(Height_cwm))^2,
         LDMC_cwv = Relative_meta_cover * (LDMC - sum(LDMC_cwm))^2)

library(Hmisc)

wtd.var(x = vector, weights = traits_temp$Total_cover)

vector <- traits_temp %>% pull(SLA)
weights <- traits_temp %>% pull(Relative_meta_cover)
sum(traits_temp$SLA_cwv)
sum(traits_temp$Height_cwv)
sum(traits_temp$LDMC_cwv)


ggplot(calculations_cwm_cwv, aes(x = SLA_cwv)) +
  geom_histogram() +
  geom_vline(xintercept = sum(traits_temp$SLA_cwv), size = 2, colour = "darkred" )


obs3 <- relative_cover_temp %>% t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Name") 

obs.stat <- stat(obs, prob, prior)
count <- 0

cwv_temp <- matrix(NA, 1000, 1)


while ( count < 1000) {
  # The loop will continue until one of the following conditions is met:
  # The computed statistic falls outside the confidence interval and the number of iterations exceeds 
  # or equals nperm. 
  # The number of iterations exceeds or equals minperms.

  count <- count + 1


    shuffled <- sample(1:n.states, n.states)
    states.perm <- states[, shuffled, drop = F]
    colnames(states.perm) <- s.names
    states.perm2 <- states.perm %>% 
      t() %>%
      as.data.frame() %>%
      select(LDMC) %>%
      rownames_to_column(var = "Name") %>%
      inner_join(obs3) %>%
      mutate(LDMC_cwm = LDMC * TR06b) %>%
      mutate(LDMC_cwv = TR06b * (LDMC - sum(LDMC_cwm))^2)
      
      
      
      cwv  <- sum(states.perm2$LDMC_cwv) 
      
      cwv_temp[count,1] <- cwv
      

  }

temp1$values

tmp <- cwv_temp %>% as.data.frame()

ggplot(tmp, aes(x = V1)) +
  geom_histogram() +
  geom_vline(xintercept = 0.13434682, size = 2, colour = "darkred" )

backup <- cwv_temp

val.temp <- values[1:count]
(length(backup[backup >= 0.13434682]) + 
             1)/(length(backup) + 1)


ggplot(calculations_cwm_cwv, aes(x = Height_cwv)) +
  geom_histogram() +
  geom_vline(xintercept = sum(traits_temp$Height_cwv), size = 2, colour = "darkred" )

ggplot(calculations_cwm_cwv, aes(x = LDMC_cwv)) +
  geom_histogram() +
  geom_vline(xintercept = sum(traits_temp$LDMC_cwv), size = 2, colour = "darkred" )


cwm_cal_temp <- calculations_cwm_cwv %>%
  filter(Plot_code == "TR07b")



relative_abundances <- calculate_abundance(chalmandrier_cover, chalmandrier_trait, 0.8)
trait_data <- chalmandrier_trait

cwm_cwv_all_traits <- as.data.frame(plot_code_vector) %>%
  rename(Plot_code = plot_code_vector)

#backup_traits <- cwm_cwv_all_traits

for(t in 1:length(trait_list)) {
  
  trait <- trait_list[t]
  cwm_cwv_one_trait <- data.frame()
  trait_cwm <- cwm_names[t]
  trait_cwv <- cwv_names[t]
  
  for(p in plot_code_vector){
    # Calculate 'community weighted means' (CWMs) by multiplying species proportional abundances by matched trait values and summing by site
    
    # Select one specific plot
    veg_perplot <- as.data.frame(relative_abundances$matrix4traits_prop) %>%
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
  
  cwm_cwv_all_traits <- full_join(cwm_cwv_all_traits, cwm_cwv_one_trait)
  
}







