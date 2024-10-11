boot_result <- boot(temp1$values, statistic = mean_statistic, R = 5000)

boot_stat_mean_u <- mean(boot_result$t)

boot_ci <- boot.ci(boot_result, type = "basic") 
boot_ci_lower_u <- boot_ci$basic[[4]]
boot_ci_upper_u <- boot_ci$basic[[5]]

mean_statistic <- function(data, indices) {
  sample_data <- data[indices]
  return(mean(sample_data))
}


boot_stat_mean_u, boot_ci_lower_u, boot_ci_upper_u, 

boot_stat_mean_m, boot_ci_lower_m, boot_ci_upper_m, 