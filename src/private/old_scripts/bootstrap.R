mean_statistic <- function(data, indices) {
  sample_data <- data[indices]
  return(mean(sample_data))
}
