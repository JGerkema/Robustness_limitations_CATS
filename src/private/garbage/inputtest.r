prep_output <- prep_orig
input <- prep_orig

plot <- plot_code_match[[i]]

model <- fit1
obs <- relative_cover_temp
nperm <- 400
nperm <- 10

i <- 1
dataset <- data_input_list[[i]]

input_veg = dataset[["input_veg"]]
input_veg_fd = dataset[["input_veg_fd"]]
input_trait = dataset[["input_trait"]]  
input_trait_fd = dataset[["input_trait_fd"]]
cutoff_value = dataset[["cutoff_value"]]
trait_detail = dataset[["trait_detail"]]
trait_information = dataset[["trait_information"]]
study_name = dataset[["study_name"]]
