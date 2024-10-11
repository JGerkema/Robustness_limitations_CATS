#------------------------- Packages and set up - 2-----------------------------#

list_of_packages <- c("tidyverse", "FD", "janitor", "progress", "readxl",
                      "ggplot2", "ggthemes", "hrbrthemes", "writexl", "readr", "progress",
                      "RColorBrewer", "gcookbook", "janitor", "ggpmisc",
                      "vegan", "foreach", "doSNOW", "docstring", "renv", "ggcorrplot",
                      "reshape2", "conflicted", 
                      "remotes", "roxygen2", "devtools", "geosphere",
                      "scales", "ecodist", "WorldFlora", "rtry", "MuMIn",
                      "AICcmodavg", "ggpubr", "car", "boot", "ggh4x")
## rstatix gives problems with the foreach loop
# Display numbers without scientific notation
options(scipen = 15)   

# Set the seed to make reproducible
set.seed(19970606)

# Indicating preferences
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("recode", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("load", "base")
conflict_prefer("filter", "dplyr")
conflict_prefer("where", "dplyr")

