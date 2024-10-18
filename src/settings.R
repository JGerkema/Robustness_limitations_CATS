#-------------------------- Packages and set up ------------------------------#
list_of_packages <- c("tidyverse", "FD", "janitor", "progress", "readxl",
                      "ggplot2", "ggthemes", "hrbrthemes", "writexl", "readr", "progress",
                      "RColorBrewer", "janitor", "ggpmisc",
                      "vegan", "foreach", "doSNOW", "docstring", "renv", 
                       "conflicted", 
                      "remotes", "roxygen2", "devtools", "geosphere",
                      "WorldFlora", "rtry", 
                      "AICcmodavg", "ggpubr", "car", "ggh4x", "ggbreak", "viridis",
                      "lmodel2", "boot", "MASS", "sf", "geodata", "raster", "patchwork", 
                      "plotbiomes")


new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(FD)
library(janitor) # Row to column names
library(progress) # Progress bar
library(readxl) # Import excel files
library(readr) # Import csv files
library(ggplot2)
library(ggthemes) # More themes for ggplot2
library(hrbrthemes) # More themes for ggplot2
library(writexl)
library(progress) # Progress bar
library(RColorBrewer) # More colors in ggplot2
library(ggpmisc) # Correlation in ggplot2
library(vegan)
library(foreach) # Package for a parallel loop
library(doSNOW)
library(roxygen2) # Needed for docstring?
library(docstring) # Provide documentation for custom build functions
library(renv)
library(conflicted) # Function conflict_prefer
library(remotes)
library(AICcmodavg) # Use the predictSE function for a mixed model
library(rstatix) # Needed for the t-test
library(ggpubr)  # Also needed for t-test and visualisation?
library(car) # Also needed for t-test and visualisation?
library(ggh4x)
library(devtools)
library(WorldFlora)
library(rtry)
library(ggbreak)
library(viridis) #Colourblind friendly
library(lmodel2) #MA regression
library(boot) #inv.logit
library(MASS) #quasibinomial linear mixed model
library(sf)
library(geodata)
library(raster)
library(patchwork)
library(plotbiomes)




# Indicating preferences
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("recode", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("load", "base")
conflict_prefer("filter", "dplyr")
conflict_prefer("where", "dplyr")



# Display numbers without scientific notation
options(scipen = 15)   

# Set the seed for reproducability
set.seed(19970606)



