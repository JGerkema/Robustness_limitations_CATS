install.packages("NCmisc")
library(NCmisc)
library(grateful)

c("FD", "tidyverse", "FD", "janitor", "progress", "readxl", "readr", "ggplot2", 
  "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
  "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
  "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest") %>%
  map(citation) %>%
  print(style = "bibtex") 

c("ggbreak") %>%
map(citation) %>%
  print(style = "bibtex")

get_pkgs_info(pkgs = c("FD", "tidyverse", "janitor", "progress", "readxl", "readr", "ggplot2", 
                       "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
                       "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
                       "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest", "ggh4x",
                       "viridis", "sf", "geodata", "raster", "patchwork", "plotbiomes", "lmodel2",
                       "rnaturalearth"), 
              out.dir = getwd(), out.format = "docx")




cite_packages(pkgs = c("ggplot2", "FD", "tidyverse", "janitor", "progress", "readxl", "readr",  
                       "ggthemes", "hrbrthemes", "writexl", "RColorBrewer", "vegan", "foreach", 
                       "doSNOW", "roxygen2", "docstring", "renv", "conflicted", "AICcmodavg", "rstatix", 
                       "ggpubr", "car", "devtools", "taxonlookup", "WorldFlora", "rtry", "lmerTest", "ggh4x"),
              out.format = "docx", out.dir = ".")

c("sf", "geodata", "raster", "patchwork", "plotbiomes", "rnaturalearth")%>%
  map(citation) %>%
  print(style = "bibtex")

cite_packages(pkgs = c("sf", "geodata", "raster", "patchwork", "plotbiomes"),
              out.format = "docx", out.dir = ".")


packageVersion("rnaturalearth")


get_pkgs_info(pkgs = c("tidyverse", "FD", "janitor", "progress", "readxl",
                       "ggplot2", "ggthemes", "hrbrthemes", "writexl", "readr", 
                       "RColorBrewer",  "ggpmisc",
                       "vegan", "foreach", "doSNOW", "docstring", "renv", 
                       "conflicted", "remotes", "roxygen2", "devtools", "geosphere",
                       "WorldFlora", "rtry", "AICcmodavg", "ggpubr", "car", "ggh4x", "ggbreak", "viridis",
                       "lmodel2", "boot", "MASS", "sf", "geodata", "raster", "patchwork", 
                       "plotbiomes"), 
              out.dir = getwd(), out.format = "docx")

cite_packages(pkgs = c("tidyverse", "FD", "janitor", "progress", "readxl",
                       "ggplot2", "ggthemes", "hrbrthemes", "writexl", "readr", 
                       "RColorBrewer", "ggpmisc",
                       "vegan", "foreach", "doSNOW", "docstring", "renv", 
                       "conflicted", "remotes", "roxygen2", "devtools", "geosphere",
                       "WorldFlora", "rtry", "AICcmodavg", "ggpubr", "car", "ggh4x", "ggbreak", "viridis",
                       "lmodel2", "boot", "MASS", "sf", "geodata", "raster", "patchwork", 
                       "plotbiomes"),
              out.format = "docx", out.dir = ".")


cite_packages(pkgs = c("boot", "MASS", "viridis", "lmodel2"),
              out.format = "docx", out.dir = ".")

