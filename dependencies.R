# run these lines one at a time, and be careful to note if there are 
# any errors put out during installation

# CRAN packages
install.packages("shiny")
install.packages("magrittr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("grr")
install.packages("ggplot2")
install.packages("printr")
install.packages("plyr")
install.packages("readr")
install.packages("rmarkdown")
install.packages("BiocInstaller")
install.packages("DT")

# Bioconductor packages
# Updated 1.11.18 due to error with Bioconductor packages
# Run this code chunk in terminal first

# Bioconductor packages
# Updated 4.7.19 due to error with Bioconductor packages
# Updated again 4.30.19

# Run this code chunk in terminal first
#
# options(repos = BiocManager::repositories())
# source("https://bioconductor.org/biocLite.R")
#

BiocManager::install("Biostrings")




