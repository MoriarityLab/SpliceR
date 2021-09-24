# run these lines one at a time, and be careful to note if there are 
# any errors put out during installation

# CRAN packages
install.packages("tidyverse")
install.packages("plyr")
install.packages("shiny")
install.packages("magrittr")
install.packages("grr")
install.packages("printr")
install.packages("rmarkdown")
install.packages("DT")
install.packages("BiocInstaller")
install.packages("httr")


Sys.setenv(LIBCURL_BUILD="winssl")
install.packages("https://github.com/jeroen/curl/archive/master.tar.gz", repos = NULL)


# Bioconductor packages
# Updated 4.7.19 due to error with Bioconductor packages
# Updated again 4.30.19

# Run this code chunk in terminal first
#
# options(repos = BiocManager::repositories())
# source("https://bioconductor.org/biocLite.R")
#
# rsconnect::deployApp("/Users/kluesner/Desktop/Research/spliceR/apps/SpliceR")

BiocManager::install("Biostrings")
