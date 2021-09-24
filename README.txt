SpliceR version 1.2.0

Hello! Welcome to SpliceR! We hope that our program can be of use for your design of base editing sgRNAs for protein disruption.

If you are noticing any errors when using the application please feel free to reach out Mitch and Branden at klues009@umn.edu and mori0164@umn.edu for troubleshooting and input on the application.

To run this program online visit this webpage:
z.umn.edu/splicer

To run this program locally please do the following:

1) Download and install R (https://www.r-project.org/)
2) Open the dependencies.R file and run the script. Please note any errors that may occur during running.
3) Copy, paste, and run the following code chunk in your command line:

library(shiny)
library(Biostrings)
library(magrittr)
library(stringi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grr)
library(printr)
library(plyr)
library(readr)
library(printr)
library(rmarkdown)
library(DT)
library(httr)
library(curl)

3) Determine the directory of the SpliceR app on your computer. It will look something like this:

/Users/Name/Downloads/SpliceR-master

4) Copy, paste, and run the following line of code in your R command line, using the specific directory of the MultiEditR app on your computer:

shiny::runApp(“/Users/Name/Downloads/SpliceR-master”)

5) The application should launch and be used interactively
