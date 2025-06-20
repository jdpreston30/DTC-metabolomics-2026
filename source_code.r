#* Dependencies
  install.packages(c("dplyr", "tidyr", "readr", "stringr"))
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
#* Import
  TFT <- read_csv("metabolomics_targeted_FT.csv")
  UFT <- read_csv("metabolomics_untargeted_FT.csv")
  tumor_pathology <- read_xlsx("tumor_pathology.xlsx")