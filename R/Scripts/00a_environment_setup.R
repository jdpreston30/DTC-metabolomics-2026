#* 0a: Environment Setup for DTC Metabolomics Analysis
#+ 0a.1: Read required packages from DESCRIPTION file 
desc_file <- "DESCRIPTION"
if (!file.exists(desc_file)) {
  stop("DESCRIPTION file not found. Please ensure you're in the project root directory.")
}
#- 0a.1.1: Read DESCRIPTION file 
desc_lines <- readLines(desc_file)
#- 0a.1.2: Extract Imports section 
imports_start <- which(grepl("^Imports:", desc_lines))
if (length(imports_start) == 0) {
  stop("No Imports section found in DESCRIPTION file.")
}
#- 0a.1.3: Find where Imports section ends (next field or end of file) 
next_field <- which(grepl("^[A-Z]", desc_lines[(imports_start + 1):length(desc_lines)])) # nolint
if (length(next_field) > 0) {
  imports_end <- imports_start + next_field[1] - 1
} else {
  imports_end <- length(desc_lines)
}
#- 0a.1.4: Extract package names 
imports_lines <- desc_lines[imports_start:imports_end]
imports_text <- paste(imports_lines, collapse = " ")
#- 0a.1.5: Clean up and extract package names 
imports_text <- gsub("Imports:", "", imports_text)
imports_text <- gsub("\\s+", " ", imports_text)
packages <- strsplit(imports_text, ",")[[1]]
required_packages <- trimws(packages)
required_packages <- required_packages[required_packages != ""]
#- 0a.1.6: Extract Bioconductor packages separately
bioc_start <- which(grepl("^Bioconductor:", desc_lines))
bioc_packages <- character(0)
if (length(bioc_start) > 0) {
  next_bioc_field <- which(grepl("^[A-Z]", desc_lines[(bioc_start + 1):length(desc_lines)]))
  if (length(next_bioc_field) > 0) {
    bioc_end <- bioc_start + next_bioc_field[1] - 1
  } else {
    bioc_end <- length(desc_lines)
  }
  bioc_lines <- desc_lines[bioc_start:bioc_end]
  bioc_text <- paste(bioc_lines, collapse = " ")
  bioc_text <- gsub("Bioconductor:", "", bioc_text)
  bioc_text <- gsub("\\s+", " ", bioc_text)
  bioc_pkgs <- strsplit(bioc_text, ",")[[1]]
  bioc_packages <- trimws(bioc_pkgs)
  bioc_packages <- bioc_packages[bioc_packages != ""]
}
#+ 0a.2: Install missing packages 
#- 0a.2.1: Check for missing CRAN packages 
missing_cran <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cran.rstudio.com/")
}
#- 0a.2.2: Check for missing Bioconductor packages
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(missing_bioc)
}
#+ 0a.3: Pre-set critical conflicts before loading packages
library(conflicted, quietly = TRUE)
conflicts_prefer(purrr::flatten)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::rename)
conflicts_prefer(ggplot2::margin)
conflicts_prefer(igraph::union)
conflicts_prefer(igraph::compose)
conflicts_prefer(purrr::map)
conflicts_prefer(readxl::read_xlsx)
conflicts_prefer(scales::alpha)
#+ 0a.4: Load all required packages 
all_packages <- c(required_packages, bioc_packages)
invisible(sapply(all_packages, library, character.only = TRUE, quietly = TRUE))
#+ 0a.5: Verify critical packages for metabolomics analysis
critical_packages <- c("dplyr", "ggplot2", "mixOmics", "vegan", "here")
missing_critical <- critical_packages[!sapply(critical_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})]
if (length(missing_critical) > 0) {
  stop("❌ Critical packages missing: ", paste(missing_critical, collapse = ", "))
}
#+ 0a.6: Set up basic R options
options(repos = c(CRAN = "https://cran.rstudio.com/"))
options(expressions = 10000)
