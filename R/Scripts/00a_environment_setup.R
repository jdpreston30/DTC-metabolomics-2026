#* 0a: Environment Setup
#' Package management is handled by renv for reproducibility.
#' All required packages and their exact versions are specified in renv.lock.
#' 
#' The renv environment is verified and restored (if needed) by restore_renv.R
#' which is sourced at the beginning of run.R.

#+ 0a.1: Load conflicted and set ALL preferences BEFORE loading other packages
library(conflicted)
# Set all conflict preferences to prevent warnings during package loading
conflicts_prefer(purrr::map)
conflicts_prefer(purrr::flatten)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::count)
conflicts_prefer(dplyr::rename)
conflicts_prefer(ggplot2::margin)
conflicts_prefer(readxl::read_xlsx)
conflicts_prefer(igraph::union)
conflicts_prefer(igraph::compose)
conflicts_prefer(scales::alpha)

#+ 0a.2: Load all packages from DESCRIPTION file
source("R/Utilities/Setup/load_packages_from_description.R")
load_packages_from_description()

#+ 0a.3: Load GitHub Packages explicitly for renv detection
library(MetaboAnalystR)

#+ 0a.4: Check system dependencies
source("R/Utilities/Setup/check_system_dependencies.R")
check_system_dependencies()

