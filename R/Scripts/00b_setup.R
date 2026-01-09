#* 0b: Configuration Setup for DTC Metabolomics Analysis
#+ 0b.1: Load utility functions first (needed for dynamic config)
utils_path <- "R/Utilities/"
purrr::walk(
  list.files(utils_path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE),
  source
)
#+ 0b.2: Load dynamic project configuration 
config <- load_dynamic_config("auto", "All_Run/config_dynamic.yaml")
.GlobalEnv$CONFIG <- config
#+ 0b.3: Set up global paths from config 
output_path <- config$paths$output  
scripts_path <- config$paths$scripts
raw_path <- config$paths$base_data_path
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}
#+ 0b.4: Set up R environment preferences 
options(
  tibble.print_max = config$analysis$tibble_options$print_max,
  tibble.print_min = config$analysis$tibble_options$print_min,
  pillar.sigfig = config$analysis$tibble_options$sigfig
)
options(
  datatable.print.class = config$analysis$datatable_options$print_class,
  datatable.print.keys = config$analysis$datatable_options$print_keys
)
.datatable.aware = config$analysis$datatable_options$aware
#+ 0b.5: Set random seed for reproducibility
set.seed(2025)
#+ 0b.6: Import All Clinical and Path Data
#- 0b.6.1: Tumor IDs
tumor_IDs <- readr::read_csv(config$data_files$tumor_ids, show_col_types = FALSE)
#- 0b.6.1: Tumor Path
tumor_pathology_raw <- readxl::read_excel(config$data_files$tumor_pathology, sheet = "Logan Update")
#- 0b.6.1.1: Demographics
demographics_raw <- read_excel(config$paths$manifest, sheet = "Manifest")
#- 0b.6.2: Load MetaboJanitoR processed CSV files  
TFT_annot_import <- readr::read_csv(config$data_files$TFT_annot, show_col_types = FALSE)
TFT_annot_key <- readr::read_csv(config$data_files$TFT_annot_key, show_col_types = FALSE)
UFT_full_import <- readr::read_csv(config$data_files$UFT_full, show_col_types = FALSE)
UFT_filtered_import <- readr::read_csv(config$data_files$UFT_filtered, show_col_types = FALSE)
