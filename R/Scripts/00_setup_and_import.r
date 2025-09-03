#* 0: Dependencies and setting seeds
#+ 0.0: Call all utility and modeling functions
  purrr::walk(
    list.files(
      here::here("R", "Utilities"),
      pattern = "\\.[rR]$",
      full.names = TRUE,
      recursive = TRUE
    ),
    source
  )
#+ 0.1: Dependencies
  #- 0.1.0: Define a vector of all dependencies (CRAN + Bioconductor mixOmics)
    pkgs <- c(
      "broom",
      "dplyr",
      "forcats",
      "ggplot2",
      "ggprism",
      "ggraph",
      "gtools",
      "igraph",
      "KEGGREST",
      "magick",
      "memoise",
      "mice",
      "mixOmics", # Bioconductor
      "patchwork",
      "permute",
      "pheatmap",
      "purrr",
      "RColorBrewer",
      "readr",
      "readxl",
      "rlang",
      "scales",
      "stringr",
      "tibble",
      "tidyr",
      "vegan",
      "conflicted"
    )
  #- 0.1.1: Install missing CRAN packages
    cran_pkgs <- setdiff(pkgs, rownames(installed.packages()))
    cran_pkgs <- cran_pkgs[cran_pkgs != "mixOmics"] # skip mixOmics here
    if (length(cran_pkgs)) {
      install.packages(cran_pkgs)
    }
  #- 0.1.2: Install mixOmics (Bioconductor)
    if (!requireNamespace("mixOmics", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("mixOmics")
    }
  #- 0.1.3: Load libraries
    invisible(lapply(sort(pkgs), function(pkg) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    }))
    library(grid)
#+ 0.2: Set conflicts
  conflicts_prefer(ggplot2::margin)
  conflicts_prefer(dplyr::select)
  conflicts_prefer(dplyr::filter)
  conflicts_prefer(igraph::union)
#+ 0.3: Set seed for reproducibility
  set.seed(2025)
#+ 0.5: Import Data
  #- 0.5.1: Load data
  final_rds <- file.path(raw_path, "all_objects.rds")
  #- 0.5.2: Define CSVs and corresponding names (no extension here)
  files_csv <- tibble::tribble(
    ~file, ~name,
    "C18neg_UFT_medsum", "UFT_C18",
    "HILICpos_UFT_medsum", "UFT_HILIC",
    "Thyroid_mapping_c18neg", "C18map",
    "Thyroid_mapping_hilicpos", "HILICmap",
    "tumor_IDs", "tumor_IDs",
    "C18neg_TFT", "TFT_C18",
    "HILICpos_TFT", "TFT_HILIC"
  )
  #- 0.5.3: Define XLSX files and sheets
  #- 0.5.2: Define all XLSX with sheets and corresponding import tibble names
    files_xlsx <- tibble::tribble(
      ~file, ~sheet, ~name,
      "tumor_pathology", "pathology", "tumor_pathology_raw",
      "variant_mummichog", "FVPTC_PTC_MFN", "FVPTC_PTC_MFN_raw",
      "variant_mummichog", "FTC_PTC_MFN", "FTC_PTC_MFN_raw",
      "variant_mummichog", "FVPTC_FTC_MFN", "FVPTC_FTC_MFN_raw",
      "variant_mummichog", "FVPTC_PTC_KEGG", "FVPTC_PTC_KEGG_raw",
      "variant_mummichog", "FTC_PTC_KEGG", "FTC_PTC_KEGG_raw",
      "variant_mummichog", "FVPTC_FTC_KEGG", "FVPTC_FTC_KEGG_raw",
      "cluster_mummichog", "EnrichNet_KEGG", "EnrichNet_KEGG_import",
      "cluster_mummichog", "cluster_KEGG", "cluster_EF_raw"
    )
  #- 0.5.4: Load data
    if (file.exists(final_rds)) {
      message(">>> Found final checkpoint, loading: ", final_rds)
      objs <- readRDS(final_rds)
      list2env(objs, envir = .GlobalEnv)
      rm(objs)
    } else {
      message(">>> No final checkpoint found, running import pipeline...")
      load_with_checkpoints(raw_path, files_csv, files_xlsx)
      promote_checkpoint(raw_path) # promote temporary checkpoint to final + cleanup
    }
