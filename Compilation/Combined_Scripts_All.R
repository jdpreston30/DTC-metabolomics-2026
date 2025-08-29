# ===== File: 00_setup_and_import.r =====
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
      "vegan"
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
#+ 0.2: Set conflicts
  conflicts_prefer(ggplot2::margin)
#+ 0.3: Set seed for reproducibility
  set.seed(2025)
#+ 0.5: Import Data
  #- 0.5.1: Define paths
  raw_path <- "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Manuscripts and Projects/Active Projects/Thyroid Metabolomics/Raw Data"
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

# ===== File: 01_cleanup.r =====
#* 1: Cleanup
#+ 1.1 Untargeted Processing
  #- 1.1.1: Name, make features, transpose for C18
    C18_U <- process_feature_table(UFT_C18, C18map,"C18") %>%
      left_join(tumor_IDs, by = "Sample_ID") %>%
      select(ID, everything(),-Sample_ID) %>%
      pivot_longer(-ID, names_to = "Feature", values_to = "Value") %>%
      pivot_wider(names_from = ID, values_from = Value) %>%
      separate(Feature, into = c("mode", "mz", "rt"), sep = "_", remove = FALSE) %>%
      relocate(mode, mz, rt, .after = Feature) %>%
      mutate(ESI = case_when(
        mode == "HILIC" ~ "pos",
        mode == "C18" ~ "neg",
        TRUE ~ NA_character_)) %>%
      relocate(ESI, .after = mode)
  #- 1.1.2: Name, make features, transpose for HILIC
    HILIC_U <- process_feature_table(UFT_HILIC, HILICmap, "HILIC") %>%
      left_join(tumor_IDs, by = "Sample_ID") %>%
      select(ID, everything(), -Sample_ID) %>%
      pivot_longer(-ID, names_to = "Feature", values_to = "Value") %>%
      pivot_wider(names_from = ID, values_from = Value) %>%
      separate(Feature, into = c("mode", "mz", "rt"), sep = "_", remove = FALSE) %>%
      relocate(mode, mz, rt, .after = Feature) %>%
      mutate(ESI = case_when(
        mode == "HILIC" ~ "pos",
        mode == "C18" ~ "neg",
        TRUE ~ NA_character_)) %>%
      relocate(ESI, .after = mode)
  #- 1.1.3: Combine HILIC and C18 UFTs
    UFT <- rbind(HILIC_U, C18_U)
#+ 1.2 Targeted Processing
  #- 1.3.1: Clean sample column names in C18 table
    sample_cols <- names(TFT_C18)[!(names(TFT_C18) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    cleaned_names <- gsub("\\.mzXML$", "", sample_cols)
    rename_vector <- setNames(C18map$Sample_ID, C18map$File_Name)
    new_names <- rename_vector[cleaned_names]
    names(TFT_C18)[match(sample_cols, names(TFT_C18))] <- new_names
    id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
  #- 1.2.2: Collapse technical replicates via median
    TFT_C18_median <- TFT_C18 %>%
      # Pivot to long format to work with sample columns
      pivot_longer(cols = -(1:5), names_to = "Sample_ID_full", values_to = "Value") %>%
      # Extract base sample ID (before the underscore)
      mutate(Sample_ID = str_extract(Sample_ID_full, "^[^_]+")) %>%
      # Group and summarize by median
      group_by(across(1:5), Sample_ID) %>%
      summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
      # Pivot back to wide
      pivot_wider(names_from = Sample_ID, values_from = Value)
  #- 1.2.3: Strip off anything after first underscore for q4June2022 samples
    sample_cols <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_sample_cols <- ifelse(
      str_starts(sample_cols, "q4June2022"),
      str_extract(sample_cols, "^[^_]+"),
      sample_cols
    )
    names(TFT_C18_median)[match(sample_cols, names(TFT_C18_median))] <- new_sample_cols
  #- 1.2.4: Final rename using tumor_IDs
    id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
    sample_cols_median <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_names_median <- id_map[sample_cols_median]
    names(TFT_C18_median)[match(sample_cols_median, names(TFT_C18_median))] <- new_names_median
  #- 1.2.5: Prep for joining
    TFT_C18_median <- TFT_C18_median %>%
      select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
      mutate(ESI = "neg") %>%
      relocate(ESI, .after = delta_ppm_vec)
  #- 1.2.6: Clean sample column names in HILIC table
    sample_cols_HILIC <- names(TFT_HILIC)[!(names(TFT_HILIC) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    cleaned_names_HILIC <- gsub("\\.mzXML$", "", sample_cols_HILIC)
    rename_vector_HILIC <- setNames(HILICmap$Sample_ID, HILICmap$File_Name)
    new_names_HILIC <- rename_vector_HILIC[cleaned_names_HILIC]
    names(TFT_HILIC)[match(sample_cols_HILIC, names(TFT_HILIC))] <- new_names_HILIC
  #- 1.2.7: Collapse technical replicates via median
    TFT_HILIC_median <- TFT_HILIC %>%
      pivot_longer(cols = -(1:5), names_to = "Sample_ID_full", values_to = "Value") %>%
      mutate(Sample_ID = str_extract(Sample_ID_full, "^[^_]+")) %>%
      group_by(across(1:5), Sample_ID) %>%
      summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = Sample_ID, values_from = Value)
  #- 1.2.8: Strip off anything after first underscore for q4June2022 samples
    sample_cols_HILIC_median <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_sample_cols_HILIC <- ifelse(
      str_starts(sample_cols_HILIC_median, "q4June2022"),
      str_extract(sample_cols_HILIC_median, "^[^_]+"),
      sample_cols_HILIC_median
    )
    names(TFT_HILIC_median)[match(sample_cols_HILIC_median, names(TFT_HILIC_median))] <- new_sample_cols_HILIC
  #- 1.2.9: Final rename using tumor_IDs
    id_map_HILIC <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
    sample_cols_HILIC_final <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
    new_names_HILIC_final <- id_map_HILIC[sample_cols_HILIC_final]
    names(TFT_HILIC_median)[match(sample_cols_HILIC_final, names(TFT_HILIC_median))] <- new_names_HILIC_final
  #- 1.2.10: Prep for joining
    TFT_HILIC_median <- TFT_HILIC_median %>%
      select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
      mutate(ESI = "pos") %>%
      relocate(ESI, .after = delta_ppm_vec)
#+ 1.3: Combine targeted datasets
  TFT <- rbind(TFT_C18_median, TFT_HILIC_median) %>%
    mutate(mode = case_when(
      ESI == "pos" ~ "HILIC",
      ESI == "neg" ~ "C18",
      TRUE ~ NA_character_
    )) %>%
    relocate(mode, .before = ESI)
#+ 1.4: Prepare Targeted data for analysis (metaboanalyst style)
  #- 1.5.1: Transform targeted dataset
    TFT_metaboanalyst_i <- TFT %>%
      mutate(Feature = paste0(Name, "_", mode, "_", mz)) %>%
      select(Feature, everything(),-c(mode, ESI, mz, time, mz.1, Name, delta_ppm_vec, nist_1:q4_6)) %>%
      column_to_rownames("Feature") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Patient_ID") %>%
      mutate(Variant = case_when(
        str_detect(Patient_ID, "^FVPTC\\d+$") ~ "FV-PTC",
        str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
        str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
        TRUE ~ "Unknown"
      )) %>%
      relocate(Variant, .after = Patient_ID) %>%
      as_tibble()
  #- 1.4.2: Remove any columns with > 20% missing values
    feature_cols_tft <- names(TFT_metaboanalyst_i)[!names(TFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
    zero_percentages_tft <- TFT_metaboanalyst_i %>%
      select(all_of(feature_cols_tft)) %>%
      summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
      pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
  #- 1.5.3: Identify features with <= 20% zeros (keep these)
    features_to_keep_tft <- zero_percentages_tft %>%
      filter(zero_pct <= 0.20) %>%
      pull(feature)
  #- 1.4.4: Filter TFT_metaboanalyst to keep only good features
    TFT_metaboanalyst_raw <- TFT_metaboanalyst_i %>%
      select(Patient_ID, Variant, all_of(features_to_keep_tft))
  #- 1.4.5: Replace 0 values which remain with 1/2 the column minimum
    TFT_metaboanalyst_halfmin <- TFT_metaboanalyst_raw %>%
      mutate(across(all_of(features_to_keep_tft), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
  #- 1.4.6: Log2 transform dataset
    TFT_metaboanalyst_log2 <- TFT_metaboanalyst_halfmin %>%
      mutate(across(all_of(features_to_keep_tft), ~ log2(.x))) %>%
      mutate(Variant = as.factor(Variant))
#+ 1.5: Set up untargeted data for analysis (metaboanalyst style)
  #- 1.6.1: Transform untargeted dataset
    UFT_metaboanalyst_i <- UFT %>%
      select(-c(mode, ESI, mz, rt, nist_1:q4_6)) %>%
      column_to_rownames("Feature") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Patient_ID") %>%
      mutate(Variant = case_when(
        str_detect(Patient_ID, "^FVPTC\\d+$") ~ "FV-PTC",
        str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
        str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
        TRUE ~ "Unknown"
      )) %>%
      relocate(Variant, .after = Patient_ID) %>%
      as_tibble()
  #- 1.6.2: Remove any columns with > 20% missing values
    feature_cols <- names(UFT_metaboanalyst_i)[!names(UFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
    zero_percentages <- UFT_metaboanalyst_i %>%
      select(all_of(feature_cols)) %>%
      summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
      pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
  #- 1.6.3: Identify features with <= 20% zeros (keep these)
    features_to_keep <- zero_percentages %>%
      filter(zero_pct <= 0.20) %>%
      pull(feature)
  #- 1.6.4: Filter UFT_metaboanalyst to keep only good features (this will proceed for OTHER analyses)
    UFT_metaboanalyst_raw <- UFT_metaboanalyst_i %>%
      select(Patient_ID, Variant, all_of(features_to_keep))
  #- 1.6.5: Replace 0 values which remain with 1/2 the column minimum
    UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
      mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
  #- 1.6.6: Log2 transform dataset
    UFT_metaboanalyst_log2 <- UFT_metaboanalyst_halfmin %>%
      mutate(across(all_of(features_to_keep), ~ log2(.x))) %>%
      mutate(Variant = as.factor(Variant))
  #- 1.6.7: Split UFT into C18 and HILIC datasets
    hilic_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "HILIC")]
    c18_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "C18")]
    UFT_HILIC_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
      select(Patient_ID, Variant, all_of(hilic_cols))
    UFT_C18_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
      select(Patient_ID, Variant, all_of(c18_cols))
  #- 1.6.8: Prepare a complete UFT for mummichog (half min and log2)
    UFT_mummichog <- suppressWarnings(
      UFT_metaboanalyst_i %>%
        mutate(across(
          -c(Patient_ID, Variant),
          ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)
        )) %>%
        mutate(across(-c(Patient_ID, Variant), ~ log2(.x))) %>%
        mutate(Variant = as.factor(Variant))
    )
    #! all 0 columns will throw warning but just become NAs
#+ 1.7: Cleanup tumor pathology data
  #- 1.7.1: Read and preprocess tumor pathology data
    tumor_pathology <- tumor_pathology_raw %>%
      mutate(LVI = ifelse(LVI == 2, NA, LVI)) %>%
      assign_T_stage(ld_col = "LD", ete_col = "ETE", units = "cm", out_col = "T_stage") %>%
      mutate(
        T_computed = factor(T_stage, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
        T_computed_bin = factor(
          dplyr::case_when(
            T_computed %in% c("T1", "T2") ~ "T1-T2",
            T_computed %in% c("T3", "T4") ~ "T3-T4",
            TRUE ~ NA_character_
          ),
          levels = c("T1-T2", "T3-T4"), ordered = TRUE
        ),
        LVI = factor(
          dplyr::case_when(
            is.na(LVI) ~ NA_character_,
            LVI == 1 ~ "+LVI",
            LVI == 0 ~ "-LVI"
          ),
          levels = c("-LVI", "+LVI")
        ),
        Sex = factor(
          dplyr::case_when(
            is.na(Sex) ~ NA_character_,
            Sex == 1 ~ "Female",
            Sex == 0 ~ "Male"
          ),
          levels = c("Female", "Male")
        ),
        MFC = factor(
          dplyr::case_when(
            is.na(MFC) ~ NA_character_,
            MFC == 1 ~ "+MFC",
            MFC == 0 ~ "-MFC"
          ),
          levels = c("-MFC", "+MFC")
        ),
        Age = as.numeric(Age)
      ) %>%
      select(Patient_ID, Sex, Age, MFC, LD, LVI, ETE, T_computed, T_computed_bin, LD, LVI, T_computed)
  #- 1.7.2: Combine with UFT and TFT datasets
    UFT_metaboanalyst_log2_path <- UFT_metaboanalyst_log2 %>%
      left_join(tumor_pathology, by = "Patient_ID") %>%
      select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
      mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
    TFT_metaboanalyst_log2_path <- TFT_metaboanalyst_log2 %>%
      left_join(tumor_pathology, by = "Patient_ID") %>%
      select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
      mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
  #- 1.7.3: Define metadata and feature columns
    metadata_cols <- c("Patient_ID","Variant", "Sex","Age", "MFC", "LD", "LVI", "ETE", "T_computed", "T_computed_bin")
    feature_cols <- grep("^(C18|HILIC)", colnames(UFT_metaboanalyst_log2_path), value = TRUE)

# ===== File: 02_combined_clustering.r =====
#* 2: Combined Clustering analysis
#+ 2.1: Run Exploratory PCAs
  #- 2.1.1: Define datasets
    variant_data <- UFT_metaboanalyst_log2_path %>% 
      select(Patient_ID, Variant, all_of(feature_cols))
    T_stage_data <- UFT_metaboanalyst_log2_path %>%
      select(Patient_ID, T_computed_bin, all_of(feature_cols))
    LVI_data <- UFT_metaboanalyst_log2_path %>%
      select(Patient_ID, LVI, all_of(feature_cols))
    Sex_data <- UFT_metaboanalyst_log2_path %>%
    select(Patient_ID, Sex, all_of(feature_cols))
  #- 2.1.2: Create variant PCAs
    variant_pca_12 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 1, comp_y = 2
    )
    variant_pca_23 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 2, comp_y = 3
    )
    variant_pca_34 <- make_PCA(
      data = variant_data,
      ellipse_colors = variant_colors,
      comp_x = 3, comp_y = 4
    )
  #- 2.1.3: Create T-stage PCAs
    T_stage_PCA_12 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 1, comp_y = 2
    )
    T_stage_PCA_23 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 2, comp_y = 3
    )
    T_stage_PCA_34 <- make_PCA(
      data = T_stage_data,
      ellipse_colors = T_stage_colors,
      comp_x = 3, comp_y = 4
    )
  #- 2.1.4: Create LVI PCAs
    LVI_PCA_12 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 1, comp_y = 2
    )
    LVI_PCA_23 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 2, comp_y = 3
    )
    LVI_PCA_34 <- make_PCA(
      data = LVI_data,
      ellipse_colors = LVI_colors,
      comp_x = 3, comp_y = 4
    )
  #- 2.1.5: Create Sex PCAs
    Sex_PCA_12 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 1, comp_y = 2
    )
    Sex_PCA_23 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 2, comp_y = 3
    )
    Sex_PCA_34 <- make_PCA(
      data = Sex_data,
      ellipse_colors = sex_colors,
      comp_x = 3, comp_y = 4
    )
#+ 2.2: Create heatmaps
  #- 2.2.1: Prepare data for heatmap
    heatmap_data_T <- UFT_metaboanalyst_log2_path %>%
    select(Patient_ID, Variant, T_computed_bin, any_of(feature_cols)) %>%
    rename("T_stage" = T_computed_bin)
  #- 2.2.2: Create heatmaps with different feature selections
    variance_1000 <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "variance",
      top_features = 1000,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_bin_colors,
      cluster_colors = cluster_colors
    )
    anova_1000 <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "anova",
      top_features = 1000,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_colors_heatmap,
      cluster_colors = cluster_colors
    )
    variance_full <- make_heatmap(
      heatmap_data_T,
      variant_colors = variant_colors,
      feature_selector = "variance",
      top_features = FALSE,
      annotate_t_stage = TRUE,
      T_stage_colors = T_stage_bin_colors,
      cluster_colors = cluster_colors
    )
  #- 2.2.3: Create cluster data for PCA analysis
    UFT_with_clusters <- UFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      mutate(Cluster = factor(paste0("Cluster ", Cluster), levels = c("Cluster 1", "Cluster 2"))) %>%
      select(Patient_ID, Cluster, any_of(feature_cols)) %>%
      arrange(Cluster)
#+ 2.4: Run PERMANOVA analysis
  #- 2.4.1: Define feature columns
    permanova_features <- rownames(variance_1000$M)
  #- 2.4.1: Extract features and prepare data
    variance_study <- UFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      select(Patient_ID, Cluster, Variant, T_computed_bin, Sex, Age, MFC, LVI, any_of(permanova_features)) %>%
      rename(T_stage = T_computed_bin) %>%
      mutate(
        across(c(Sex, LVI, Variant, MFC), as.factor),
        Age = as.numeric(Age),
        T_stage = factor(T_stage, levels = c("T1-T2", "T3-T4")),
        Cluster = factor(if_else(Cluster == 1, "Cluster 1", "Cluster 2"), 
                        levels = c("Cluster 1", "Cluster 2"))
      ) %>%
      arrange(Cluster)
  #- 2.4.2: Define PERMANOVA variables
    permanova_variables <-  c("T_stage", "Sex", "LVI", "Variant", "MFC", "Age", "Cluster")
  #- 2.4.3: Impute missing values
    #! IMPUTED for LVI variable, one sample missing data
    meta_i <- variance_study %>% 
      select(all_of(permanova_variables))
    meta <- complete(mice(
      meta_i,
      m = 1,
      method = replace(setNames(rep("", ncol(meta_i)), colnames(meta_i)), "LVI", "logreg"),
      predictorMatrix = {
        p <- matrix(0, ncol(meta_i), ncol(meta_i),
          dimnames = list(colnames(meta_i), colnames(meta_i))
        )
        p["LVI", c("T_stage", "Sex", "Variant", "MFC", "Age")] <- 1
        p
      },
      seed = 123,
      printFlag = FALSE
    ), 1)
  #- 2.4.4: Extract feature data for analysis
    features_1000 <- variance_study %>%
      select(any_of(permanova_features))
  #- 2.4.5: Prepare metadata for individual tests
    meta_use <- variance_study %>%
      select(all_of(permanova_variables)) %>%
      mutate(
        across(c(T_stage, Sex, LVI, Variant, MFC, Cluster), as.factor),
        Age = as.numeric(Age)
      )
  #- 2.4.6: Run PERMANOVA for each variable
    permanova_results_1000 <- bind_rows(lapply(
      permanova_variables, 
      get_permanova,
      feature_data = features_1000,
      meta_data = meta_use,
      ctrl = permute::how(nperm = 9999),
      seed = 2025)) %>%
      arrange(p_value) %>%
      mutate(Variable = if_else(Variable == "T_stage", "T stage", Variable))
  #- 2.4.7: Create PERMANOVA visualiation data
    permanova_viz <- permanova_results_1000 %>%
      mutate(
        Significance = case_when(
          p_value < 0.001 ~ "p < 0.001",
          p_value < 0.01 ~ "p < 0.01", 
          p_value < 0.05 ~ "p < 0.05",
          p_value < 0.1 ~ "p < 0.1",
          TRUE ~ "n.s."
        ),
        Significance = factor(Significance, levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p < 0.1", "n.s.")),
        p_label = case_when(
          p_value < 0.001 ~ "p < 0.001",
          p_value < 0.05 ~ paste0("p = ", round(p_value, 3)),
          TRUE ~ paste0("p = ", round(p_value, 2))
        )
      )
#+ 2.5: Final Figure Creations
  #- 2.5.1: Create heatmap for Figure 1
    heatmap_1A <- patchwork::wrap_elements(variance_1000$heatmap_plot$gtable)
  #- 2.5.2 PCAs - Safe execution (uses unified theme)
    tryCatch(
      {
        pca_1D <- make_PCA(
          data = variant_data,
          ellipse_colors = variant_colors,
          point_size = 1
        )$plot + theme_pub_pca()

        pca_1F <- make_PCA(
          data = T_stage_data,
          ellipse_colors = T_stage_bin_colors,
          point_size = 1
        )$plot + theme_pub_pca()

        pca_1E <- make_PCA(
          data = Sex_data,
          ellipse_colors = sex_colors,
          point_size = 1
        )$plot + theme_pub_pca()

        pca_1C <- make_PCA(
          data = UFT_with_clusters,
          ellipse_colors = cluster_colors,
          point_size = 1
        )$plot + theme_pub_pca()

        cat("✅ All PCA plots created successfully!\n")
      },
      error = function(e) {
        cat("❌ Error creating PCA plots:", e$message, "\n")
        cat("Check your data and color definitions.\n")
      }
    )
  #- 2.4.8: Create PERMANOVA bar plot (vertical with horizontal bars)
    permanova_1B <- ggplot(permanova_viz, aes(y = reorder(Variable, -p_value), x = R2)) +
      geom_col(
        width = 0.72,
        fill = "black",
        color = "black",
        alpha = 1,
        linewidth = 0,
        na.rm = TRUE
      ) +
      geom_text(aes(label = p_label), hjust = -0.1, vjust = 0.5, size = 3.0, fontface = "bold") +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1)) +
      coord_cartesian(xlim = c(0, 0.7), clip = "off") +
      theme_pub_simple() +
      labs(y = NULL, x = expression(bold("R"^2 * " (Explained Variance)")))

# ===== File: 03_pathway_enrichment.r =====
#+ 3.1: Run t-tests for pairwise variant comparisons
  #- 3.1.1: FV-PTC vs PTC comparison
    variant_assignments_FVPTC_PTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "FTC") %>%
      dplyr::select(Patient_ID, Variant)
    FVPTC_vs_PTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "FTC"),
      group_assignments = variant_assignments_FVPTC_PTC,
      group_column = "Variant",
      output_filename = "FVPTC_vs_PTC.csv",
      group1_value = "PTC",
      group2_value = "FV-PTC"
    )
  #- 3.1.2: FTC vs PTC comparison
    variant_assignments_FTC_PTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "FV-PTC") %>%
      dplyr::select(Patient_ID, Variant)
    FTC_vs_PTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "FV-PTC"),
      group_assignments = variant_assignments_FTC_PTC,
      group_column = "Variant",
      output_filename = "FTC_vs_PTC.csv",
      group1_value = "PTC",
      group2_value = "FTC"
    )
  #- 3.1.3: FTC vs FV-PTC comparison
    variant_assignments_FTC_FVPTC <- UFT_mummichog %>%
      dplyr::filter(Variant != "PTC") %>%
      dplyr::select(Patient_ID, Variant)
    FTC_vs_FVPTC <- mummichog_ttests(
      data = UFT_mummichog %>% dplyr::filter(Variant != "PTC"),
      group_assignments = variant_assignments_FTC_FVPTC,
      group_column = "Variant",
      output_filename = "FTC_vs_FVPTC.csv",
      group1_value = "FV-PTC",
      group2_value = "FTC"
    )
  #! All performed in metaboanalyst, saved to spreadsheet, then imported (in 00)
#+ 3.2: Run t-tests for clusters
  #- 3.2.1: Pull cluster assignments data
    cluster_data <- variance_study %>%
      select(-any_of(feature_cols)) %>%
      left_join(tumor_pathology %>% select(Patient_ID, T_computed), by = "Patient_ID") %>%
      select(Patient_ID, Cluster)
  #- 3.2.2: Add cluster assignments and run ttests
    cluster_assignments <- UFT_mummichog %>%
      left_join(cluster_data, by = "Patient_ID") %>%
      dplyr::select(Patient_ID, Cluster)
    cluster_comparison <- mummichog_ttests(
      data = UFT_mummichog,
      group_assignments = cluster_assignments,
      group_column = "Cluster",
      output_filename = "cluster_v_cluster.csv",
      group1_value = "Cluster 1",
      group2_value = "Cluster 2"
    )
  #! All performed in metaboanalyst, saved to spreadsheet, then imported (in 00)
#+ 3.3: Pathway Enrichment Plot (MFN Variants)
  #- 3.3.1: Import results from mummichog
    FVPTC_PTC <- FVPTC_PTC_MFN_raw %>%
      mutate(
        Comparisons = "FVPTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FTC_PTC <- FTC_PTC_MFN_raw %>%
      mutate(
        Comparisons = "FTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FVPTC_FTC <- FVPTC_FTC_MFN_raw %>%
      mutate(Comparisons = "FVPTC_FTC",
             p_fisher = as.numeric(`P(Fisher)`),
             enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor)
  #- 3.3.2: Bind rows then filter to important variables
    variant_enrichment <- bind_rows(FVPTC_PTC, FTC_PTC, FVPTC_FTC) %>%
      tidyr::complete(pathway_name, Comparisons) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        Comparisons = dplyr::case_when(
          Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
          Comparisons == "FTC_PTC" ~ "PTC vs. FTC",
          Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
          TRUE ~ Comparisons
        ),
        Comparisons = factor(Comparisons, levels = c("PTC vs. FTC", "PTC vs. FV-PTC","FTC vs. FV-PTC")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 5)) %>%
      mutate(pathway_name = clean_pathway_names(pathway_name)) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = {
            # Get pathways ordered by FTC vs. FV-PTC enrichment factor
            ftc_fvptc_order <- filter(., Comparisons == "FTC vs. FV-PTC" & !is.na(enrichment_factor)) %>%
              arrange(desc(enrichment_factor)) %>%
              pull(pathway_name) %>%
              unique()
            
            # Get any additional pathways not in FTC vs. FV-PTC
            all_pathways <- unique(.$pathway_name)
            remaining_pathways <- setdiff(all_pathways, ftc_fvptc_order)
            
            # Combine: FTC vs. FV-PTC order first, then remaining pathways
            c(ftc_fvptc_order, remaining_pathways)
          }
        )
      )
  #- 3.3.3: Plot
    conflicts_prefer(ggplot2::margin)
    variant_enrichment_plot_MFN <- ggplot(
      variant_enrichment,
      aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
    ) +
      # One dummy row per facet -> avoid the warning
      geom_tile(
        data = data.frame(x = 0.5, y = 0.5),
        aes(x = x, y = y),
        width = 1, height = 1,
        fill = "white", colour = "grey80", linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_point(
        alpha = 0.95, shape = 16, stroke = 0,
        na.rm = TRUE, show.legend = TRUE
      ) +
      facet_grid(
        rows = vars(pathway_name),
        cols = vars(Comparisons),
        switch = "y", drop = FALSE
      ) +
      coord_fixed(clip = "off") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

      # Keep limits ascending; reverse legend order via guide
      scale_size_continuous(
        range = c(5, 10),
        limits = c(0, 5),
        breaks = c(5, 3, 1), # Labels show in this order…
        name = "Enrichment factor",
        guide = guide_legend(reverse = TRUE) # …because we reverse the legend
      ) +

      # Keep p limits ascending; reverse colorbar via guide
      scale_color_gradient(
        low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
        limits = c(0.01, 0.05),
        oob = scales::squish,
        name = "p-value\n",
        guide = guide_colorbar(
          reverse = TRUE, # 0.01 at top, 0.05 at bottom
          barheight = unit(5, "cm"),
          barwidth = unit(0.9, "cm")
        )
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_family = "Arial") +
      theme(
        text = element_text(family = "Arial"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(0, "pt"),
        strip.placement = "outside",
        strip.text.x.top = element_text(
          angle = 90, vjust = 0.3, hjust = 0,
          face = "bold", family = "Arial", size = 12, margin = margin(l = -14, b = 5)
        ),
        strip.text.y.left = element_text(
          angle = 0, hjust = 1,
          face = "bold", family = "Arial", size = 12,
          margin = margin(r = 6)
        ),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
      ) +
      coord_cartesian(clip = "off")
  #- 3.3.4: Save for later grob import
    ggsave(
      "Outputs/Grob/variant_enrichment_plot_MFN.png",
      variant_enrichment_plot_MFN,
      width = length(unique(variant_enrichment$Comparisons)) * 0.3 + 7,
      height = length(unique(variant_enrichment$pathway_name)) * 0.3 + 2,
      units = "in",
      dpi = 600
    )
#+ 3.4: Pathway Enrichment Plot (KEGG Variants)
  #- 3.4.1: Clean results from mummichog
    FVPTC_PTC_KEGG <- FVPTC_PTC_KEGG_raw %>%
      mutate(
        Comparisons = "FVPTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FTC_PTC_KEGG <- FTC_PTC_KEGG_raw %>%
      mutate(
        Comparisons = "FTC_PTC",
        p_fisher = as.numeric(`P(Fisher)`),
        enrichment_factor = Hits.sig / Expected) %>%
        select(Comparisons, pathway_name, p_fisher, enrichment_factor)
    FVPTC_FTC_KEGG <- FVPTC_FTC_KEGG_raw %>%
      mutate(Comparisons = "FVPTC_FTC",
             p_fisher = as.numeric(`P(Fisher)`),
             enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor)
  #- 3.4.2: Bind rows then filter to important variables
    variant_enrichment_KEGG <- bind_rows(FVPTC_PTC_KEGG, FTC_PTC_KEGG, FVPTC_FTC_KEGG) %>%
      tidyr::complete(pathway_name, Comparisons) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        Comparisons = dplyr::case_when(
          Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
          Comparisons == "FTC_PTC" ~ "PTC vs. FTC",
          Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
          TRUE ~ Comparisons
        ),
        Comparisons = factor(Comparisons, levels = c("PTC vs. FTC", "PTC vs. FV-PTC","FTC vs. FV-PTC")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 8)) %>%
      mutate(pathway_name = clean_pathway_names(pathway_name)) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = {
            # Get pathways ordered by FTC vs. FV-PTC enrichment factor
            ftc_fvptc_order <- filter(., Comparisons == "FTC vs. FV-PTC" & !is.na(enrichment_factor)) %>%
              arrange(desc(enrichment_factor)) %>%
              pull(pathway_name) %>%
              unique()
            
            # Get any additional pathways not in FTC vs. FV-PTC
            all_pathways <- unique(.$pathway_name)
            remaining_pathways <- setdiff(all_pathways, ftc_fvptc_order)
            
            # Combine: FTC vs. FV-PTC order first, then remaining pathways
            c(ftc_fvptc_order, remaining_pathways)
          }
        )
      )
  #- 3.4.3: Plot
    conflicts_prefer(ggplot2::margin)
    variant_enrichment_plot_KEGG <- ggplot(
      variant_enrichment_KEGG,
      aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
    ) +
      # One dummy row per facet -> avoid the warning
      geom_tile(
        data = data.frame(x = 0.5, y = 0.5),
        aes(x = x, y = y),
        width = 1, height = 1,
        fill = "white", colour = "grey80", linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      geom_point(
        alpha = 0.95, shape = 16, stroke = 0,
        na.rm = TRUE, show.legend = TRUE
      ) +
      facet_grid(
        rows = vars(pathway_name),
        cols = vars(Comparisons),
        switch = "y", drop = FALSE
      ) +
      coord_fixed(clip = "off") +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

      # Keep limits ascending; reverse legend order via guide
      scale_size_continuous(
        range = c(4, 8),
        limits = c(0, 8),
        breaks = c(8, 6, 4, 2), # Labels show in this order…
        name = "Enrichment factor",
        guide = guide_legend(reverse = TRUE) # …because we reverse the legend
      ) +

      # Keep p limits ascending; reverse colorbar via guide
      scale_color_gradient(
        low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
        limits = c(0.01, 0.05),
        oob = scales::squish,
        name = "p-value\n",
        guide = guide_colorbar(
          reverse = TRUE, # 0.01 at top, 0.05 at bottom
          barheight = unit(5, "cm"),
          barwidth = unit(0.9, "cm")
        )
      ) +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_family = "Arial") +
      theme(
        text = element_text(family = "Arial"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing.x = unit(0, "pt"),
        panel.spacing.y = unit(0, "pt"),
        strip.placement = "outside",
        strip.text.x.top = element_text(
          angle = 90, vjust = 0.3, hjust = 0,
          face = "bold", family = "Arial", size = 12, margin = margin(l = -14, b = 5)
        ),
        strip.text.y.left = element_text(
          angle = 0, hjust = 1,
          face = "bold", family = "Arial", size = 12,
          margin = margin(r = 6)
        ),
        legend.title = element_text(size = 12, face = "bold", family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
      ) +
      coord_cartesian(clip = "off") +
      theme(legend.position = "none")
  #- 3.4.4: Save for later grob import
    ggsave(
      "Outputs/Grob/variant_enrichment_plot_KEGG.png",
      variant_enrichment_plot_KEGG,
      width = length(unique(variant_enrichment_KEGG$Comparisons)) * 0.3 + 4.2,
      height = length(unique(variant_enrichment_KEGG$pathway_name)) * 0.3 + 2,
      units = "in",
      dpi = 600
    )
#+ 3.5: Enrichment Network Plot
  #- 3.5.1: Get enrichment factors
    cluster_enrichment_factors <- cluster_EF_raw %>%
      mutate(enrichment_factor = Hits.sig / Expected) %>%
      mutate(p_fisher = as.numeric(`P(Fisher)`)) %>%
      select(pathway_name, enrichment_factor) %>%
      arrange(desc(enrichment_factor))
  #- 3.5.2: Import significant FET pathways and 
    KEGG_enrich_raw <- EnrichNet_KEGG_import %>%
      filter(FET < 0.05) %>%
        left_join(cluster_enrichment_factors, by = "pathway_name") %>%
        mutate(neg_log_p = -log10(FET)) %>%
        select(-FET) %>%
        mutate(
          pathway_name = clean_pathway_names(pathway_name)
        ) %>%
        arrange(enrichment_factor)
  #- 3.5.3: Create Plot
    p <- build_enrichment_network(KEGG_enrich_raw, edge_thresh = 0.05, prefer_hsa = TRUE, save_path = "Outputs/Grob/enrichment_network.png", plot_title = "Cluster 1 vs. 2 Enrichment Network", width = 10, height = 10, dpi = 600, show_enrichment = FALSE)

# ===== File: 04_cluster_var_path.r =====
#* 4: Pathology analysis by cluster
#+ 4.1 Prepare Data
  cluster_path_analysis <- variance_study %>%
    select(-any_of(feature_cols)) %>%
    left_join(tumor_pathology %>% select(Patient_ID, T_computed), by = "Patient_ID") %>%
    select(T_computed, Cluster, LVI, MFC, T_stage, Variant)
#+ 4.2: Faceted Stacked Bar Chart
  #- 4.2.1 Prepare Data for Faceted Stacked Bar Plot
    stacked_bar_data_facet <- cluster_path_analysis %>%
      count(Cluster, Variant, T_computed, name = "count") %>%
      mutate(
        Cluster = factor(Cluster, levels = c("Cluster 1", "Cluster 2")),
        Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")),
        T_computed = factor(T_computed, levels = c("T4", "T3", "T2", "T1"))
      )
  #- 4.2.2 Create Faceted Stacked Bar Chart
    facet_3A <- ggplot(stacked_bar_data_facet, aes(x = Variant, y = count, fill = T_computed)) +
      geom_col(width = 0.7, color = "black", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
      scale_fill_manual(values = T_stage_cluster_colors, name = "T Stage") +
      guides(fill = guide_legend(reverse = TRUE)) +
      # add minor breaks so minor grid can draw
      scale_y_continuous(
        limits = c(0, 21), breaks = seq(0, 20, 5),
        minor_breaks = seq(0, 20, 1), expand = c(0, 0)
      ) +
      facet_grid(~Cluster, space = "free_x", scales = "free_x") +
      labs(x = NULL, y = "Patients (n)") +
      # theme_publication_bars() +
      theme(
        # don't double-draw border + axis lines
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),

        # show horizontal major + minor like your bar themes
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),

        # keep vertical grids off
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = grid::unit(c(8, 5, 5, 5), "pt"),
        axis.ticks = element_line(linewidth = 0.3),
        axis.ticks.length = grid::unit(0.2, "cm"),
        axis.title.y = element_text(size = 9, face = "bold", family = "Arial"),
        axis.title.x = element_text(size = 9, face = "bold", family = "Arial"),
        axis.text.y = element_text(size = 8, face = "bold", family = "Arial"),
        axis.text.x = element_text(size = 8, face = "bold", family = "Arial", angle = 45, hjust = 1, vjust = 0.97),
        panel.background = element_blank(),
        panel.ontop = FALSE,
        strip.text = element_text(size = 13, face = "bold", family = "Arial"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.key.size = grid::unit(0.4, "cm"),
        legend.key.width = grid::unit(0.4, "cm"),
        legend.key.height = grid::unit(0.4, "cm"),
        legend.text = element_text(size = 8, face = "bold", family = "Arial"),
        plot.title = element_text(size = 12, face = "bold", family = "Arial", hjust = 0.5),
        strip.background = element_blank()
      )

# ===== File: 05_targeted.r =====
#* 5: Targeted  
#+ 5.1: Merge TFT with path data
  #- 5.1.2: Join with TFT_metaboanalyst_log2
    TFT_metaboanalyst_log2_path_analysis <- TFT_metaboanalyst_log2_path %>%
      left_join(variance_1000$cluster_df, by = "Patient_ID") %>%
      select(Patient_ID, T_computed, Variant, Cluster, Variant, everything(), -c(Age, Sex, MFC, ETE, LVI, LD, T_computed_bin))
#+ 5.2: Statistical comparisons using targeted_metabolite_comparison function
  #- 5.2.1: ANOVA comparison between T_computed stages
    t_computed_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "T_computed",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "Cluster", "Variant"),
      exclude_isotopes = TRUE
    )
  #- 5.2.3: ANOVA comparison between Variants
    variant_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "Variant",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "T_computed", "Cluster"),
      exclude_isotopes = TRUE
    )
  #- 5.2.3: ANOVA comparison between Variants
    cluster_results <- targeted_metabolite_comparison(
      data = TFT_metaboanalyst_log2_path_analysis,
      grouping_var = "Cluster",
      test_type = "anova",
      exclude_cols = c("Patient_ID", "T_computed", "Variant"),
      exclude_isotopes = TRUE
    )
#+ 5.3: Filter metabolites and create top 10 dataset
  #- 5.3.1: Get top 10 metabolites from variant comparison (by FDR p-value)
    top_10_variant_metabolites <- variant_results %>%
      filter(significant_fdr == TRUE) %>%
      arrange(p_value_fdr) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.2: Create top 5 dataset with metadata for variants
    TFT_top5_dataset_variant <- TFT_metaboanalyst_log2_path_analysis %>%
      select(Variant, any_of(top_10_variant_metabolites)) %>%
      select(-c("Cinnamoylglycine_HILIC_206.0811588","L-Citrulline_HILIC_176.1024546", "Pyroglutamic acid_C18_128.0353782")) %>%
      rename_with(~ sub("_.*", "", .x)) %>%
      rename("Citrulline" = "L-Citrulline") %>%
      rename("Indolelactic Acid*" = "indolelactic acid") %>%
      rename("Oxoproline*" = "Oxoproline") %>%
      select(1:6)
  #- 5.3.3: Create top 10 metabolites from cluster comparison (by FDR p-value)
    top_10_cluster_metabolites <- cluster_results %>%
      filter(significant_fdr == TRUE) %>%
      arrange(p_value_fdr) %>%
      slice_head(n = 10) %>%
      pull(Metabolite)
  #- 5.3.4: Create top 5 dataset with metadata for clusters
    TFT_top5_dataset_cluster <- TFT_metaboanalyst_log2_path_analysis %>%
      select(Cluster, any_of(top_10_cluster_metabolites)) %>%
      select(-c("Niacinamide_HILIC_123.0552695")) %>%
      rename_with(~ sub("_.*", "", .x)) %>%
      rename("Nicotinamide*" = "Nicotinamide") %>%
      select(1:6)
#+ 5.4: Data prep for visualization
  #- 5.4.1: Prepare data for plotting
    #_Get metabolite columns maintaining current order
    metab_cols <- setdiff(names(TFT_top5_dataset_variant), "Variant")
    #_Convert to long format
    df_long <- TFT_top5_dataset_variant %>%
      mutate(Variant = factor(Variant, levels = names(variant_colors))) %>%
      pivot_longer(
        cols = all_of(metab_cols),
        names_to = "Metabolite",
        values_to = "Value"
      ) %>%
      mutate(Metabolite = factor(Metabolite, levels = metab_cols)) %>%
      mutate(Metabolite = forcats::fct_recode(
        Metabolite,
        "Indolelactic\nAcid*" = "Indolelactic Acid*"
      ))
    #_Compute means per metabolite × variant for bar heights
    sum_df <- df_long %>%
      group_by(Metabolite, Variant) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")
  #- 5.4.4: Cluster distribution bar plot
    #_Prepare cluster dataset for plotting (like 5.4.1)
    metab_cols_cluster <- setdiff(names(TFT_top5_dataset_cluster), "Cluster")
    df_long_cluster <- TFT_top5_dataset_cluster %>%
      mutate(Cluster = factor(paste0("Cluster ", Cluster), levels = names(cluster_colors))) %>%
      pivot_longer(
        cols = all_of(metab_cols_cluster),
        names_to = "Metabolite",
        values_to = "Value"
      ) %>%
      mutate(Metabolite = factor(Metabolite, levels = metab_cols_cluster)) %>%
      mutate(Metabolite = forcats::fct_recode(
        Metabolite,
        "Pantothenic\nAcid" = "Pantothenic acid"
      ))
    #_Compute means per metabolite × cluster for bar heights
    sum_df_cluster <- df_long_cluster %>%
      group_by(Metabolite, Cluster) %>%
      summarise(mean_value = mean(Value, na.rm = TRUE), .groups = "drop")      
#+ 5.5: Create dot plot with floating mean bars
  #- 5.5.1: Clusters based targeted top 5
    clusters_targ_3B <- ggplot() +
      # mean bars (no SD)
      geom_col(
        data = sum_df_cluster,
        aes(
          x = Metabolite, y = mean_value,
          fill = Cluster, color = Cluster, group = Cluster
        ),
        position = position_dodge2(width = 0.8, preserve = "single"),
        width = 0.72,
        linewidth = 0.5,
        na.rm = TRUE
      ) +
      # mean dots on bar centers
      geom_point(
        data = sum_df_cluster,
        aes(x = Metabolite, y = mean_value, color = Cluster, group = Cluster),
        position = position_dodge(width = 0.8),
        size = 1.2
      ) +
      # jittered raw points
      geom_jitter(
        data = df_long_cluster,
        aes(x = Metabolite, y = Value, color = Cluster),
        position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
        size = 0.5, alpha = 0.8, show.legend = FALSE
      ) +
      # palettes (use lighter fills, solid outlines)
      scale_fill_manual(values = cluster_light, drop = FALSE) +
      scale_color_manual(values = cluster_colors, drop = FALSE) +
      # keep bars alive with 0 baseline; zoom view for readability
      scale_y_pub(lims = c(0, 36), expand = c(0, 0)) +
      coord_cartesian(ylim = c(9, 36), clip = "on") +
      theme_pub_dotbar() +
      theme_pub_barplot() +
      guides(
        color = guide_legend(nrow = 1, byrow = TRUE),
        fill  = guide_legend(nrow = 1, byrow = TRUE)
      ) +
      labs(
        x = NULL,
        y = expression(bold(log[2] ~ "(Spectral Intensity)")),
        fill = "Cluster", color = "Cluster"
      ) +
      theme(
        legend.key.size   = grid::unit(0.4, "cm"),
        legend.key.width  = grid::unit(0.4, "cm"),
        legend.key.height = grid::unit(0.4, "cm"),
        legend.text       = element_text(size = 8, face = "bold", family = "Arial")
      )
  #- 5.5.2: Variant based targeted top 5
    variant_targ_3C <- ggplot() +
    geom_col(
      data = sum_df,
      aes(x = Metabolite, y = mean_value, fill = Variant, color = Variant, group = Variant),
      position = position_dodge2(width = 0.8, preserve = "single"),
      width = 0.72, linewidth = 0.5, na.rm = TRUE
    ) +
    geom_point(
      data = sum_df,
      aes(x = Metabolite, y = mean_value, color = Variant, group = Variant),
      position = position_dodge(width = 0.8),
      size = 1.2
    ) +
    geom_jitter(
      data = df_long,
      aes(x = Metabolite, y = Value, color = Variant),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
      size = 0.5, alpha = 0.8, show.legend = FALSE
    ) +
    scale_fill_manual(values = variant_light, drop = FALSE) +
    scale_color_manual(values = variant_colors, drop = FALSE) +
    scale_y_pub(lims = c(0, 36), expand = c(0, 0)) +
    coord_cartesian(ylim = c(9, 36), clip = "on") +
    theme_pub_dotbar() +
    theme_pub_barplot() +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE),
      fill  = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    labs(
      x = NULL,
      y = expression(bold(log[2] ~ "(Spectral Intensity)")),
      fill = "Variant", color = "Variant"
    ) +
    theme(
      legend.key.size   = grid::unit(0.4, "cm"),
      legend.key.width  = grid::unit(0.4, "cm"),
      legend.key.height = grid::unit(0.4, "cm"),
      legend.text       = element_text(size = 8, face = "bold", family = "Arial")
    )

# ===== File: 06_figure_creation.r =====
#* 6: Figure Creation
#+ 6.0: Setup
  #- 6.1.1: Setup a blank plot for spacing
    blank_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(tag = NULL) + # this removes the letter completely
      ggplot2::theme(
        plot.tag = ggplot2::element_text(color = "white") # fallback if patchwork still assigns tags
      )
#+ 6.1: Figure 1 - Clustering Analysis Layout
  #- 6.1.1: Fix the border widths
    permanova_1B <- permanova_1B + theme_pub_simple(border_linewidth = 0.5)
    pca_1C <- pca_1C + theme_pub_pca(border_linewidth = 0.5)
    pca_1D <- pca_1D + theme_pub_pca(border_linewidth = 0.5)
    pca_1E <- pca_1E + theme_pub_pca(border_linewidth = 0.5)
    pca_1F <- pca_1F + theme_pub_pca(border_linewidth = 0.5)
  #- 6.1.4: Assemble Figure 1 (Heatmap top 50%, vertical PERMANOVA + 4 PCAs bottom 50%)
    Figure_1 <- patchwork::wrap_plots(
      # Top 50%: Heatmap spanning full width
      heatmap_1A,
      # Bottom 50%: Vertical PERMANOVA (left) + 4 PCAs (2x2 grid, right)
      permanova_1B, pca_1C, pca_1D, pca_1E, pca_1F,
      design = "AAAA\nAAAA\nBBCD\nBBEF",  # Heatmap top 50%, vertical bar + 2x2 PCAs bottom 50%
      heights = c(1, 1, 0.5, 0.5)  # 50% heatmap, 50% bottom split
    ) +
    patchwork::plot_annotation(
      title = "Figure 1",
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
        plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
      )
    ) &
    ggplot2::theme(
      plot.tag.position = c(0, 0.98),
      plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
    )
  #- 6.1.6: Save Figure 1 to PDF and PNG
    print_to_png(Figure_1, "Figure 1")
#+ 6.2: Figure 2 - Pathway Enrichment
  #- 6.2.1: Load in KEGG and enrichment plots
    kegg_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/variant_enrichment_plot_KEGG.png")), interpolate = TRUE)
    enrichment_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/enrichment_network.png")), interpolate = TRUE)
  #- 6.2.2: Convert grob to ggplot for KEGG
    kegg_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(as.raster(magick::image_read("Outputs/Grob/variant_enrichment_plot_KEGG.png")),
          interpolate = TRUE
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "in"))
  #- 6.2.4: Convert grob to ggplot for enrichment
    enrichment_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(
        grid::rasterGrob(as.raster(magick::image_read("Outputs/Grob/enrichment_network.png")),
          interpolate = TRUE
        ),
        xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
      ) +
      ggplot2::theme_void() +
      theme(plot.margin = margin(t = -0, r = 0, b = 0, l = 0))
  #- 6.2.3: Add the fixed title band via patchwork::plot_annotation
    Figure_2 <- patchwork::wrap_plots(
      kegg_as_plot + ggplot2::labs(tag = "A"),
      enrichment_as_plot + ggplot2::labs(tag = "B"),
      blank_plot + ggplot2::labs(tag = NULL),
      blank_plot + ggplot2::labs(tag = NULL),
      design = "
      AB
      CD
      "
    ) +
      patchwork::plot_annotation(
        title = "Figure 2",
        theme = ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(
            hjust = 0, face = "bold", family = "Arial", size = 16
          ),
          plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in")
        )
      ) &
      ggplot2::theme(
        plot.tag = ggplot2::element_text(
          size = 14, face = "bold", family = "Arial", color = "black"
        ),
        plot.tag.position = c(0, 0.98)
      )
  #- 6.2.4: Print
    print_to_png(Figure_2, "Figure 2", dpi = 600)
#+ 6.3: Figure 3 - Targeted Metabolomics
  #- 6.3.1 Assemble Figure 3
    Figure_3 <- patchwork::wrap_plots(
      facet_3A, clusters_targ_3B,                # Top 25%: A, B
      variant_targ_3C,
      design = "A\nB\nC",
      heights = c(0.25, 0.25, 0.25, 0.25)
    ) +
    patchwork::plot_annotation(
      title = "Figure 3",
      tag_levels = "A",
      theme = ggplot2::theme(
        plot.title.position = "plot",
        plot.title = ggplot2::element_text(
          hjust = 0, face = "bold", family = "Arial", size = 16
        ),
        plot.margin = grid::unit(c(0.3, 2, 0.3, 2), "in") 
      )
    ) &
    ggplot2::theme(
      plot.tag.position = c(0, 0.98),
      plot.tag = ggplot2::element_text(size = 14, face = "bold", vjust = 0, hjust = 0, family = "Arial", color = "black")
    )
  #- 6.3.2 Save Figure 3 to PNG
    print_to_png(Figure_3, "Figure 3")
#+ 6.4: Supplemental Figure 1 - MFN Pathway Enrichment
  #- 6.4.1: Load in MFN plot as grob
    mfn_grob <- rasterGrob(as.raster(image_read("Outputs/Grob/variant_enrichment_plot_MFN.png")), interpolate = TRUE)
  #- 6.4.2: Convert grob to ggplot
    mfn_as_plot <- ggplot2::ggplot() +
      ggplot2::annotation_custom(mfn_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = grid::unit(c(0, 0, 0, 0), "in"))
  #- 6.4.3: Add the fixed title band via patchwork::plot_annotation
    Supplemental_Figure_1 <- mfn_as_plot +
      patchwork::plot_annotation(
        title = "Supplemental Figure 1",
        theme = ggplot2::theme(
          plot.title.position = "plot",
          plot.title = ggplot2::element_text(hjust = 0, face = "bold", family = "Arial", size = 16),
          plot.margin = grid::unit(c(0.3, 0.5, 0.3, 0.5), "in") # SAME as Fig 1 & 3
        )
      )
  #- 6.4.4: Print
    print_to_png(Supplemental_Figure_1, "Supplemental Figure 1")
#+ 6.5: Build combined PDF
{
png_files <- c(
  here::here("Figures", "Figure 1.png"),
  here::here("Figures", "Figure 2.png"),
  here::here("Figures", "Figure 3.png"),
  here::here("Figures", "Supplemental Figure 1.png")
)

# Read them in
imgs <- lapply(png_files, image_read)

# Combine into a single PDF
pdf_file <- here::here("Figures", "Figures_Compiled.pdf")
image_write(image_join(imgs), path = pdf_file, format = "pdf")

message("✅ Combined PDF saved: ", pdf_file)
}

