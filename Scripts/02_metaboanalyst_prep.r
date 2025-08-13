#* 2: Prepare data for Metaboanalyst
  #+ 2.1: Set up targeted data for metaboanalyst
    #- 2.1.1: Transform targeted dataset
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
    #- 2.1.2: Remove any columns with > 20% missing values
      feature_cols_tft <- names(TFT_metaboanalyst_i)[!names(TFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
      zero_percentages_tft <- TFT_metaboanalyst_i %>%
        select(all_of(feature_cols_tft)) %>%
        summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
        pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
      #_Identify features with <= 20% zeros (keep these)
      features_to_keep_tft <- zero_percentages_tft %>%
        filter(zero_pct <= 0.20) %>%
        pull(feature)
      #_Filter TFT_metaboanalyst to keep only good features
      TFT_metaboanalyst_raw <- TFT_metaboanalyst_i %>%
        select(Patient_ID, Variant, all_of(features_to_keep_tft))
    #- 2.1.3: Replace 0 values which remain with 1/2 the column minimum
      TFT_metaboanalyst_halfmin <- TFT_metaboanalyst_raw %>%
        mutate(across(all_of(features_to_keep_tft), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
    #- 2.1.4: Log2 transform dataset
      TFT_metaboanalyst_log2 <- TFT_metaboanalyst_halfmin %>%
        mutate(across(all_of(features_to_keep_tft), ~ log2(.x))) %>%
        mutate(Variant = as.factor(Variant))
  #+ 2.2: Set up untargeted data for metaboanalyst
    #- 2.2.1: Transform untargeted dataset
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
    #- 2.2.2: Remove any columns with > 20% missing values
      feature_cols <- names(UFT_metaboanalyst_i)[!names(UFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
      zero_percentages <- UFT_metaboanalyst_i %>%
        select(all_of(feature_cols)) %>%
        summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
        pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
      #_Identify features with <= 20% zeros (keep these)
      features_to_keep <- zero_percentages %>%
        filter(zero_pct <= 0.20) %>%
        pull(feature)
      #_Filter UFT_metaboanalyst to keep only good features
      UFT_metaboanalyst_raw <- UFT_metaboanalyst_i %>%
        select(Patient_ID, Variant, all_of(features_to_keep))
    #- 2.2.3: Replace 0 values which remain with 1/2 the column minimum
      UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
        mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
    #- 2.2.4: Log2 transform dataset
      UFT_metaboanalyst_log2 <- UFT_metaboanalyst_halfmin %>%
        mutate(across(all_of(features_to_keep), ~ log2(.x))) %>%
        mutate(Variant = as.factor(Variant))
    #- 2.2.5: Split UFT into C18 and HILIC datasets
      hilic_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "HILIC")]
      c18_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "C18")]
      UFT_HILIC_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
        select(Patient_ID, Variant, all_of(hilic_cols))
      UFT_C18_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
        select(Patient_ID, Variant, all_of(c18_cols))
    #+ 2.3: Export for metaboanalyst
      write.csv(TFT_metaboanalyst_log2, "Raw_Data/Metaboanalyst/TFT_metaboanalyst_log2.csv", row.names = FALSE)
      write.csv(UFT_metaboanalyst_raw, "Raw_Data/Metaboanalyst/UFT_metaboanalyst_raw.csv", row.names = FALSE)
      write.csv(UFT_C18_metaboanalyst_log2, "Raw_Data/Metaboanalyst/UFT_C18_metaboanalyst_log2.csv", row.names = FALSE)
      write.csv(UFT_HILIC_metaboanalyst_log2, "Raw_Data/Metaboanalyst/UFT_HILIC_metaboanalyst_log2.csv", row.names = FALSE)
      #! Performed HCA with heatmaps in metaboanalyst