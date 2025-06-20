#* Dependencies
  install.packages(c("dplyr", "tidyr", "readr", "stringr"))
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
#* Import
  UFT_C18 <- read_csv("Raw_Data/C18neg_UFT_medsum.csv")
  UFT_HILIC <- read_csv("Raw_Data/HILICpos_UFT_medsum.csv")
  C18map <- read_csv("Raw_Data/Thyroid_mapping_c18neg.csv")
  HILICmap <- read_csv("Raw_Data/Thyroid_mapping_hilicpos.csv")
  tumor_IDs <- read_csv("Raw_Data/tumor_IDs.csv")
  TFT_C18 <- read_csv("Raw_Data/C18neg_TFT.csv")
  TFT_HILIC <- read_csv("Raw_Data/HILICpos_TFT.csv")  
#* Functions
  process_feature_table <- function(data_tibble, sequence_tibble, data_type) {
    data_tibble %>%
      mutate(Feature = paste(data_type, mz, time, sep = "_")) %>%
      select(Feature, everything(), -c(mz, time)) %>%
      rename_with(~ gsub("\\.mzXML$", "", .)) %>%
      pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
      pivot_wider(names_from = Feature, values_from = Value) %>%
      left_join(sequence_tibble, by = "File_Name") %>%
      select(Sample_ID, everything(), -File_Name) %>%
      mutate(Sample_ID = sub("_.*", "", Sample_ID)) %>%
      arrange(Sample_ID) %>%
      select(-Batch)
  }
  t_and_name_FT <- function(data_tibble, sequence_tibble) {
    data_tibble %>%
      pivot_longer(cols = -Feature, names_to = "File_Name", values_to = "Value") %>%
      pivot_wider(names_from = Feature, values_from = Value) %>%
      left_join(sequence_tibble, by = "File_Name") %>% # Join based on File_Name
      select(Sample_ID, everything(), -File_Name) %>% # Make Sample_ID the first column and remove File_Name
      arrange(Sample_ID)
  }
#* Run for untargeted
  #+ Name, make features, transpose
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
    UFT <- rbind(HILIC_U, C18_U)
    write.csv(UFT, "metabolomics_untargeted_FT.csv", row.names = FALSE)
#* Run for targeted
  #+ C18
    #- Step 1: Clean sample column names in C18 table
      sample_cols <- names(TFT_C18)[!(names(TFT_C18) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      cleaned_names <- gsub("\\.mzXML$", "", sample_cols)
      rename_vector <- setNames(C18map$Sample_ID, C18map$File_Name)
      new_names <- rename_vector[cleaned_names]
      names(TFT_C18)[match(sample_cols, names(TFT_C18))] <- new_names
      id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
    #- Step 2: Collapse technical replicates via median
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
    #- Step 3: Strip off anything after first underscore for q4June2022 samples
      sample_cols <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      new_sample_cols <- ifelse(
        str_starts(sample_cols, "q4June2022"),
        str_extract(sample_cols, "^[^_]+"),
        sample_cols
      )
      names(TFT_C18_median)[match(sample_cols, names(TFT_C18_median))] <- new_sample_cols
    #- Step 4: Final rename using tumor_IDs
      id_map <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
      sample_cols_median <- names(TFT_C18_median)[!(names(TFT_C18_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      new_names_median <- id_map[sample_cols_median]
      names(TFT_C18_median)[match(sample_cols_median, names(TFT_C18_median))] <- new_names_median
    #- Prep for joining
    TFT_C18_median <- TFT_C18_median %>%
      select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
      mutate(ESI = "neg") %>%
      relocate(ESI, .after = delta_ppm_vec)
  #+ HILIC
    #- Step 1: Clean sample column names in HILIC table
      sample_cols_HILIC <- names(TFT_HILIC)[!(names(TFT_HILIC) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      cleaned_names_HILIC <- gsub("\\.mzXML$", "", sample_cols_HILIC)
      rename_vector_HILIC <- setNames(HILICmap$Sample_ID, HILICmap$File_Name)
      new_names_HILIC <- rename_vector_HILIC[cleaned_names_HILIC]
      names(TFT_HILIC)[match(sample_cols_HILIC, names(TFT_HILIC))] <- new_names_HILIC
    #- Step 2: Collapse technical replicates via median
      TFT_HILIC_median <- TFT_HILIC %>%
        pivot_longer(cols = -(1:5), names_to = "Sample_ID_full", values_to = "Value") %>%
        mutate(Sample_ID = str_extract(Sample_ID_full, "^[^_]+")) %>%
        group_by(across(1:5), Sample_ID) %>%
        summarise(Value = median(Value, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = Sample_ID, values_from = Value)
    #- Step 3: Strip off anything after first underscore for q4June2022 samples
      sample_cols_HILIC_median <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      new_sample_cols_HILIC <- ifelse(
        str_starts(sample_cols_HILIC_median, "q4June2022"),
        str_extract(sample_cols_HILIC_median, "^[^_]+"),
        sample_cols_HILIC_median
      )
      names(TFT_HILIC_median)[match(sample_cols_HILIC_median, names(TFT_HILIC_median))] <- new_sample_cols_HILIC
    #- Step 4: Final rename using tumor_IDs
      id_map_HILIC <- setNames(tumor_IDs$ID, tumor_IDs$Sample_ID)
      sample_cols_HILIC_final <- names(TFT_HILIC_median)[!(names(TFT_HILIC_median) %in% c("mz", "time", "mz.1", "Name", "delta_ppm_vec"))]
      new_names_HILIC_final <- id_map_HILIC[sample_cols_HILIC_final]
      names(TFT_HILIC_median)[match(sample_cols_HILIC_final, names(TFT_HILIC_median))] <- new_names_HILIC_final
    #- Prep for joining
      TFT_HILIC_median <- TFT_HILIC_median %>%
        select(mz, time, mz.1, Name, delta_ppm_vec, everything()) %>%
        mutate(ESI = "pos") %>%
        relocate(ESI, .after = delta_ppm_vec)
  #+ Bind 
    TFT <- rbind(TFT_C18_median, TFT_HILIC_median) %>%
      mutate(mode = case_when(
        ESI == "pos" ~ "HILIC",
        ESI == "neg" ~ "C18",
        TRUE ~ NA_character_
      )) %>%
      relocate(mode, .before = ESI)
    write.csv(TFT, "metabolomics_targeted_FT.csv", row.names = FALSE)
