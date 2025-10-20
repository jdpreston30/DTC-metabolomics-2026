#* 0d: Cleanup
#+ 0d.1: Targeted Processing
#!  Will do this later on with confirmed library for manuscript
#+ 0d.5: Set up untargeted data for analysis (metaboanalyst style)
#- 0d.5.1: Transform untargeted dataset
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
#- 0d.5.2: Remove any columns with > 20% missing values
feature_cols <- names(UFT_metaboanalyst_i)[!names(UFT_metaboanalyst_i) %in% c("Patient_ID", "Variant")]
zero_percentages <- UFT_metaboanalyst_i %>%
  select(all_of(feature_cols)) %>%
  summarise(across(everything(), ~ sum(.x == 0, na.rm = TRUE) / length(.x))) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "zero_pct")
#- 0d.5.3: Identify features with <= 20% zeros (keep these)
features_to_keep <- zero_percentages %>%
  filter(zero_pct <= 0.20) %>%
  pull(feature)
#- 0d.5.4: Filter UFT_metaboanalyst to keep only good features (this will proceed for OTHER analyses)
UFT_metaboanalyst_raw <- UFT_metaboanalyst_i %>%
  select(Patient_ID, Variant, all_of(features_to_keep))
#- 0d.5.5: Replace 0 values which remain with 1/2 the column minimum
UFT_metaboanalyst_halfmin <- UFT_metaboanalyst_raw %>%
  mutate(across(all_of(features_to_keep), ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)))
#- 0d.5.6: Log2 transform dataset
UFT_metaboanalyst_log2 <- UFT_metaboanalyst_halfmin %>%
  mutate(across(all_of(features_to_keep), ~ log2(.x))) %>%
  mutate(Variant = as.factor(Variant))
#- 0d.5.7: Split UFT into C18 and HILIC datasets
hilic_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "HILIC")]
c18_cols <- names(UFT_metaboanalyst_log2)[str_starts(names(UFT_metaboanalyst_log2), "C18")]
UFT_HILIC_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
  select(Patient_ID, Variant, all_of(hilic_cols))
UFT_C18_metaboanalyst_log2 <- UFT_metaboanalyst_log2 %>%
  select(Patient_ID, Variant, all_of(c18_cols))
#- 0d.5.8: Prepare a complete UFT for mummichog (half min and log2)
UFT_mummichog <- suppressWarnings(
  UFT_metaboanalyst_i %>%
    mutate(across(
      -c(Patient_ID, Variant),
      ~ ifelse(.x == 0, 0.5 * min(.x[.x > 0], na.rm = TRUE), .x)
    )) %>%
    mutate(across(-c(Patient_ID, Variant), ~ log2(.x))) %>%
    mutate(Variant = as.factor(Variant))
)
# ! all 0 columns will throw warning but just become NAs
#- 0c.1.2: Combine with UFT and TFT datasets
  UFT_metaboanalyst_log2_path <- UFT_metaboanalyst_log2 %>%
    left_join(tumor_pathology, by = "Patient_ID") %>%
    select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
    mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
  TFT_metaboanalyst_log2_path <- TFT_metaboanalyst_log2 %>%
    left_join(tumor_pathology, by = "Patient_ID") %>%
    select(Patient_ID, Variant, colnames(tumor_pathology), everything()) %>%
    mutate(Variant = factor(Variant, levels = c("PTC", "FV-PTC", "FTC")))
#- 0c.1.3: Define metadata and feature columns
  metadata_cols <- c("Patient_ID", "Variant", "Sex", "Age", "MFC", "LD", "LVI", "ETE", "T_computed", "T_computed_bin")
  feature_cols <- grep("^(C18|HILIC)", colnames(UFT_metaboanalyst_log2_path), value = TRUE)