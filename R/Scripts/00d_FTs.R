#* 0d: Cleanup
#+ 0d.1: Targeted Processing
#- 0d.1.1: Remove QC from TFTs; tag with IDs; join with clinical metadata
TFT_annot <- process_feature_table(TFT_annot_import) |>
  select(-c(Variant:Stage)) |>
  remove_qc(threshold = 0.2)
#+ 0d.2: Untargeted Processing
#- 0d.2.1: Remove QC from UFTs (full); tag with IDs; join with clinical metadata
UFT_full <- process_feature_table(UFT_full_import) |>
  select(-c(Variant:Stage))
#- 0d.2.2: Remove QC from UFTs (filtered); tag with IDs; join with clinical metadata; apply remove_qc)
UFT_filtered <- process_feature_table(UFT_filtered_import) |>
  select(-c(Variant:Stage)) |>
  # remove features with >20% identical values for QC (half min replace)
  remove_qc(threshold = 0.2)
#+ 0d.3: Define any possible feature columns in vector
feature_cols <- grep("^(C18|HILIC)", colnames(UFT_full), value = TRUE)