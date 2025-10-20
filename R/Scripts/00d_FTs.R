#* 0d: Cleanup
#+ 0d.1: Targeted Processing
#- 0d.1.1: Remove QC from TFTs; tag with IDs; join with clinical metadata
TFT_annot <- process_feature_table(TFT_annot_import)
#!  Will do this later on with confirmed library for manuscript
#+ 0d.2: Untargeted Processing
#- 0d.2.1: Remove QC from UFTs (full); tag with IDs; join with clinical metadata
UFT_full <- process_feature_table(UFT_full_import)
#- 0d.2.2: Remove QC from UFTs (filtered); tag with IDs; join with clinical metadata
UFT_filtered <- process_feature_table(UFT_filtered_import)
#+ 0d.3: Define any possible feature columns in vector
feature_cols <- grep("^(C18|HILIC)", colnames(UFT_full), value = TRUE)
