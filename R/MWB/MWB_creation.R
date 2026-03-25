#* MWB Creation

#+ Subset mzXML
#- Pull metadata from UFT
metadata_relevant <- UFT |>
  select(Patient:any_PGD)
#- Subset to standards/controls and study samples; join metadata
  mzxml_filtered <- zip_contents_mzxml |>
    filter(
      str_starts(sample_id, "nist") | 
      str_starts(sample_id, "q") | 
      (str_starts(sample_id, "H") & str_detect(sample_id, "S0_\\d+$"))
    ) |>
    mutate(sample_id = str_replace(sample_id, "H46SS0", "H46S0")) |>
    mutate(Patient = case_when(
      str_starts(sample_id, "H") ~ str_extract(sample_id, "H\\d+"),
      TRUE ~ NA_character_
    )) |>
    mutate(order = str_extract(file_name_base, "\\d+$")) |>
    mutate(replicate = str_extract(sample_id, "\\d+$")) |>
    mutate(sample_type = case_when(
      str_starts(sample_id, "q") | str_starts(sample_id, "nist") ~ "pooled reference standard",
      str_starts(sample_id, "H") ~ "sample",
      TRUE ~ NA_character_
    )) |>
    mutate(sample_base = str_remove(sample_id, "_\\d+$")) |> 
    arrange(Batch, order) |>
    left_join(metadata_relevant, by = "Patient")
#- Create c18 and hilic versions
mzxml_filtered_c18 <- mzxml_filtered |>
  filter(str_detect(folder, "c18neg"))
mzxml_filtered_hilicpos <- mzxml_filtered |>
  filter(str_detect(folder, "hilicpos"))
#- Write outputs
write.csv(mzxml_filtered_c18, "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/raw_zips/subsets/AJT/mzxml_filtered_AJT_c18.csv")
write.csv(mzxml_filtered_hilicpos, "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/raw_zips/subsets/AJT/mzxml_filtered_AJT_hilicpos.csv")
#+ Create MWB tsv metadata
#- Convert into MWB format
mzxml_MWB <- mzxml_filtered |>
  mutate(sample_source = "Plasma") |>
  select(Subject_ID = Patient, Sample_ID = sample_id, sample_source, severe_PGD, Batch, RAW_FILE_NAME = file_name, PGD_grade_tier, sample_base, replicate, sample_type)
#- Write output
write_tsv(mzxml_MWB, "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/LCMS/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/raw_zips/MWB/AJT/mzxml_MWB_AJT.txt")
#+ Create MWB TFT
#- Import raw feature tables
c18_raw <- read_tsv("/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/LCMS/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/xMSanalyzer/c18neg/Stage3b/RAW_mzcalibrated_untargeted_featuretable.txt")
hilic_raw <- read_tsv("/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/LCMS/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/xMSanalyzer/hilicpos/Stage3b/RAW_mzcalibrated_untargeted_featuretable.txt")
#- Create renaming vectors from RAW_FILE_NAME to Sample_ID
rename_vector <- setNames(mzxml_MWB$Sample_ID, mzxml_MWB$RAW_FILE_NAME)
#- Mutate c18 into MWB format; subset; pull sample names
c18_mwb <- c18_raw |>
  mutate(Feature = paste0(mz, "_", time)) |>
  select(Feature, any_of(mzxml_MWB$RAW_FILE_NAME)) |>
  rename_with(~ rename_vector[.x], .cols = -Feature)
#- Mutate hilic into MWB format; subset; pull sample names
hilic_mwb <- hilic_raw |>
  mutate(Feature = paste0(mz, "_", time)) |>
  select(Feature, any_of(mzxml_MWB$RAW_FILE_NAME)) |>
  rename_with(~ rename_vector[.x], .cols = -Feature)
#- Write outputs
write_tsv(c18_mwb, "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/LCMS/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/raw_zips/MWB/AJT/UFTs/c18neg_uft.txt")
write_tsv(hilic_mwb, "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/MS_raw_data/LCMS/Chan_studies_MS_raw/PGD_studies/Chan_PGD_9_2025/raw_zips/MWB/AJT/UFTs/hilicpos_uft.txt")
