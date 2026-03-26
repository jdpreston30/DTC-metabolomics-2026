#* 5: Study Design File Creation

library(dplyr)
library(stringr)
library(readr)

#+ 5.1: Import Tumor Metadata
tumor_metadata <- read.csv("R/MWB/data/tumor_pathology_table.csv") |>
  dplyr::select(Patient_ID, Sample_ID, Variant, stage_bin, Stage)
#+ 5.2: Create study design tibble
study_design <- zip_contents_mzxml |>
  dplyr::mutate(
    sample_source = if_else(
      str_starts(sample_id, regex("nist|q4", ignore_case = TRUE)),
      "Plasma",
      "Tumor"
    ),
    sample_base = str_remove(sample_id, "_\\d+$"),
    replicate   = as.integer(str_extract(sample_id, "\\d+$"))
  ) |>
  dplyr::left_join(tumor_metadata, by = c("sample_base" = "Sample_ID")) |>
  dplyr::mutate(
    Subject_ID  = Patient_ID,
    sample_type = if_else(sample_source == "Plasma", "pooled reference standard", "sample"),
    Variant     = if_else(sample_source == "Plasma", NA_character_, Variant)
  ) |>
  dplyr::select(Subject_ID, Sample_ID = sample_id, sample_source, stage_bin, RAW_FILE_NAME = file_name, Stage, Variant, sample_base, replicate, sample_type)
#+ 5.3: Write study design CSV
write_tsv(study_design, "R/MWB/data/study_design.tsv")
