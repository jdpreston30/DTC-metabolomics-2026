#* 5: Study Design File Creation

library(dplyr)
library(stringr)
library(readr)

#! DE-IDENTIFIED. Sample_ID and sample_base carry Patient_ID (F1/P1/FVPTC2), never the specimen
#! accession from the bank manifest. Those accessions run sequentially -- across this cohort they
#! correlate with collection year at Spearman r ~= 0.95 -- so they are not random codes and must not
#! reach a public repository or deposition. The accession survives only as the join key inside this
#! script; it is never written out. Do not reintroduce it, and do not cite examples of it in comments.
#+ 5.1: Import Tumor Metadata
tumor_metadata <- read.csv("R/MWB/data/tumor_pathology_table.csv") |>
  dplyr::select(Patient_ID, Sample_ID, Variant, stage_bin, Stage, Age, Sex)
#+ 5.2: Create study design tibble
study_design <- zip_contents_mzxml |>
  dplyr::mutate(
    sample_source = if_else(
      str_starts(sample_id, regex("nist|q4", ignore_case = TRUE)),
      "Plasma",
      "Tumor"
    ),
    accession   = str_remove(sample_id, "_\\d+$"),   # join key only -- not exported
    replicate   = as.integer(str_extract(sample_id, "\\d+$"))
  ) |>
  dplyr::left_join(tumor_metadata, by = c("accession" = "Sample_ID")) |>
  dplyr::mutate(
    Subject_ID  = Patient_ID,
    sample_type = if_else(sample_source == "Plasma", "pooled reference standard", "sample"),
    Variant     = if_else(sample_source == "Plasma", NA_character_, Variant),
#! Plasma standards have no accession, so their sample_id is already de-identified and is kept as-is.
    sample_base = if_else(sample_source == "Plasma", accession, Patient_ID),
    Sample_ID   = paste0(sample_base, "_", replicate)
  ) |>
  dplyr::select(Subject_ID, Sample_ID, sample_source, stage_bin, RAW_FILE_NAME = file_name,
                Stage, Variant, Age, Sex, sample_base, replicate, sample_type)
#+ 5.3: Guard -- no accession may survive into the output, and the file mapping must stay 1:1
stopifnot(
  !any(study_design$Sample_ID   %in% tumor_metadata$Sample_ID),
  !any(study_design$sample_base %in% tumor_metadata$Sample_ID),
  !any(duplicated(study_design$Sample_ID)),
  !any(duplicated(study_design$RAW_FILE_NAME)),
  !any(is.na(study_design$Sample_ID)),
  all(study_design$sample_source != "Tumor" | !is.na(study_design$Age))
)
#+ 5.4: Write study design CSV
write_tsv(study_design, "R/MWB/data/study_design.tsv")
