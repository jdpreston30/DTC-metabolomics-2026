#* 7: Feature Table Import

library(yaml)
library(readr)
library(dplyr)

#+ 7.1: Load config
#- 7.1.1: Read config YAML
config <- yaml::read_yaml("R/MWB/config.yaml")

#+ 7.2: Import feature tables and study design
#- 7.2.1: Import C18 negative feature table
ft_c18neg    <- readr::read_tsv(config$c18neg_feature_table)
#- 7.2.2: Import HILIC positive feature table
ft_hilicpos  <- readr::read_tsv(config$hilicpos_feature_table)
#- 7.2.3: Import study design
study_design <- readr::read_tsv("R/MWB/data/study_design.tsv")

#+ 7.3: Clean feature tables and map sample names
#- 7.3.1: Build RAW_FILE_NAME -> Sample_ID lookup
sample_map <- setNames(study_design$Sample_ID, study_design$RAW_FILE_NAME)
#- 7.3.2: Clean C18 negative table and rename sample columns
ft_c18neg_clean <- ft_c18neg |>
  dplyr::mutate(Feature = paste(mz, time, sep = "_")) |>
  dplyr::select(Feature, dplyr::ends_with(".mzXML")) |>
  dplyr::rename_with(~ sample_map[.x], dplyr::ends_with(".mzXML"))
#- 7.3.3: Clean HILIC positive table and rename sample columns
ft_hilicpos_clean <- ft_hilicpos |>
  dplyr::mutate(Feature = paste(mz, time, sep = "_")) |>
  dplyr::select(Feature, dplyr::ends_with(".mzXML")) |>
  dplyr::rename_with(~ sample_map[.x], dplyr::ends_with(".mzXML"))

#+ 7.4: Validate sample name mapping
#- 7.4.1: Extract raw .mzXML column names from each table
c18neg_raw_cols   <- grep("\\.mzXML$", colnames(ft_c18neg),   value = TRUE)
hilicpos_raw_cols <- grep("\\.mzXML$", colnames(ft_hilicpos), value = TRUE)
#- 7.4.2: Check for unmapped columns
c18neg_unmapped   <- c18neg_raw_cols[!c18neg_raw_cols   %in% names(sample_map)]
hilicpos_unmapped <- hilicpos_raw_cols[!hilicpos_raw_cols %in% names(sample_map)]
#- 7.4.3: Stop on any unmapped columns, otherwise report success
if (length(c18neg_unmapped) > 0)
  stop("C18neg: ", length(c18neg_unmapped), " column(s) not found in study_design:\n",
       paste(c18neg_unmapped, collapse = "\n"))
if (length(hilicpos_unmapped) > 0)
  stop("HILIC pos: ", length(hilicpos_unmapped), " column(s) not found in study_design:\n",
       paste(hilicpos_unmapped, collapse = "\n"))
message("Mapping PASSED: ",
        length(c18neg_raw_cols),   " C18neg and ",
        length(hilicpos_raw_cols), " HILIC pos sample columns mapped 1:1.")

#+ 7.5: Write MWB-formatted feature tables
#- 7.5.1: Create output directory if it does not exist
dir.create(config$mwb_uploads_dir, recursive = TRUE, showWarnings = FALSE)
#- 7.5.2: Write C18 negative feature table
readr::write_tsv(ft_c18neg_clean,   file.path(config$mwb_uploads_dir, "c18neg_untargeted.tsv"))
#- 7.5.3: Write HILIC positive feature table
readr::write_tsv(ft_hilicpos_clean, file.path(config$mwb_uploads_dir, "hilicpos_untargeted.tsv"))

