#* 04_clean.R  (Step 4)
#* Reads the sequence file and zip_contents.csv, cleans and joins them.
#* Run outside renv: R --no-init-file → source("R/MWB/Scripts/04_clean.R")

library(readr)
library(dplyr)
library(stringr)

#+ Config
cfg <- yaml::read_yaml("R/MWB/config.yaml")

#+ Import
sequence     <- read_tsv(cfg$sequence_file)
zip_contents <- read_tsv(file.path(cfg$raw_zips_dir, "zip_contents.csv"))

#+ Clean Sequence
sequence_clean <- sequence |>
  rename(file_name_base = `File Name`, sample_id = `Sample ID`)
#+ Clean up zip contents
#- Extract base file names
zip_contents_clean <- zip_contents |>
  mutate(file_name_base = str_extract(file_name, "^[^.]+"))
#- Split into mzXML and raw files; join with sequence
zip_contents_mzxml <- zip_contents_clean |>
  filter(str_ends(file_name, ".mzXML")) |>
  left_join(sequence_clean, by = "file_name_base")

#+ Validation: 1:1 match between sequence and zip_contents_mzxml
seq_bases  <- sort(sequence_clean$file_name_base)
zip_bases  <- sort(zip_contents_mzxml$file_name_base)

in_seq_not_zip  <- setdiff(seq_bases, zip_bases)
in_zip_not_seq  <- setdiff(zip_bases, seq_bases)
has_unmatched   <- nrow(filter(zip_contents_mzxml, is.na(sample_id))) > 0

if (length(in_seq_not_zip) == 0 && length(in_zip_not_seq) == 0 && !has_unmatched) {
  message("Validation PASSED: all ", length(seq_bases), " files matched 1:1 between sequence and zip contents.")
} else {
  if (length(in_seq_not_zip) > 0)
    warning("In sequence but NOT in zip: ", paste(in_seq_not_zip, collapse = ", "))
  if (length(in_zip_not_seq) > 0)
    warning("In zip but NOT in sequence: ", paste(in_zip_not_seq, collapse = ", "))
  if (has_unmatched)
    warning("zip_contents_mzxml has rows with no matching sequence entry (NA sample_id).")
  stop("Validation FAILED — review warnings above before proceeding.")
}
