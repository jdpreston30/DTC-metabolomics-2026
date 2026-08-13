#* 05_custom_info.R — MWB Study Summary
#! Derived from study_design.tsv, NOT the pathology table: the pathology table is the full cohort (60),
#! but 5 subjects have no LC-MS runs (P5, P19, F5, FVPTC1, F3), so the deposit covers 55. Summarising the
#! cohort would report subject/sex counts that do not describe the uploaded data. RUN AFTER 06.
#+ Read the deposited study design
study_design <- readr::read_tsv("R/MWB/data/study_design.tsv", show_col_types = FALSE)
#+ One row per deposited subject (QC/plasma rows carry NA Subject_ID and are excluded)
deposited <- study_design |>
  dplyr::filter(!is.na(Subject_ID), Subject_ID != "NA", !is.na(Age)) |>
  dplyr::distinct(Subject_ID, .keep_all = TRUE)
stopifnot(nrow(deposited) > 0, !any(duplicated(deposited$Subject_ID)))
#+ Study summary from the deposited subjects
mwb_subject_summary <- c(
  "Number of Groups"  = as.character(dplyr::n_distinct(deposited$stage_bin)),
  "Total Subjects"    = as.character(nrow(deposited)),
  "Number of Males"   = as.character(sum(deposited$Sex == "Male",   na.rm = TRUE)),
  "Number of Females" = as.character(sum(deposited$Sex == "Female", na.rm = TRUE)),
  "Age or Age Range"  = paste0(min(deposited$Age, na.rm = TRUE), "-", max(deposited$Age, na.rm = TRUE))
)
lines <- paste0(formatC(paste0(names(mwb_subject_summary), ":"), width = -22), mwb_subject_summary)
cat(lines, sep = "\n")
#+ Write to data/
out_path <- "R/MWB/data/study_summary.md"
writeLines(c("# MWB Study Summary", "", lines), out_path)
message("Summary written to ", out_path)
