#* 06_custom_info.R — MWB Study Summary
#+ Read tumor_pathology_table from MWB data
tumor_pathology_table <- read.csv("R/MWB/data/tumor_pathology_table.csv", check.names = FALSE)
tumor_pathology_table$stage_bin <- as.factor(tumor_pathology_table$stage_bin)
#+ Study summary from tumor_pathology_table
mwb_subject_summary <- c(
  "Number of Groups"  = as.character(nlevels(tumor_pathology_table$stage_bin)),
  "Total Subjects"    = as.character(nrow(tumor_pathology_table)),
  "Number of Males"   = as.character(sum(tumor_pathology_table$Sex == "Male",   na.rm = TRUE)),
  "Number of Females" = as.character(sum(tumor_pathology_table$Sex == "Female", na.rm = TRUE)),
  "Age or Age Range"  = paste0(min(tumor_pathology_table$Age, na.rm = TRUE), "-", max(tumor_pathology_table$Age, na.rm = TRUE))
)
lines <- paste0(formatC(paste0(names(mwb_subject_summary), ":"), width = -22), mwb_subject_summary)
cat(lines, sep = "\n")
#+ Write to data/
out_path <- "R/MWB/data/study_summary.md"
writeLines(c("# MWB Study Summary", "", lines), out_path)
message("Summary written to ", out_path)
