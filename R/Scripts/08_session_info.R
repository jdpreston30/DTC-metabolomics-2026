#* 8: Session Information
#' Captures the exact R and package versions used for this analysis.
#' This script only runs if the entire pipeline completed successfully.
#+ 8.1: Generate and save session info
#- 8.1.1: Call session_info function to generate comprehensive documentation
session_info()
cat("\n✓ Pipeline completed successfully!\n")
