#!/usr/bin/env Rscript
#' Check renv Environment Status
#' 
#' This script provides a quick overview of the renv environment status
#' and helps diagnose any package management issues.

cat("\n=== DTC Metabolomics Pipeline - Environment Status ===\n\n")

# Load renv
if (!("renv" %in% loadedNamespaces())) {
  source("renv/activate.R")
}

# Check R version
cat("📊 R Version:\n")
cat("   Current:", as.character(getRversion()), "\n")
if (file.exists("renv.lock")) {
  lock <- jsonlite::fromJSON("renv.lock", simplifyVector = FALSE)
  if (!is.null(lock$R$Version)) {
    cat("   Lockfile:", lock$R$Version, "\n")
  }
}

# Check renv status
cat("\n📦 Package Library Status:\n")
status_output <- capture.output(renv::status(), type = "message")
cat(paste(status_output, collapse = "\n"))

# Check critical packages
cat("\n\n🔍 Critical Package Availability:\n")
critical_packages <- c(
  "dplyr", "ggplot2", "purrr", "readr", "tidyr",
  "mixOmics", "vegan", "here", "yaml",
  "MetaboAnalystR"
)

all_available <- TRUE
for (pkg in critical_packages) {
  is_available <- requireNamespace(pkg, quietly = TRUE)
  status_icon <- if (is_available) "✅" else "❌"
  cat("   ", status_icon, pkg, "\n")
  if (!is_available) all_available <- FALSE
}

# System dependencies
cat("\n🔧 System Dependencies:\n")
source("R/Utilities/Helpers/check_system_dependencies.R")

# Summary
cat("\n=== Summary ===\n")
if (all_available) {
  cat("✅ Environment is ready to run the analysis!\n")
  cat("   Run: source('all_run/run.R')\n\n")
} else {
  cat("⚠️  Some packages are missing.\n")
  cat("   Run: Rscript setup_dependencies.R\n")
  cat("   Or: renv::restore()\n\n")
}
