#!/usr/bin/env Rscript
#' Setup and Install All Dependencies
#' 
#' This script installs all required packages for the DTC metabolomics pipeline.
#' It should be run once to set up the renv environment.
#' 
#' Usage: Rscript setup_dependencies.R

cat("\n=== DTC Metabolomics Pipeline - Dependency Setup ===\n\n")

# Ensure renv is active
if (!("renv" %in% loadedNamespaces())) {
  cat("📦 Loading renv...\n")
  source("renv/activate.R")
}

cat("✅ renv activated\n\n")

# Step 1: Install CRAN packages from DESCRIPTION
cat("📦 Step 1: Installing CRAN packages...\n")
cran_packages <- c(
  "broom", "conflicted", "cowplot", "dplyr", "forcats", 
  "ggplot2", "ggprism", "ggraph", "gtools", "here", 
  "igraph", "jsonlite", "magick", "memoise", "mice", 
  "permute", "pheatmap", "purrr", "RColorBrewer", "readr", 
  "readxl", "rlang", "scales", "stringr", "tibble", "tidyr", 
  "vegan", "yaml"
)

missing_cran <- cran_packages[!sapply(cran_packages, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("   Installing", length(missing_cran), "CRAN packages...\n")
  for (pkg in missing_cran) {
    cat("   -", pkg, "...\n")
    tryCatch({
      renv::install(pkg)
    }, error = function(e) {
      cat("     ⚠️ Failed:", pkg, "\n")
    })
  }
} else {
  cat("   ✅ All CRAN packages already installed\n")
}

# Step 2: Install BiocManager if needed
cat("\n📦 Step 2: Installing BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  renv::install("BiocManager")
  cat("   ✅ BiocManager installed\n")
} else {
  cat("   ✅ BiocManager already installed\n")
}

# Step 3: Install Bioconductor packages
cat("\n📦 Step 3: Installing Bioconductor packages...\n")
bioc_packages <- c("mixOmics", "KEGGREST")
missing_bioc <- bioc_packages[!sapply(bioc_packages, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) {
  cat("   Installing", length(missing_bioc), "Bioconductor packages...\n")
  for (pkg in missing_bioc) {
    cat("   -", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    }, error = function(e) {
      cat("     ⚠️ Failed:", pkg, "\n")
    })
  }
} else {
  cat("   ✅ All Bioconductor packages already installed\n")
}

# Step 4: Install GitHub packages
cat("\n📦 Step 4: Installing GitHub packages...\n")
github_packages <- list(
  list(repo = "xia-lab/MetaboAnalystR", name = "MetaboAnalystR")
)

for (pkg_info in github_packages) {
  if (!requireNamespace(pkg_info$name, quietly = TRUE)) {
    cat("   -", pkg_info$name, "from", pkg_info$repo, "...\n")
    tryCatch({
      renv::install(pkg_info$repo)
      cat("     ✅ Installed\n")
    }, error = function(e) {
      cat("     ⚠️ Failed to install", pkg_info$name, "\n")
      cat("     Error:", conditionMessage(e), "\n")
      cat("     You may need to install this package manually.\n")
      cat("     Try: remotes::install_github('", pkg_info$repo, "')\n", sep = "")
    })
  } else {
    cat("   ✅", pkg_info$name, "already installed\n")
  }
}

# Step 5: Create snapshot
cat("\n📸 Step 5: Creating renv snapshot...\n")
tryCatch({
  renv::snapshot(prompt = FALSE)
  cat("✅ Snapshot created successfully!\n")
}, error = function(e) {
  cat("⚠️ Snapshot creation had issues:\n")
  cat("   ", conditionMessage(e), "\n")
  cat("   You can create the snapshot manually later with: renv::snapshot()\n")
})

cat("\n=== Setup Complete ===\n")
cat("\nYou can now run the analysis with:\n")
cat("  source('all_run/run.R')\n\n")
