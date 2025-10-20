# Configuration-Based Setup and Import

This document describes the new hybrid configuration approach for managing dependencies and data imports in the DTC Metabolomics analysis pipeline.

## Overview

The setup and import process has been refactored to use:
- **DESCRIPTION file** for package dependencies (standard R package approach)
- **YAML configuration** for data import specifications
- **Modular setup scripts** (00a, 00b pattern)

This approach combines the best of R package standards with flexible configuration management.

## Configuration Files

### `DESCRIPTION`
Standard R package DESCRIPTION file containing:
- **Package dependencies**: CRAN packages in `Imports:` section
- **Bioconductor packages**: Listed in separate `Bioconductor:` section
- **Project metadata**: Name, description, version, authors
- **Repository information**: URL, bug reports

### `config_imports.yaml`
YAML configuration for data import specifications:
- **CSV file configurations**: File names, variable names, descriptions
- **Excel file configurations**: File names, sheet names, variable names, descriptions
- **Checkpoint settings**: RDS file specifications

## Setup Scripts

### `00a_environment_setup.R`
Handles package installation and loading:
- Reads DESCRIPTION file to extract package lists
- Installs missing CRAN packages
- Installs missing Bioconductor packages
- Loads all packages with conflict-safe loading

### `00b_setup.R` 
Handles configuration and environment setup:
- Loads utility functions
- Loads dynamic configuration (paths, computer detection)
- Sets up R options and preferences
- Configures package conflict resolution
- Sets random seed for reproducibility

## Usage

### Main Setup Script
The `00_setup_and_import.r` script now follows this pattern:
```r
#+ 0.0: Environment setup (packages from DESCRIPTION file)
  source(here::here("R", "Scripts", "00a_environment_setup.R"))

#+ 0.1: Configuration setup (paths, conflicts, options)  
  source(here::here("R", "Scripts", "00b_setup.R"))

#+ 0.5: Import Data (from YAML configuration)
  import_config <- load_import_config("config_imports.yaml")
  # ... rest of import logic
```

## Benefits

### 1. Documentation
- Each package has a clear description of its purpose
- Data files include descriptions explaining their contents
- Configuration is self-documenting

### 2. Maintainability
- Easy to add/remove packages without editing R code
- Centralized dependency management
- Clear separation of configuration from logic

### 3. Flexibility
- Different configurations for different environments
- Easy to create project-specific setups
- Version-controlled configuration changes

### 4. Reproducibility
- Explicit package versions can be added to configuration
- Seed and project metadata tracked
- Clear audit trail of dependencies

## Configuration Structure

### Dependencies Section
```yaml
dependencies:
  description: "Core packages required for metabolomics analysis"
  cran_packages:
    - "dplyr"    # Data manipulation
    - "ggplot2"  # Data visualization
    # ... more packages
  bioconductor_packages:
    - name: "mixOmics"
      description: "Multivariate analysis for omics data"
  conflicts:
    - package: "dplyr"
      function: "select"
      description: "Use dplyr select over base select"
```

### Data Import Section
```yaml
data_import:
  csv_files:
    - file: "C18neg_UFT_medsum"
      name: "UFT_C18"
      description: "C18 negative mode untargeted feature table"
  xlsx_files:
    - file: "tumor_pathology"
      sheet: "pathology"
      name: "tumor_pathology_raw"
      description: "Raw tumor pathology data"
```

## Migration Notes

### What Changed
- Replaced hardcoded package vector with YAML configuration
- Replaced hardcoded tibble definitions with configuration-driven approach
- Added comprehensive documentation and descriptions
- Centralized conflict resolution

### What Stayed the Same
- Overall script flow and structure
- Data loading functions (`load_with_checkpoints`, etc.)
- Variable names and downstream compatibility
- Checkpoint and caching system

## Future Enhancements

1. **Package Versioning**: Add specific version requirements to configuration
2. **Environment Profiles**: Different configurations for dev/prod environments
3. **Validation**: Schema validation for configuration files
4. **Auto-documentation**: Generate documentation from configuration
5. **Dependency Graphs**: Visualize package dependencies

## Utility Functions

### `load_setup_config()`
- Loads and validates YAML configuration
- Returns structured configuration object
- Includes metadata about loading context

### `install_and_load_packages()`
- Installs missing packages (CRAN and Bioconductor)
- Loads packages with conflict resolution
- Provides verbose progress reporting

### `create_import_tibbles()`
- Converts YAML configuration to tibbles
- Maintains compatibility with existing import functions
- Preserves file descriptions for documentation