# DTC Metabolomics Analysis
This repository contains code and data used for the analysis presented in the submitted manuscript 'Title TBD' by Glosser et al. 2025.

To clone this repository to your local machine:
```bash
git clone https://github.com/jdpreston30/DTC-metabolomics-2025
cd DTC-metabolomics-2025
```

Make sure you have [Git](https://git-scm.com/) installed. After cloning, you can open all files in RStudio or your preferred IDE to begin the analysis.

---

## Requirements
- R version ≥ 4.3.1 is recommended for running all analyses in this repository
- **Required Packages:**
  - All packages can be installed and loaded under the subheader "#+ 0.1: Dependencies" in `R/Scripts/00_setup_and_import.r`

---

## Repository Structure

### Top-level
- `end_to_end.R` – Wrapper script to run the entire pipeline
- `all_objects.rds` – Serialized checkpoint of imported raw data (auto-generated on first run)
- `JDP_R_style_guide.md` – Internal coding/style notes
- `copilot-instructions.md` – Internal development notes

### R/

Core R code is organized into two main directories:

#### R/Scripts/
Sequential scripts that implement the full analysis pipeline:
- `00_setup_and_import.r` – package setup, utility sourcing, data import (or checkpoint restore)
- `01_cleanup.r` – feature table preprocessing (targeted and untargeted)
- `02_combined_clustering.r` – PCA, clustering, heatmaps, PERMANOVA
- `03_pathway_enrichment.r` – Mummichog enrichment analyses and visualization
- `04_cluster_var_path.r` – Pathology/cluster integration and plots
- `05_targeted.r` – Targeted metabolite statistical comparisons
- `06_figure_creation.r` – Assembly of final figures for manuscript

#### R/Utilities/
Modular function library (sourced automatically in `00_setup_and_import.r`):
- **Analysis/** – statistical and enrichment functions
- **Preprocessing/** – feature table and pathology preprocessing utilities
- **Visualization/** – plotting themes and figure export helpers

### Outputs/
- **Grob/** – ggplot2 objects exported as PNGs for later figure assembly
- **Mummichog Inputs/** – input CSVs for enrichment analysis
- **Mummichog Outputs/** – Excel outputs from Mummichog

### Figures/
- Final figures for the manuscript (`Figure 1.png`, `Figure 2.png`, `Figure 3.png`, `Supplemental Figure 1.png`)

---

## How to Run

There are two ways to execute the complete analysis:

### Option 1: Run the entire pipeline at once
```r
source("end_to_end.R")
```

### Option 2: Run scripts sequentially
Execute the scripts in order:
```r
source("R/Scripts/00_setup_and_import.r")  # Setup and data import
source("R/Scripts/01_cleanup.r")           # Data preprocessing
source("R/Scripts/02_combined_clustering.r") # PCA and clustering
source("R/Scripts/03_pathway_enrichment.r")  # Pathway enrichment
source("R/Scripts/04_cluster_var_path.r")    # Pathology integration
source("R/Scripts/05_targeted.r")            # Targeted analysis
source("R/Scripts/06_figure_creation.r")     # Figure assembly
```

**Note:** On first run, raw data will be imported and saved as a checkpoint (`all_objects.rds`). Subsequent runs will load from this checkpoint for faster execution.

---

## Accessing Raw Data

Raw data files are stored securely via OneDrive.
- Please request access by contacting me directly (joshua.preston@emory.edu)
- OR contact the corresponding author for permission
