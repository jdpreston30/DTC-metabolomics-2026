# Metabolic reprogramming in differentiated thyroid cancer progression: a tumor metabolomics analysis

Analysis code for an untargeted LC–MS metabolomics study of 55 differentiated thyroid cancer (DTC)
specimens, comparing early-stage (AJCC I/II, n = 49) with advanced-stage (III/IV, n = 6) tumors.

## Citation

> Preston JD, Glosser LD, Jarrell ZR, Senanayake MS, Tran V, Weinberg J, Robertson JM, Weber CJ,
> Smith MR, Liang Y, Zhan J, Safley SA, Jackson AS, Patel SG, Sharma J, Saunders ND, Go Y-M,
> Jones DP, Szabo Yamashita T. Metabolic reprogramming in differentiated thyroid cancer progression:
> a tumor metabolomics analysis. *Surgery*. In press, 2026.

## Data availability

**No data files are included in this repository.** Clinical metadata and specimen identifiers are
withheld: the source biorepository retains a re-identification key, so specimen accession numbers are
not publicly shareable. Subjects are referred to throughout by study-assigned identifiers
(`F1`, `P1`, `FVPTC2`, …) that cannot be linked back to the biorepository without that key.

Raw spectral data (`.mzXML`) and untargeted feature tables are deposited at the NIH Common Fund's
National Metabolomics Data Repository, **Metabolomics Workbench** (Study ID and DOI pending release).
The deposition carries per-sample age, sex, tumor type and AJCC stage, keyed to the same
study-assigned identifiers used here.

To rerun the pipeline you must supply the input paths yourself — see `all_run/config_dynamic.yaml`.

## Requirements

- **R ≥ 4.5.1**
- **System tools**, verified by `check_system_dependencies()`: ImageMagick (`magick`),
  Ghostscript (`gs`), and a LaTeX engine (`tinytex` / `xelatex`) for supplementary table rendering
- A standard desktop machine; no special hardware

Package versions are pinned with **renv** (`renv.lock`), and dependencies are declared in
`DESCRIPTION` under `Imports:`, `Bioconductor:` and `Remotes:` with
`Config/renv/snapshot/type: explicit`.

## Running the analysis

```bash
git clone https://github.com/jdpreston30/DTC-metabolomics-2026.git
cd DTC-metabolomics-2026
```

```r
# renv activates automatically via .Rprofile
renv::restore()                 # first time only, ~10-20 min

source("R/Utilities/Setup/check_system_dependencies.R")
check_system_dependencies()

# edit all_run/config_dynamic.yaml so the input paths match your system, then:
source("all_run/run.R")
```

`run.R` executes the analysis scripts (`00a`–`08`) in order. The Metabolomics Workbench deposition
scripts are run separately, on demand — see below.

## Pipeline

### Analysis (`R/Scripts/`)

| Script | Purpose |
|---|---|
| `00a_environment_setup.R` | Conflict preferences, package loading from `DESCRIPTION`, system-dependency check |
| `00b_setup.R` | Sources every utility under `R/Utilities/` |
| `00c_clinical_metadata.R` | Clinical and pathology metadata, stage assignment |
| `00d_FTs.R` | Feature table import and preprocessing |
| `01_volcano_heatmap.R` | Differential features; volcano plot and heatmap |
| `02_pathway_enrichment.R` | Mummichog pathway enrichment (KEGG, MFN) |
| `03_annotated.R` | Annotated metabolite plots |
| `04_assign_plots.R` | Assigns panels to figures |
| `05_render_figures.R` | Renders publication figures |
| `06_data_not_shown.R` | Supporting analyses not shown in the manuscript |
| `07_tables.R` | Manuscript and supplementary tables |
| `08_session_info.R` | Writes `session_info.txt` |

### Metabolomics Workbench deposition (`R/MWB/`)

Not part of `run.R` — run these on demand when preparing a deposition. Shell steps (`00`–`03`)
organise, compress and stage the `.mzXML` archives; R steps build the deposition metadata. Paths are
configured separately in `R/MWB/config.yaml`.

Run the R steps in the order **04 → 06 → 05 → 07**: `05` summarises the deposited cohort by reading
the `study_design.tsv` that `06` writes, so running it earlier reports a stale cohort.

| Script | Purpose |
|---|---|
| `Scripts/00_setup.sh` | Checks/installs p7zip |
| `Scripts/01_organize.sh` | Sorts `.mzXML` into `hilicpos/` and `c18neg/` |
| `Scripts/02_compress.sh` | Builds `hilicpos.7z` / `c18neg.7z` and `zip_contents.csv` |
| `Scripts/03_copy.sh` | Stages archives for upload |
| `Scripts/04_clean.R` | Joins the acquisition sequence to the archive contents |
| `Scripts/06_study_design.R` | Builds the de-identified `SUBJECT_SAMPLE_FACTORS` table |
| `Scripts/05_custom_info.R` | Study summary — **runs after 06**, since it reads `study_design.tsv` |
| `Scripts/07_feature_tables.R` | Renames feature-table columns to study identifiers for upload |

## Repository layout

```
├── DESCRIPTION                 # dependencies (CRAN, Bioconductor, GitHub/R-universe)
├── renv.lock                   # pinned package versions
├── session_info.txt            # session record from the manuscript run
├── all_run/
│   ├── config_dynamic.yaml     # input paths and analysis options — edit before running
│   └── run.R                   # pipeline entry point
├── R/
│   ├── Scripts/                # analysis workflow, 00a-08
│   ├── MWB/                    # Metabolomics Workbench deposition
│   └── Utilities/              # Analysis, Helpers, Other, Preprocessing,
│                               #   Setup, Tabulation, Terminal, Visualization
├── Outputs/                    # Annotation, Correlation, Figures, Tables, mummichog
└── Supplementary/              # Components (incl. Tables), Build_Logs
```

## Contact

**Joshua D. Preston** — joshua.preston@emory.edu ·
[ORCID 0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017) ·
Medical Scientist Training Program, Emory University School of Medicine

**Corresponding author: Thomas Szabo Yamashita** — thomas.szabo.yamashita@emory.edu ·
[ORCID 0000-0002-2786-6678](https://orcid.org/0000-0002-2786-6678) ·
Division of General and GI Surgery, Department of Surgery, Emory University School of Medicine
