# Metabolomic Profiling Identifies Functional Tumor Subtypes in Differentiated Thyroid Cancer

## 📖 Citation

This code is associated with the analysis presented in the following manuscript:
> Preston et al. (2026). Metabolomic Profiling Identifies Functional Tumor Subtypes in Differentiated Thyroid Cancer. *[Journal Name - Pending Publication]*.

## 🚀 Quick Start for Reproduction

**⚠️ Data Availability Notice**: 
- **No data files** (raw data, processed feature tables, or clinical metadata) are included in this repository
- **All instructions below assume you have obtained data files or are using your own data**
- **To reproduce this analysis**: Contact the author (Joshua D. Preston, joshua.preston@emory.edu) to obtain the data files—this is the easiest and recommended approach
- **Public data access**: [Pending - will be deposited to public repository]
- **To run analyses with your own data or provided data files**: Update file paths in `all_run/config_dynamic.yaml` to match your system

### Option 1: Using Docker (Recommended for Exact Reproducibility)

**Prerequisites**: 
- Install [Docker Desktop](https://www.docker.com/products/docker-desktop)
- (Optional) Create free [Docker Hub](https://hub.docker.com) account

#### Method A: Pull Pre-built Image (Fastest - Recommended)

```bash
# 1. Clone the repository
git clone https://github.com/jdpreston30/DTC-metabolomics-2026.git
cd DTC-metabolomics-2026

# 2. Pull the pre-built Docker image (~5-10 minutes)
docker pull jdpreston30/dtc-metabolomics:latest  # [PENDING]

# 3. Run the complete analysis pipeline
docker run -v $(pwd):/analysis jdpreston30/dtc-metabolomics:latest

# Windows users: Replace $(pwd) with:
# - CMD: docker run -v %cd%:/analysis jdpreston30/dtc-metabolomics:latest
# - PowerShell: docker run -v ${PWD}:/analysis jdpreston30/dtc-metabolomics:latest
```

#### Method B: Build Image Locally

```bash
# 1. Clone the repository
git clone https://github.com/jdpreston30/DTC-metabolomics-2026.git
cd DTC-metabolomics-2026

# 2. Build the Docker image from Dockerfile (~30-45 minutes)
docker build -t dtc-metabolomics .  # [PENDING - Dockerfile to be created]

# 3. Run the complete analysis pipeline
docker run -v $(pwd):/analysis dtc-metabolomics

# Windows users: Replace $(pwd) with:
# - CMD: docker run -v %cd%:/analysis dtc-metabolomics
# - PowerShell: docker run -v ${PWD}:/analysis dtc-metabolomics
```

All outputs (figures, tables, pathway results) will be saved to your local workspace.

#### Testing the Container

To verify the Docker image was built correctly before running the full analysis:

```bash
# Quick verification (< 1 minute)
docker run --rm dtc-metabolomics Rscript -e "packageVersion('mixOmics'); packageVersion('vegan')"
```

This should display package versions. If it succeeds, the container is ready for the full analysis.

#### Troubleshooting

**Build fails or is very slow**: 
- This is normal for the first build (30-45 minutes) due to compiling complex packages
- If build is interrupted, run `docker build` again - it will resume from where it stopped

**"No space left on device" error**:
- The final image is ~3-5 GB
- Check available disk space: `docker system df`
- Clean up old images: `docker image prune`

**Container runs but produces no output**:
- Check the Outputs directory exists and is writable
- Verify the mount path: `ls Outputs/`

### Option 2: Manual Installation (Without Docker)

**Prerequisites**: 
- R >= 4.5.1
- Git (to clone repository)

**Note**: This project uses `renv` for package management to ensure reproducibility. The `renv.lock` file contains exact versions of all packages used in the manuscript.

```r
# 1. Clone the repository
# (from terminal)
git clone https://github.com/jdpreston30/DTC-metabolomics-2026.git
cd DTC-metabolomics-2026

# 2. Start R in the project directory
# (renv automatically activates via .Rprofile)

# 3. Restore all packages at exact versions (first time only, ~10-20 minutes)
renv::restore()

# 4. Check system dependencies
source("R/Utilities/Setup/check_system_dependencies.R")
check_system_dependencies()

# 5. Run the complete analysis pipeline
source("all_run/run.R")
```

**What happens during `renv::restore()`**:
- Installs all R packages at exact versions from `renv.lock`
- Installs CRAN packages (e.g., ggplot2, dplyr, broom, conflicted)
- Installs Bioconductor packages (e.g., mixOmics, KEGGREST)
- Installs GitHub packages (MetaboAnalystR)
- Creates isolated project library (doesn't affect your system R packages)
- Only needed once per computer; subsequent runs use installed packages
- Packages are automatically loaded from `DESCRIPTION` file during pipeline execution

## 📁 Project Structure

```
├── DESCRIPTION                 # R package dependencies (CRAN, Bioconductor, GitHub)
├── Dockerfile                  # [PENDING] Docker container for reproducible environment
├── renv.lock                   # Exact package versions for reproducibility
├── all_run/                    # Pipeline execution
│   ├── config_dynamic.yaml     # Analysis configuration (update paths for your system)
│   └── run.R                   # Main pipeline execution script
├── R/                          # Analysis code
│   ├── Scripts/                # Analysis workflow scripts (00a-06)
│   ├── Utilities/              # Custom analysis functions
│   │   ├── Analysis/           # Statistical and pathway analysis
│   │   ├── Helpers/            # Helper functions
│   │   ├── Preprocessing/      # Data preprocessing functions
│   │   ├── Visualization/      # Plotting functions
│   │   └── Other/              # Other utilities
│   └── legacy/                 # Legacy code (not part of main pipeline)
├── Figures/                    # Static/reference figures
├── Outputs/                    # Generated results
│   ├── Figures/                # Publication figures (PDF, EPS, PNG)
│   └── mummichog/              # Pathway enrichment results (KEGG, MFN)
└── R/Utilities/Setup/          # Setup and configuration scripts
    ├── setup_dependencies.R    # One-time package installation
    └── check_environment.R     # Environment status checker
```

## 🔬 Analysis Workflow

The complete pipeline executes in sequence:

1. **00a-00d**: Environment setup, clinical metadata, feature tables
2. **01**: Clustering analysis
3. **02**: Pathway enrichment (mummichog analysis)
4. **03**: Annotated bar plots and visualizations
5. **04**: Assignment of plots to figure panels
6. **05**: Render final figures
7. **06**: Generate abstract/presentation data

## 💻 System Requirements

### Computational Requirements
- **R**: Version 4.5.1 or higher
- **Platform**: Developed on macOS but should work on Windows or Linux
- **Note**: Standard modern computer sufficient; no special hardware required

### System Dependencies
- **Ghostscript**: PDF to EPS conversion
- **ImageMagick**: Image processing

*Note: All system dependencies will be automatically installed in the Docker container. For manual installation, run `check_system_dependencies()` for platform-specific instructions.*

## 📦 Package Dependencies

All R package dependencies are specified in `DESCRIPTION`. Key packages include:

### CRAN Packages
- **Data manipulation**: tidyverse (dplyr, tidyr, purrr, readr, tibble, stringr, forcats)
- **Visualization**: ggplot2, ggraph, ggprism, pheatmap, cowplot, RColorBrewer, scales, magick
- **Statistical modeling**: vegan, mice, permute
- **Utilities**: here, yaml, conflicted, broom, jsonlite, readxl, memoise, gtools

### Bioconductor Packages
- **Metabolomics analysis**: mixOmics
- **Pathway databases**: KEGGREST

### GitHub Packages
- `xia-lab/MetaboAnalystR`: Metabolomics analysis and pathway enrichment

*See `DESCRIPTION` file for complete list of all dependencies.*

## � Package Management with renv

This project uses `renv` for reproducible package management. All package dependencies and their exact versions are tracked in `renv.lock`.

### How It Works

The project follows a structured approach to package management:

1. **DESCRIPTION file** - Contains all package dependencies:
   - `Imports:` - CRAN packages
   - `Bioconductor:` - Bioconductor packages  
   - `Remotes:` - GitHub packages

2. **Helper utilities** in `R/Utilities/Setup/`:
   - `restore_renv.R` - Ensures renv environment is ready
   - `load_packages_from_description.R` - Loads packages from DESCRIPTION
   - `check_system_dependencies.R` - Verifies system requirements

3. **00a_environment_setup.R** - Sets up the R environment:
   - Sets conflict preferences (using `conflicted` package)
   - Loads all packages from DESCRIPTION
   - Explicitly loads GitHub packages for renv detection
   - Checks system dependencies

### Explicit Snapshot Mode

This project uses `Config/renv/snapshot/type: explicit` in DESCRIPTION. This means:

- ✅ Only packages listed in DESCRIPTION are tracked
- ✅ Cleaner, more maintainable dependency list
- ✅ Easier to audit what packages are actually needed
- ✅ GitHub packages must be loaded explicitly with `library()` for renv to detect them

### Manual Package Management

#### Add a new package

1. Add to DESCRIPTION:
   ```
   Imports:
       ...,
       newpackage
   ```

2. Install and snapshot:
   ```r
   renv::install("newpackage")
   renv::snapshot()
   ```

#### Update packages

```r
# Update specific package
renv::update("dplyr")

# Update all packages
renv::update()

# Create new snapshot after updates
renv::snapshot()
```

#### Restore packages

If packages are missing or out of sync:

```r
# Check status
renv::status()

# Restore from lockfile
renv::restore()
```

#### Remove unused packages

```r
# Remove packages not in DESCRIPTION
renv::clean()
```

### Troubleshooting

#### "Package not found" errors

```r
# Check what's missing
renv::status()

# Restore from lockfile
renv::restore()
```

#### GitHub package installation fails

GitHub packages may require authentication. Set up a GitHub Personal Access Token:

```r
# Install usethis if needed
install.packages("usethis")

# Set up GitHub credentials
usethis::create_github_token()
gitcreds::gitcreds_set()
```

Then retry installation:
```r
renv::install("xia-lab/MetaboAnalystR")
```

#### Slow renv activation

If you see "renv took longer than expected to activate the sandbox", you can disable it:

Add to `.Renviron`:
```
RENV_CONFIG_SANDBOX_ENABLED=FALSE
```

#### Starting fresh

To completely rebuild the environment:

```r
# Remove renv library
unlink("renv/library", recursive = TRUE)

# Restore everything
renv::restore()
```

## �🔄 Reproducibility Features

This project implements best practices for computational reproducibility:

- ✅ **Version Control**: Complete analysis code on GitHub
- ✅ **Package Management**: `renv` with `renv.lock` pinning all packages to exact versions
- ✅ **Dependency Declaration**: All dependencies specified in `DESCRIPTION` with automatic loading
- ✅ **Containerization**: [PENDING] Docker image available at `jdpreston30/dtc-metabolomics:latest`
- ✅ **Conflict Resolution**: `conflicted` package ensures predictable function behavior
- ✅ **Docker Hub Distribution**: [PENDING] Pre-built image at [jdpreston30/dtc-metabolomics](https://hub.docker.com/r/jdpreston30/dtc-metabolomics)
- ✅ **Configuration-Driven**: All parameters in `config_dynamic.yaml`
- ✅ **System Dependency Checking**: Automated validation via `check_system_dependencies()`
- ✅ **Documentation**: Comprehensive function documentation and workflow comments

## 📧 Contact

**Author & Repository Maintainer**: Joshua D. Preston
- **Email**: joshua.preston@emory.edu  
- **ORCID**: [0000-0001-9834-3017](https://orcid.org/0000-0001-9834-3017)  
- **Institution**: Department of Surgery, Emory University School of Medicine

**Corresponding Author**: [Pending]

---

**Repository**: https://github.com/jdpreston30/DTC-metabolomics-2026  
**Docker Hub**: [PENDING] https://hub.docker.com/r/jdpreston30/dtc-metabolomics  
**Zenodo Archive**: [PENDING] https://doi.org/10.5281/zenodo.XXXXXXX
