# MWB Submission — Instructions

> **Mac only.** Study-specific paths live in `config.yaml`, not here.  
> This is a reusable recipe for Metabolomics Workbench submissions.

> **Prerequisites:** mzXML files must already exist (converted on Windows via msconvert/ProteoWizard). This workflow does not perform conversion.

---

## Step 0: Setup

```bash
bash R/MWB/Scripts/00_setup.sh
```

Checks/installs p7zip (required for compression).

---

## Step 1: Organize raw and mzXML files

> **Note:** mzXML files are large. It is strongly advised to keep them on an external drive.

Set `mzxml_dir` in `config.yaml` (the study mzXML folder — may be flat or pre-sorted), then:

```bash
bash R/MWB/Scripts/01_organize.sh
```

Checks for `.raw` files in the study folder and `.mzXML` files in `mzxml_dir`. For each:
- If already sorted into `hilicpos/` and `c18neg/` subfolders: reports and moves on.
- If flat (unsorted): prompts to sort by odd/even sequential acquisition order and moves files into subfolders.

---

## Step 2: Compress mzXML files

```bash
bash R/MWB/Scripts/02_compress.sh
```

**Requires** `hilicpos/` and `c18neg/` subfolders to exist in `mzxml_dir` (run Step 1 first). Compresses each into a `.7z` archive written as a sibling of the mzXML folder and generates `zip_contents.csv`.

---

## Step 3: Copy archives to cloud destination

Set `raw_zips_dir` in `config.yaml` (the standalone `raw_zips/` folder on OneDrive/cloud), then:

```bash
bash R/MWB/Scripts/03_copy.sh
```

Copies the `.7z` archives, `zip_contents.csv`, and any `.sld` files to `raw_zips_dir`. Original files are not modified.

---

## Step 4: Pull custom study metadata

Modify `05_custom_info.R` to match your study's metadata tibbles, then run it from a plain zsh terminal (**outside renv**):

```bash
R --no-init-file
source("R/MWB/Scripts/06_study_design.R")   # must run first — writes study_design.tsv
source("R/MWB/Scripts/05_custom_info.R")
```

> `05_custom_info.R` summarises the **deposited** cohort by reading `study_design.tsv`, so `06` has to
> run before it. The web form is filled in this order, but the scripts run 04 → 06 → 05 → 07.

Writes a summary markdown to `R/MWB/data/` with study-level stats (age range, sex n's, group breakdowns, etc.) for use when filling in the MWB web form.

### MWB web form — Basic Data Entry

Go to **Metabolomics Workbench → Upload/Manage Studies → New Submission** and enter:

- **Raw data files**: `c18neg.7z, hilicpos.7z`
- **Protocol/methods filename**: `methods.pdf` (or your study's filename)
- **MS instrument manufacturer**: `Thermo Scientific`
- **MS instrument model**: `Orbitrap Fusion Tribrid Mass Spectrometer` (or `Orbitrap ID-X Tribrid Mass Spectrometer` if using IDX)
- **Binary data format**: `.RAW`
- **Open source text format**: `mzXML`

### MWB web form — Project Information

You can largely reuse templates from previous projects. Key entries:
- **Project type**: `MS1 untargeted analysis`
- For n's and group breakdowns, use the summary stats written by `05_custom_info.R`.

---

## Step 5: Clean and validate

From a plain zsh terminal (**outside renv**):

```bash
R --no-init-file
source("R/MWB/Scripts/04_clean.R")
```

Reads the sequence file and `zip_contents.csv`, cleans both, and joins on `file_name_base`. Runs a validation check — must pass (all files matched 1:1) before proceeding. Fix any mismatches in the sequence file or re-run Step 2 if files are missing.

---

## Step 6: Build study design file

Modify `06_study_design.R` to join your study metadata and export the study design file, then run:

```bash
R --no-init-file
source("R/MWB/Scripts/06_study_design.R")
```

Produces `study_design.tsv` in `R/MWB/data/`. Paste this directly into the **Study Design** box in MWB. After pasting, click **View/Check Study Design** and assign columns:
- `Subject_ID`, `Sample_ID`, `sample_source`, `RAW_FILE_NAME` → their matching MWB categories
- Primary factor (e.g. `stage_bin`) → `Factor`
- All other columns (e.g. `Stage`, `Variant`, `replicate`, `sample_type`) → `Other`

---

## Step 7: Feature tables

Set `c18neg_feature_table`, `hilicpos_feature_table`, and `mwb_uploads_dir` in `config.yaml`, then:

```bash
R --no-init-file
source("R/MWB/Scripts/07_feature_tables.R")
```

Takes in xMSanalyzer-format untargeted feature tables (`RAW_mzcalibrated_untargeted_featuretable.txt`) for both methods, reshapes them into MWB-ready tab-delimited files, and writes them to `mwb_uploads_dir`. When uploading to MWB:

- **Units of measure**: `peak area`
- **Feature names contain m/z values**: yes
- **Feature names contain retention time values**: yes, in **seconds**

---

## Step 8: Sample Information

In the MWB web form, enter free text for each of the following sections:

**Collection Summary** — describe sample source and storage, e.g.:
> "Thyroid tumor samples were obtained from an institutional biorepository that collected intraoperative specimens from patients undergoing total thyroidectomy, subtotal thyroidectomy, or lobectomies. Thyroid tissues were stored at −80 °C following collection until analysis."

**Treatment Summary** — e.g.:
> "Samples were collected intraoperatively from patients undergoing total thyroidectomy, subtotal thyroidectomy, or lobectomies, and no treatments were applied."

**Sample Preparation** — e.g.:
> "Sample weights were recorded, and a solution containing two parts acetonitrile (with ~2% internal standard [IS]) and one part water was added to achieve a 15:1 volume:weight ratio. Tumor samples were then homogenized. Pooled reference plasma (Equitech-Bio SHP45 and NIST SRM1950) was also prepared and analyzed in parallel (50 μL plasma:100 μL ACN + ~2% IS). All samples were then incubated on ice for 30 minutes, followed by centrifugation at 20,817 × g for 15 minutes. Supernatants were transferred into autosampler vials and kept at −20 °C until analysis."

---

## Step 9: Chromatography

Enter chromatography details for both methods. You can copy from a previous study or use the entries below (HFQE and Fusion methods are the same).

### HILIC (positive ion mode)

- **Chromatography type**: HILIC
- **Instrument**: Thermo Dionex Ultimate 3000
- **Column**: Waters ACQUITY UPLC BEH Amide (100 × 2.1 mm, 1.7 μm)
- **Ion mode**: Positive
- **Injection volume**: 10 μL
- **Solvent A**: 100% LC-MS grade water
- **Solvent B**: 100% LC-MS grade acetonitrile
- **Solvent C**: 2% formic acid (v/v) in LC-MS grade water
- **Flow gradient**: Initial conditions: 22.5% A, 75% B, 2.5% C, hold for 1.5 min; linear gradient to 77.5% A, 20% B, 2.5% C at 4 min, hold for 1 min; total run time 5 min. HILIC column equilibrated with 77.5% A, 20% B, 2.5% C during reverse phase flushing phase.
- **Flow rate**: 0.35 mL/min (0–1.5 min); increased to 0.4 mL/min at 4 min, held for 1 min
- **Column temperature**: 60 °C
- **Sample injection**: 10 μL
- **Analytical time**: 5 min
- **Sample loop size**: 15 μL
- **Sample syringe size**: 100 μL

### C18 (negative ion mode)

- **Chromatography type**: Reversed phase
- **Instrument**: Thermo Dionex Ultimate 3000
- **Column**: Thermo Hypersil GOLD aQ (100 × 2.1 mm, 1.9 μm)
- **Ion mode**: Negative
- **Injection volume**: 10 μL
- **Solvent A**: 100% LC-MS grade water
- **Solvent B**: 100% LC-MS grade acetonitrile
- **Solvent C**: 10 mM ammonium acetate in LC-MS grade water
- **Flow gradient**: Initial conditions: 60% A, 35% B, 5% C, hold for 0.5 min; linear gradient to 0% A, 95% B, 5% C at 1.5 min, hold for 3.5 min; total run time 5 min. C18 column equilibrated with 0% A, 95% B, 5% C until 2.5 min, followed by re-equilibration at 60% A, 35% B, 5% C for 2.5 min during HILIC analytical separation.
- **Flow rate**: 0.4 mL/min (0–1.5 min); increased to 0.5 mL/min at 2 min, held for 3 min
- **Column temperature**: 60 °C
- **Sample injection**: 10 μL
- **Analytical time**: 5 min
- **Sample loop size**: 15 μL
- **Sample syringe size**: 100 μL
