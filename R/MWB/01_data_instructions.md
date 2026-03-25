# MWB Submission — Instructions

> **Mac only.** Study-specific paths live in `00_config.yaml`, not here.
This is a recipe for metabolomics workbench submissions that can apply to any study.

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

Set `mzxml_dir` in `00_config.yaml` (the study mzXML folder — may be flat or pre-sorted), then:

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

**Requires** `hilicpos/` and `c18neg/` subfolders to exist in `mzxml_dir` (run Step 1 first). Compresses each into a `.7z` archive written as a sibling of the mzXML folder, and generates `zip_contents.csv`.

---

## Step 3: Copy archives to cloud destination

Set `raw_zips_dir` in `00_config.yaml` (the standalone `raw_zips/` folder on OneDrive/cloud), then:

```bash
bash R/MWB/Scripts/03_copy.sh
```

Creates `raw_zips/` inside `destination_dir` and copies the `.7z` archives, `zip_contents.csv`, and any `.sld` files there from the external drive. Original files are not modified.

---

## Step 4: Clean and validate

From a plain zsh terminal (**outside renv**):

```bash
R --no-init-file
source("R/MWB/Scripts/04_clean.R")
```

Reads the sequence file and `zip_contents.csv`, cleans both, and joins on `file_name_base`. Runs a validation check — must pass (all files matched 1:1) before proceeding. Fix any mismatches in the sequence file or re-run Step 2 if files are missing.

---

## Step 5: Create study design CSV
This is going to be context dependent for your study, but you now have zip_contents_mzxml which is a tibble with columns folder, subfolder, file_name, file_name_base, sample_id, and batch. Modify the 05_study_design.R script to join data an export a CSV of the study design which can be pasted into Metabolomics Workbench (i.e., including Subject_ID, Sample_ID,	Sample, source,	severe_PGD,	Batch,	RAW_FILE_NAME,	Grade,	sample_base,	replicate,	sample_type). Once this script is run, you direclty paste the CSV output into the Study Design box in Metabolomics Workbench.

```bash
R --no-init-file
source("R/MWB/Scripts/05_study_design.R")
```

---

## Step 6: Pull any other custom info needed for study metadata
If needed, run a custom info script to pull study metadata (i.e., age range, sex n's, study groups).. change this according to your study/tibble design from your main pipeline. This will also write the results into a markdown in data/ subfolder:

```bash
R --no-init-file
source("R/MWB/Scripts/06_custom_info.R")
```

