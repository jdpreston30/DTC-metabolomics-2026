#* 7: Tables
#+ 7.1: Table 1
#- 7.1.0: Create vector of included samples
used_samples <- UFT_filtered |>
  pull(ID)
#- 7.1.1: Filter to used samples; order rows; clean data
tumor_pathology_table <- tumor_pathology_full |>
  rename(ID = Patient_ID) |>
  filter(ID %in% used_samples) |>
  mutate(
    ETE = case_when(
      ETE == 0 ~ "Negative",
      ETE == 1 ~ "Minimal",
      ETE == 2 ~ "Extensive",
      TRUE ~ as.character(ETE)
    ),
    LVI = case_when(
      LVI == 0 ~ "N",
      LVI == 1 ~ "Y",
      LVI == 2 ~ "N", #Indeterminate treated as No
      TRUE ~ as.character(LVI)
    ),
    MFC = case_when(
      MFC == 0 ~ "No",
      MFC == 1 ~ "Yes",
      TRUE ~ as.character(MFC)
    ),
    ETE = factor(ETE, levels = c("Negative", "Minimal", "Extensive")),
    Variant = factor(Variant, levels = c("FTC", "FV-PTC", "PTC"))
  ) |>
  select(`Age (y)` = Age, Sex, stage_bin, `Pathologic Stage` = Stage, `T Category` = T_stage_comp, `N Category` = N, `M Category` = M, Variant, `Extrathyroidal Extension` = ETE, `Lymphovascular Invasion (Positive)` = LVI, `Multifocality (Positive)` = MFC)
#- 7.1.2: Build Table 1
T1 <- ternG(
  data = tumor_pathology_table,
  group_var = "stage_bin",
  show_p = FALSE,
  methods_doc = FALSE,
  output_docx = "Outputs/Tables/T1.docx",
  consider_normality = TRUE,
  table_font_size = 11,
  zero_to_dash = TRUE,
  table_caption = "Table I. Clinical and pathological characteristics of patients and analyzed tumors.",
  table_footnote = "All values are displayed as mean ± SD for ratio continuous variables or n (%) for categorical variables.",
  category_start = c(
    "Patient Characteristics" = "Age (y)",
    "Staging"                 = "Pathologic Stage",
    "Histologic Features"     = "Variant"
  ),
  variable_footnote = c(
    "Pathologic Stage|T Category" = "Thyroid cancer stage and T category were assigned using AJCC 8th-edition criteria regardless of what was listed on original operative pathology reports.",
    "Lymphovascular Invasion (Positive)" = "n = 1 patient had indeterminate lymphovascular invasion but was treated as \"negative\" analytically."
  ),
  abbreviation_footnote = "AJCC, American Joint Committee on Cancer; FTC, follicular thyroid carcinoma; FV-PTC, follicular variant of papillary thyroid carcinoma; M, metastasis; N, lymph node; PTC, papillary thyroid carcinoma; SD, standard deviation; T, tumor"
)
#+ 7.2: Supplementary Table 1
#- 7.2.1: Define group order
group_order <- c(
  "Nucleotide Flux",
  "Protein/Amino Acid Turnover", 
  "Membrane Integrity",
  "Lipid Remodeling",
  "Epigenetic Signaling",
  "Steroid Metabolism",
  "Redox Homeostasis",
  "Fatty Acid Oxidation",
  "Bioenergetic Flux",
  "Glycan Defense",
  "Immune/Signaling"
)
#- 7.2.2: Prepare base ST1 data (pre-sort by p-value before formatting)
ST1_base <- readxl::read_xlsx(config$data_files$QC, sheet = "QC") |>
  mutate(
    mode_ESI = case_when(
      Mode == "C18" ~ "C18-",
      Mode == "HILIC" ~ "HILIC+",
      TRUE ~ Mode
    ),
    # Remap group names to match group_order
    main_group = case_when(
      main_group == "Protein/AA Turnover" ~ "Protein/Amino Acid Turnover",
      TRUE ~ main_group
    ),
    # Create factor for group ordering
    Group = factor(main_group, levels = group_order)
  ) |>
  arrange(Group, p_value) |>
  select(
    Name = table_name,
    Group = main_group,
    Subgroup = Subcategory,
    p_value_raw = p_value,
    q_value_raw = p_value_fdr,
    log2FC_raw = log2FC,
    mz,
    `RT (s)` = `Retention Time`,
    `Column/ESI` = mode_ESI,
    Adduct,
    `KEGG ID` = KEGG,
    CID
  )
#- 7.2.3: Build hierarchical structure with Group: Name nesting
ST1_list <- list()
all_groups <- group_order[group_order %in% unique(ST1_base$Group)]

for (group_idx in seq_along(all_groups)) {
  group <- all_groups[group_idx]
  
  # Add GROUP header row
  ST1_list[[length(ST1_list) + 1]] <- tibble(
    Display_Name = group,
    row_type = "GROUP",
    Subgroup = NA_character_,
    `P Value` = NA_character_,
    `q Value` = NA_character_,
    `log2(FC)` = NA_character_,
    mz = NA_character_,
    `RT (s)` = NA_character_,
    `Column/ESI` = NA_character_,
    Adduct = NA_character_,
    `KEGG ID` = NA_character_,
    CID = NA_character_
  )

  # Get metabolites for this group (already sorted by p-value within group)
  group_data <- ST1_base |> filter(Group == group)
  
  # Check if we need to insert "Redox Homeostasis (Continued)" header
  # This should appear before DH-Monapterin triphosphate
  if (group == "Redox Homeostasis" && "DH-Monapterin triphosphate" %in% group_data$Name) {
    # Split into two parts: before and including DH-Monapterin onwards
    dh_index <- which(group_data$Name == "DH-Monapterin triphosphate")
    
    # Add metabolites before DH-Monapterin
    if (dh_index > 1) {
      metabolites_part1 <- group_data |>
        slice(1:(dh_index - 1)) |>
        mutate(
          Display_Name = paste0("    ", Name),
          row_type = "Metabolite",
          # Format p and q values
          `P Value` = case_when(
            p_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", p_value_raw)),
            TRUE ~ sprintf("%.3f", p_value_raw)
          ),
          `q Value` = case_when(
            q_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", q_value_raw)),
            TRUE ~ sprintf("%.3f", q_value_raw)
          ),
          `log2(FC)` = sprintf("%.2f", log2FC_raw),
          mz = sprintf("%.4f", mz),
          `RT (s)` = as.character(`RT (s)`),
          `Column/ESI` = as.character(`Column/ESI`),
          Adduct = as.character(Adduct),
          `KEGG ID` = as.character(`KEGG ID`),
          CID = as.character(CID)
        ) |>
        select(Display_Name, row_type, Subgroup, `P Value`, `q Value`, `log2(FC)`, mz, `RT (s)`, `Column/ESI`, Adduct, `KEGG ID`, CID)
      
      ST1_list[[length(ST1_list) + 1]] <- metabolites_part1
    }
    
    # Add "Redox Homeostasis (Continued)" header
    ST1_list[[length(ST1_list) + 1]] <- tibble(
      Display_Name = "Redox Homeostasis (Continued)",
      row_type = "GROUP",
      Subgroup = NA_character_,
      `P Value` = NA_character_,
      `q Value` = NA_character_,
      `log2(FC)` = NA_character_,
      mz = NA_character_,
      `RT (s)` = NA_character_,
      `Column/ESI` = NA_character_,
      Adduct = NA_character_,
      `KEGG ID` = NA_character_,
      CID = NA_character_
    )
    
    # Add metabolites from DH-Monapterin onwards
    metabolites_part2 <- group_data |>
      slice(dh_index:n()) |>
      mutate(
        Display_Name = paste0("    ", Name),
        row_type = "Metabolite",
        # Format p and q values
        `P Value` = case_when(
          p_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", p_value_raw)),
          TRUE ~ sprintf("%.3f", p_value_raw)
        ),
        `q Value` = case_when(
          q_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", q_value_raw)),
          TRUE ~ sprintf("%.3f", q_value_raw)
        ),
        `log2(FC)` = sprintf("%.2f", log2FC_raw),
        mz = sprintf("%.4f", mz),
        `RT (s)` = as.character(`RT (s)`),
        `Column/ESI` = as.character(`Column/ESI`),
        Adduct = as.character(Adduct),
        `KEGG ID` = as.character(`KEGG ID`),
        CID = as.character(CID)
      ) |>
      select(Display_Name, row_type, Subgroup, `P Value`, `q Value`, `log2(FC)`, mz, `RT (s)`, `Column/ESI`, Adduct, `KEGG ID`, CID)
    
    ST1_list[[length(ST1_list) + 1]] <- metabolites_part2
    
  } else {
    # Normal processing for all other groups
    # Add metabolites with 4-space indent for visual separation
    metabolites <- group_data |>
      mutate(
        Display_Name = paste0("    ", Name),
        row_type = "Metabolite",
        # Format p and q values
        `P Value` = case_when(
          p_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", p_value_raw)),
          TRUE ~ sprintf("%.3f", p_value_raw)
        ),
        `q Value` = case_when(
          q_value_raw < 0.001 ~ gsub("e([+-])0+(\\d)", "E\\1\\2", sprintf("%.1e", q_value_raw)),
          TRUE ~ sprintf("%.3f", q_value_raw)
        ),
        `log2(FC)` = sprintf("%.2f", log2FC_raw),
        mz = sprintf("%.4f", mz),
        `RT (s)` = as.character(`RT (s)`),
        `Column/ESI` = as.character(`Column/ESI`),
        Adduct = as.character(Adduct),
        `KEGG ID` = as.character(`KEGG ID`),
        CID = as.character(CID)
      ) |>
      select(Display_Name, row_type, Subgroup, `P Value`, `q Value`, `log2(FC)`, mz, `RT (s)`, `Column/ESI`, Adduct, `KEGG ID`, CID)
    
    ST1_list[[length(ST1_list) + 1]] <- metabolites
  }

  
  # Add blank row after each GROUP (except the last)
  if (group_idx < length(all_groups)) {
    ST1_list[[length(ST1_list) + 1]] <- tibble(
      Display_Name = "",
      row_type = "Spacer",
      Subgroup = "",
      `P Value` = "",
      `q Value` = "",
      `log2(FC)` = "",
      mz = "",
      `RT (s)` = "",
      `Column/ESI` = "",
      Adduct = "",
      `KEGG ID` = "",
      CID = ""
    )
  }
}
#- 7.2.4: Combine and format final tibble
# Only replace NA with "–" for Metabolite rows (actual missing data)
# GROUP and Spacer rows should have blank cells
ST1_tibble <- bind_rows(ST1_list) |>
  mutate(across(everything(), as.character)) |>
  mutate(across(
    -c(Display_Name, row_type, Subgroup),
    ~ if_else(row_type == "Metabolite" & is.na(.x), "–", .x)
  )) |>
  mutate(across(everything(), ~ replace_na(.x, ""))) |>  # Remaining NAs become blank
  select(-row_type)  # Remove helper column
#- 7.2.5: Source table builder and generate LaTeX
source("R/Utilities/Tabulation/build_ST1_latex.R")
build_ST1_latex(
  data = ST1_tibble,
  group_names = group_order,
  output_path = "Supplementary/Components/Tables/ST1.tex"
)
#- 7.2.6: Render supplementary tables PDF
{
  components_dir <- here::here("Supplementary", "Components")
  intermediates_dir <- here::here("Supplementary", "Build_Logs")
  output_dir <- here::here("Supplementary")
  
  rmarkdown::render(
    input = file.path(components_dir, "supplementary_tables.Rmd"),
    output_dir = output_dir,
    output_file = "Supplementary Table I.pdf",
    intermediates_dir = intermediates_dir,
    clean = TRUE
  )
  
  # Open the PDF
  output_pdf <- file.path(output_dir, "Supplementary Table I.pdf")
  if (Sys.info()["sysname"] == "Darwin") {
    system(paste("open", shQuote(output_pdf)))
  }
} 
