#* 8: Tables
#+ 8.1: Table 1
#- 8.1.0: Create vector of included samples
used_samples <- UFT_filtered |>
  pull(ID)
#- 8.1.1: Filter to used samples; order rows; clean data
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
    )
  ) |>
  select(Age, Sex, Variant, stage_bin, Stage, T = T_stage_comp, N, M, `Extrathyroidal extension` = ETE, `Lymphovascular invasion` = LVI, Multifocality = MFC)
#- 8.1.2: Build Table 1
T1 <- ternG(
  data = tumor_pathology_table,
  group_var = "stage_bin",
  descriptive = TRUE,
  output_docx = "Outputs/Tables/T1.docx",
  consider_normality = TRUE,
  show_test = FALSE,
  round_intg = FALSE,
  insert_subheads = TRUE
)
#+ 8.2: Supplementary Table 1
#- 8.2.1: Build Supplementary Table 1
pdf generation of st1
