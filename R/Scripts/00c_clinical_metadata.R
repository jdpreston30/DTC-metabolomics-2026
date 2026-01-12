#* 0c: Clinical metadata processing
#+ 0c.1: Cleanup tumor pathology data; compute T stage; compute overall stage
#- 0c.1.1: Create full version
tumor_pathology_full <- tumor_pathology_raw |>
  select(-T) |>
  assign_T_stage(ld_col = "LD", ete_col = "ETE", units = "cm", out_col = "T_stage_comp") |>
  mutate(
    Variant = factor(case_when(
      str_detect(Patient_ID, "^FVPTC\\d+$") ~ "FV-PTC",
      str_detect(Patient_ID, "^F\\d+$") ~ "FTC",
      str_detect(Patient_ID, "^P\\d+$") ~ "PTC",
      TRUE ~ "Unknown"
    ), levels = c("PTC", "FV-PTC", "FTC")),
    T_stage_comp = factor(T_stage_comp, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
    N = factor(case_when(
      N == 0 ~ "N0", 
      N == 1 ~ "N1", 
      TRUE ~ NA_character_
    )),
    M = factor(case_when(
      M == 0 ~ "M0", 
      M == 1 ~ "M1", 
      TRUE ~ NA_character_
    )),
    Sex = factor(
      dplyr::case_when(
        is.na(Sex) ~ NA_character_,
        Sex == 1 ~ "Female",
        Sex == 0 ~ "Male"
      ),
      levels = c("Female", "Male")
    ),
    Age = as.numeric(Age)
  ) |>
  assign_AJCC8_stage(
    t_col = "T_stage_comp",
    n_col = "N", 
    m_col = "M",
    age_col = "Age",
    out_col = "Stage"
  ) |>
  mutate(stage_bin = factor(case_when(
    Stage %in% c("I", "II") ~ "Early",
    Stage %in% c("III", "IV") ~ "Advanced",
    TRUE ~ NA_character_
  ), levels = c("Early", "Advanced")))
#- 0c.1.2: Create cleanup version
tumor_pathology <- tumor_pathology_full |>
  select(ID = Patient_ID, Variant, Sex, Age, T_stage_comp, Stage, stage_bin) |>
  arrange(Stage)