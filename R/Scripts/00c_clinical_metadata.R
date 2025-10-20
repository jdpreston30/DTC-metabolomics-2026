#* 0c: Clinical metadata processing
#+ 0c.1: Cleanup tumor pathology data; compute T stage; compute overall stage
tumor_pathology <- tumor_pathology_raw %>%
  select(-T) %>%
  assign_T_stage(ld_col = "LD", ete_col = "ETE", units = "cm", out_col = "T_stage_comp") %>%
  mutate(
    T_stage_comp = factor(T_stage_comp, levels = c("T1", "T2", "T3", "T4"), ordered = TRUE),
    N = case_when(
      N == 0 ~ "N0", 
      N == 1 ~ "N1", 
      TRUE ~ NA_character_
    ),
    M = case_when(
      M == 0 ~ "M0", 
      M == 1 ~ "M1", 
      TRUE ~ NA_character_
    ),
    Sex = factor(
      dplyr::case_when(
        is.na(Sex) ~ NA_character_,
        Sex == 1 ~ "Female",
        Sex == 0 ~ "Male"
      ),
      levels = c("Female", "Male")
    ),
    Age = as.numeric(Age)
  ) %>%
  assign_AJCC8_stage(
    t_col = "T_stage_comp",
    n_col = "N", 
    m_col = "M",
    age_col = "Age",
    out_col = "Stage"
  ) %>%
  select(Patient_ID, Sex, Age, ETE, LD, T_stage_comp, N, M, Stage) %>%
  arrange(Stage)
