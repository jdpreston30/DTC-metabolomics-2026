#* 6: Tables
#+ 6.1: Demographics Table
#- 6.1.1: Summarize Age
age <- demographics |>
  group_by(variant) |>
  summarise(Mean_Stdev = paste0(round(mean(age_collection), 1), " ± ", round(sd(age_collection), 1))) |>
  pivot_wider(names_from = variant, values_from = Mean_Stdev) |>
  mutate(Variable = "Age") |>
  bind_cols(
    tibble(Total = paste0(round(mean(demographics$age_collection), 1), " ± ", round(sd(demographics$age_collection), 1)))
  ) |>
  select(Variable, Follicular, FVPTC, Papillary, Total)
#- 6.1.2: Summarize Sex
sex <- demographics |>
  group_by(variant) |>
  summarise(n_female = sum(sex == "Female"), total = n(), .groups = "drop") |>
  mutate(
    percent = round((n_female / total) * 100),
    count_percent = paste0(n_female, " (", percent, "%)")
  ) |>
  select(variant, count_percent) |>
  pivot_wider(names_from = variant, values_from = count_percent) |>
  mutate(Variable = "Sex (Female)") |>
  bind_cols(
    demographics |>
      summarise(n_female = sum(sex == "Female"), total = n()) |>
      mutate(Total = paste0(n_female, " (", round((n_female / total) * 100), "%)")) |>
      select(Total)
  ) |>
  select(Variable, Follicular, FVPTC, Papillary, Total)
#- 6.1.3: Build Table
table_1 <- build_table_1(
  data = rbind(age, sex),
  export_path = "Outputs/Tables/T1.xlsx"
)