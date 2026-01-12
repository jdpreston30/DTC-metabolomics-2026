process_feature_table <- function(data) {
  data |>
    filter(!str_starts(Sample_ID, "nist|q4")) |>
    left_join(tumor_IDs, by = "Sample_ID") |>
    select(ID, everything(), -Sample_ID) |>
    arrange(ID) |>
    left_join(tumor_pathology, by = "ID") |>
    select(ID, colnames(tumor_pathology), everything())
}
