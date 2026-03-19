#' Assign 2025 ATA Risk of Recurrence (Nuanced Histology-Specific)
#' 
#' @param df The tumor_pathology_table tibble
assign_ATA_2025 <- function(df) {
  df %>%
    mutate(ATA_Risk = case_when(
      # --- 1. UNIVERSAL HIGH RISK (>30%) ---
      `M Category` == "M1" ~ "High",
      `T Category` == "T4" ~ "High",
      `Extrathyroidal Extension` == "Extensive" ~ "High",
      
      # --- 2. PTC SPECIFIC RULES ---
      Variant == "PTC" & (
        `T Category` == "T3" |                          # pT3a (>4cm)
        `Lymphovascular Invasion (Positive)` == "Y" |   # Any LVI in PTC
        `N Category` == "N1"                            # pN1 is generally Int-High for PTC
        # Bilateral Multifocality would go here if you have laterality data
      ) ~ "Intermediate-High",
      
      Variant == "PTC" & (
        `Extrathyroidal Extension` == "Minimal" |
        `Multifocality (Positive)` == "Yes"
      ) ~ "Low-Intermediate",

      # --- 3. FOLLICULAR-PATTERNED RULES (FTC & IEFVPTC) ---
      Variant %in% c("FTC", "IEFVPTC") & (
        `T Category` == "T3" |
        # In 2025, extensive LVI (>=4) is Int-High. 
        # Since we don't have counts, we classify 'Y' as Int-High to be safe.
        `Lymphovascular Invasion (Positive)` == "Y"
      ) ~ "Intermediate-High",
      
      # For FTC, multifocality isn't an ATA risk driver, 
      # but Minimal ETE or focal LVI would be Low-Intermediate.
      Variant %in% c("FTC", "IEFVPTC") & `Extrathyroidal Extension` == "Minimal" ~ "Low-Intermediate",

      # --- 4. DEFAULT TO LOW ---
      TRUE ~ "Low"
    )) %>%
    mutate(ATA_Risk = factor(ATA_Risk, 
                            levels = c("Low", "Low-Intermediate", "Intermediate-High", "High"), 
                            ordered = TRUE))
}
