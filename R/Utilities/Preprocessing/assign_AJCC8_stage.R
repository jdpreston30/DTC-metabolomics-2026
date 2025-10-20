#' Assign AJCC 8th Edition Overall Stage for Differentiated Thyroid Cancer
#'
#' Calculates Stage I–IV for papillary and follicular thyroid carcinoma
#' based on AJCC/UICC TNM 8th edition criteria.
#'
#' @param df Data frame with columns for T, N, M, and Age (in years).
#' @param t_col Character, column name containing T stage (e.g. "T_stage_comp").
#' @param n_col Character, column name containing N stage (e.g. "N").
#' @param m_col Character, column name containing M stage (e.g. "M").
#' @param age_col Character, column name with patient age (default = "age").
#' @param out_col Character, name of output stage column (default = "AJCC8_Stage").
#'
#' @details
#' For patients **<55 years**:
#' \itemize{
#'   \item Any T, Any N, M0 → Stage I
#'   \item Any T, Any N, M1 → Stage II
#' }
#' For patients **≥55 years**:
#' \itemize{
#'   \item T1–T2, N0/Nx, M0 → Stage I
#'   \item T1–T2 with N1 or T3 (any N), M0 → Stage II
#'   \item T4, Any N, M0 → Stage III
#'   \item Any T, Any N, M1 → Stage IV
#' }
#'
#' @return Input data frame with a new ordered factor column for AJCC Stage.
#' @export
assign_AJCC8_stage <- function(
    df,
    t_col = "T_stage_comp",
    n_col = "N",
    m_col = "M",
    age_col = "age",
    out_col = "AJCC8_Stage") {
  T_stage <- df[[t_col]]
  N_stage <- df[[n_col]]
  M_stage <- df[[m_col]]
  Age <- df[[age_col]]

  df[[out_col]] <- dplyr::case_when(
    # ---- Age <55 ----
    Age < 55 & M_stage == "M0" ~ "I",
    Age < 55 & M_stage == "M1" ~ "II",

    # ---- Age ≥55 ----
    Age >= 55 & M_stage == "M1" ~ "IV",
    Age >= 55 & T_stage %in% c("T4", "T4a", "T4b") & M_stage == "M0" ~ "III",
    Age >= 55 & T_stage %in% c("T1", "T2") & N_stage %in% c("N1", "N1a", "N1b") & M_stage == "M0" ~ "II",
    Age >= 55 & T_stage %in% c("T3") & M_stage == "M0" ~ "II",
    Age >= 55 & T_stage %in% c("T1", "T2") & N_stage %in% c("N0", "NX", NA) & M_stage == "M0" ~ "I",
    TRUE ~ NA_character_
  )

  df[[out_col]] <- factor(df[[out_col]],
    levels = c("I", "II", "III", "IV"),
    ordered = TRUE
  )
  return(df)
}
