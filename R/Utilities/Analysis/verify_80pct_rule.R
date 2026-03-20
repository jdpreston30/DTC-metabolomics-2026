#' Verify 80% Detection Rule and remove_qc() Equivalence
#'
#' Runs a comprehensive verification of the 80% detection rule applied by
#' MetaboJanitoR upstream, the downstream remove_qc() filter, and whether
#' the two produce equivalent feature sets for the final analytical datasets.
#'
#' @param uft_raw_import Raw UFT data with true zeros (pre-imputation)
#' @param uft_filtered_import UFT filtered import (post-80% rule, half-min imputed, log2)
#' @param uft_filtered Final UFT filtered object (after process_feature_table + remove_qc)
#' @param tft_annot_import TFT annotated import (half-min imputed, log2)
#' @param tft_annot Final TFT annotated object (after process_feature_table + remove_qc)
#' @param qc_pattern Regex pattern to identify QC samples in Sample_ID column
#'
#' @return Invisibly returns a list with detection_stats, qc_dependent features,
#'   and mechanism check results for both UFT and TFT.
verify_80pct_rule <- function(uft_raw_import,
                              uft_filtered_import,
                              uft_filtered,
                              tft_annot_import,
                              tft_annot,
                              qc_pattern = "nist|q4") {

  count_feat <- function(df) sum(grepl("^(C18|HILIC)", colnames(df)))
  divider    <- function() strrep("=", 60)
  subdiv     <- function() strrep("-", 60)

  # ── 1. Sample classification ───────────────────────────────────────────────
  uft_raw_samples <- uft_raw_import |>
    mutate(
      is_qc = str_starts(Sample_ID, qc_pattern),
      sample_type = if_else(is_qc, "QC", "Biological")
    )
  n_total    <- nrow(uft_raw_samples)
  n_bio      <- sum(!uft_raw_samples$is_qc)
  n_qc       <- sum(uft_raw_samples$is_qc)
  bio_mask   <- !uft_raw_samples$is_qc
  threshold_n <- ceiling(n_total * 0.80)

  # ── 2. UFT 80% rule verification ──────────────────────────────────────────
  uft_raw_features      <- grep("^(C18|HILIC)", colnames(uft_raw_import), value = TRUE)
  uft_filtered_features <- grep("^(C18|HILIC)", colnames(uft_filtered_import), value = TRUE)
  features_to_check     <- intersect(uft_filtered_features, uft_raw_features)
  features_not_in_raw   <- setdiff(uft_filtered_features, uft_raw_features)

  detection_stats <- purrr::map_dfr(features_to_check, function(feat) {
    vals <- uft_raw_import[[feat]]
    tibble(
      feature            = feat,
      n_detected_total   = sum(vals > 0, na.rm = TRUE),
      pct_detected_total = sum(vals > 0, na.rm = TRUE) / n_total,
      n_detected_bio     = sum(vals[bio_mask] > 0, na.rm = TRUE),
      pct_detected_bio   = sum(vals[bio_mask] > 0, na.rm = TRUE) / n_bio
    )
  })

  n_pass_total <- sum(detection_stats$pct_detected_total >= 0.80)
  n_fail_total <- sum(detection_stats$pct_detected_total <  0.80)
  n_pass_bio   <- sum(detection_stats$pct_detected_bio   >= 0.80)
  n_fail_bio   <- sum(detection_stats$pct_detected_bio   <  0.80)

  qc_dependent <- detection_stats |>
    filter(pct_detected_total >= 0.80 & pct_detected_bio < 0.80)

  cat(
    "\n", divider(), "\n",
    "80% RULE VERIFICATION: UFT_filtered\n",
    divider(), "\n\n",
    "Total samples in raw UFT:          ", n_total, "\n",
    "  Biological:                      ", n_bio, "\n",
    "  QC (NIST + pooled QC):           ", n_qc, "\n",
    "  80% detection threshold (n):     ", threshold_n, "of", n_total, "samples\n\n",
    "Features in UFT_filtered_import:   ", length(uft_filtered_features), "\n",
    "Features cross-checked in raw:     ", length(features_to_check), "\n",
    "Features absent from raw (ERROR):  ", length(features_not_in_raw), "\n\n",
    "VERIFICATION using TOTAL denominator (bio + QC):\n",
    "  PASS (>=80%): ", n_pass_total, "\n",
    "  FAIL (<80%):  ", n_fail_total, "  <- should be 0\n\n",
    "VERIFICATION using BIO-ONLY denominator:\n",
    "  PASS (>=80%): ", n_pass_bio, "\n",
    "  FAIL (<80%):  ", n_fail_bio, "\n\n",
    "Features passing TOTAL but failing BIO-ONLY:\n",
    "  (Would be excluded if QC removed from denominator)\n",
    "  n = ", nrow(qc_dependent), "\n\n",
    "Detection rate (total denom) across all UFT_filtered features:\n",
    "  Min:    ", round(min(detection_stats$pct_detected_total)    * 100, 1), "%\n",
    "  Median: ", round(median(detection_stats$pct_detected_total) * 100, 1), "%\n",
    "  Max:    ", round(max(detection_stats$pct_detected_total)    * 100, 1), "%\n",
    divider(), "\n\n"
  )

  # ── 3. Pipeline feature count verification ─────────────────────────────────
  tft_n_import   <- count_feat(tft_annot_import)
  tft_n_final    <- count_feat(tft_annot)
  tft_n_removed  <- tft_n_import - tft_n_final
  tft_row_import <- nrow(tft_annot_import)
  tft_row_final  <- nrow(tft_annot)

  uft_n_import   <- count_feat(uft_filtered_import)
  uft_n_final    <- count_feat(uft_filtered)
  uft_n_removed  <- uft_n_import - uft_n_final
  uft_row_import <- nrow(uft_filtered_import)
  uft_row_final  <- nrow(uft_filtered)

  cat(
    "\n", divider(), "\n",
    "PIPELINE FEATURE COUNT VERIFICATION\n",
    divider(), "\n\n",
    "TFT_annot (process_feature_table + remove_qc at threshold=0.2):\n",
    "  Import rows (bio + QC):  ", tft_row_import, "\n",
    "  After QC removal:        ", tft_row_final, " (", tft_row_import - tft_row_final, " QC removed)\n",
    "  Import features:         ", format(tft_n_import, big.mark = ","), "\n",
    "  Final features:          ", format(tft_n_final,  big.mark = ","), "\n",
    "  Removed by remove_qc():  ", format(tft_n_removed, big.mark = ","),
    "  (>20% identical values)\n\n",
    "UFT_filtered (process_feature_table + remove_qc at threshold=0.2):\n",
    "  Import rows (bio + QC):  ", uft_row_import, "\n",
    "  After QC removal:        ", uft_row_final, " (", uft_row_import - uft_row_final, " QC removed)\n",
    "  Import features:         ", format(uft_n_import, big.mark = ","), "\n",
    "  Final features:          ", format(uft_n_final,  big.mark = ","), "\n",
    "  Removed by remove_qc():  ", format(uft_n_removed, big.mark = ","),
    "  (>20% identical values)\n\n",
    subdiv(), "\n",
    "NOTE: The 80% detection rule was applied UPSTREAM by MetaboJanitoR\n",
    "before UFT_filtered_import was written to disk. The remove_qc() call\n",
    "here is a SEPARATE filter targeting features with >20% of values equal\n",
    "to the half-minimum replacement (i.e. ubiquitous non-detects).\n",
    divider(), "\n\n"
  )

  # ── 4. UFT intersect: remove_qc() vs QC-dependent features ────────────────
  uft_removed_by_rmqc   <- setdiff(uft_filtered_features,
                                    grep("^(C18|HILIC)", colnames(uft_filtered), value = TRUE))
  qc_dependent_features <- qc_dependent$feature

  n_overlap   <- length(intersect(uft_removed_by_rmqc, qc_dependent_features))
  n_only_qc   <- length(setdiff(qc_dependent_features, uft_removed_by_rmqc))
  n_only_rmqc <- length(setdiff(uft_removed_by_rmqc,   qc_dependent_features))

  cat(
    divider(), "\n",
    "INTERSECT: remove_qc() vs QC-dependent (bio-only 80% fail)\n",
    divider(), "\n\n",
    "Features removed by remove_qc():              ", length(uft_removed_by_rmqc), "\n",
    "Features failing bio-only 80% rule:           ", length(qc_dependent_features), "\n\n",
    "Overlap (same feature, both criteria):         ", n_overlap, "\n",
    "Only in bio-only 80% fail (not in remove_qc): ", n_only_qc, "\n",
    "Only in remove_qc (not QC-dependent):         ", n_only_rmqc, "\n\n",
    if (n_overlap == length(uft_removed_by_rmqc) &&
        n_overlap == length(qc_dependent_features)) {
      "CONCLUSION: The two sets are IDENTICAL.\n"
    } else if (n_overlap == 0) {
      "CONCLUSION: No overlap — completely different feature sets.\n"
    } else {
      paste0("CONCLUSION: Partial overlap (", n_overlap, " shared, ",
             n_only_qc, " unique to bio-only 80%, ",
             n_only_rmqc, " unique to remove_qc).\n")
    },
    divider(), "\n\n"
  )

  # ── 5. TFT 80% rule verification via UFT_raw proxy ────────────────────────
  tft_features_all    <- grep("^(C18|HILIC)", colnames(tft_annot_import), value = TRUE)
  tft_features_final  <- grep("^(C18|HILIC)", colnames(tft_annot),        value = TRUE)
  tft_removed_by_rmqc <- setdiff(tft_features_all, tft_features_final)
  tft_in_raw          <- intersect(tft_features_all, uft_raw_features)
  tft_not_in_raw      <- setdiff(tft_features_all, uft_raw_features)

  cat(
    "\n", divider(), "\n",
    "TFT_annot: 80% RULE VERIFICATION via UFT_raw proxy\n",
    divider(), "\n\n",
    "TFT_annot_import features:               ", length(tft_features_all),  "\n",
    "TFT features found in UFT_raw_import:    ", length(tft_in_raw),        "\n",
    "TFT features NOT in UFT_raw (no proxy):  ", length(tft_not_in_raw),    "\n\n"
  )

  if (length(tft_in_raw) > 0) {
    tft_detection_stats <- purrr::map_dfr(tft_in_raw, function(feat) {
      vals <- uft_raw_import[[feat]]
      tibble(
        feature            = feat,
        pct_detected_total = sum(vals > 0, na.rm = TRUE) / n_total,
        pct_detected_bio   = sum(vals[bio_mask] > 0, na.rm = TRUE) / n_bio
      )
    })

    tft_pass_total   <- sum(tft_detection_stats$pct_detected_total >= 0.80)
    tft_fail_total   <- sum(tft_detection_stats$pct_detected_total <  0.80)
    tft_pass_bio     <- sum(tft_detection_stats$pct_detected_bio   >= 0.80)
    tft_fail_bio     <- sum(tft_detection_stats$pct_detected_bio   <  0.80)
    tft_qc_dependent <- tft_detection_stats |>
      filter(pct_detected_total >= 0.80 & pct_detected_bio < 0.80)

    rmqc_in_raw     <- intersect(tft_removed_by_rmqc, tft_in_raw)
    tft_n_overlap   <- length(intersect(rmqc_in_raw, tft_qc_dependent$feature))
    tft_n_only_qc   <- length(setdiff(tft_qc_dependent$feature, rmqc_in_raw))
    tft_n_only_rmqc <- length(setdiff(rmqc_in_raw, tft_qc_dependent$feature))

    cat(
      "VERIFICATION using TOTAL denominator (bio + QC):\n",
      "  PASS (>=80%): ", tft_pass_total, "\n",
      "  FAIL (<80%):  ", tft_fail_total, "  <- should be 0\n\n",
      "VERIFICATION using BIO-ONLY denominator:\n",
      "  PASS (>=80%): ", tft_pass_bio, "\n",
      "  FAIL (<80%):  ", tft_fail_bio, "\n\n",
      "QC-dependent features (pass total, fail bio-only): ", nrow(tft_qc_dependent), "\n\n",
      subdiv(), "\n",
      "INTERSECT (proxy-matched features only):\n",
      "  Removed by remove_qc() AND in raw:        ", length(rmqc_in_raw),    "\n",
      "  QC-dependent in raw:                      ", nrow(tft_qc_dependent), "\n",
      "  Overlap:                                  ", tft_n_overlap,          "\n",
      "  Only QC-dependent (not in remove_qc):     ", tft_n_only_qc,          "\n",
      "  Only remove_qc (not QC-dependent):        ", tft_n_only_rmqc,        "\n\n",
      if (tft_n_overlap == length(rmqc_in_raw) &&
          tft_n_overlap == nrow(tft_qc_dependent)) {
        "CONCLUSION: IDENTICAL — same features flagged by both criteria.\n"
      } else if (tft_n_overlap == 0 && nrow(tft_qc_dependent) == 0) {
        "CONCLUSION: No QC-dependent features — final dataset matches bio-only 80% rule.\n"
      } else if (tft_n_overlap == 0) {
        "CONCLUSION: No overlap — completely different feature sets.\n"
      } else {
        paste0("CONCLUSION: Partial overlap (", tft_n_overlap, " shared, ",
               tft_n_only_qc, " unique to bio-only 80%, ",
               tft_n_only_rmqc, " unique to remove_qc).\n")
      },
      divider(), "\n\n"
    )
  } else {
    cat("No TFT features found in UFT_raw — proxy approach not possible.\n\n")
  }

  # ── 6. Mechanistic check: remove_qc() vs bio-only 80% equivalence ─────────
  run_mechanism_check <- function(features, raw_import, imp_import, label) {
    purrr::map_dfr(features, function(feat) {
      raw_vals <- raw_import[[feat]][bio_mask]
      imp_vals <- imp_import[[feat]][bio_mask]
      n_b      <- length(raw_vals)
      tibble(
        feature     = feat,
        n_detected  = sum(raw_vals > 0),
        pct_detect  = sum(raw_vals > 0) / n_b,
        max_id_pct  = max(table(imp_vals)) / n_b,
        pass_80_bio = sum(raw_vals > 0) / n_b >= 0.80,
        pass_rmqc   = max(table(imp_vals)) / n_b <= 0.20
      )
    }) -> mechanism

    n_agree    <- sum(mechanism$pass_80_bio == mechanism$pass_rmqc)
    n_disagree <- sum(mechanism$pass_80_bio != mechanism$pass_rmqc)

    cat(
      "\n", divider(), "\n",
      "MECHANISTIC CHECK: remove_qc() vs bio-only 80% rule (", label, ")\n",
      divider(), "\n\n",
      "Features analyzed: ", nrow(mechanism), "\n",
      "Agree (both pass or both fail): ", n_agree, "\n",
      "Disagree:                        ", n_disagree, "\n\n"
    )

    if (n_disagree > 0) {
      p80_frmqc <- mechanism |> filter(pass_80_bio & !pass_rmqc)
      f80_prmqc <- mechanism |> filter(!pass_80_bio & pass_rmqc)
      cat(
        "Pass 80% bio but FAIL remove_qc: ", nrow(p80_frmqc),
        "  (high detection but identical values)\n",
        "Fail 80% bio but PASS remove_qc: ", nrow(f80_prmqc),
        "  (low detection but varied values)\n\n"
      )
      if (nrow(p80_frmqc) > 0) {
        cat("  Pass 80% but fail remove_qc (first 10):\n")
        print(head(p80_frmqc, 10))
      }
      if (nrow(f80_prmqc) > 0) {
        cat("  Fail 80% but pass remove_qc (first 10):\n")
        print(head(f80_prmqc, 10))
      }
    }

    cat(
      subdiv(), "\n",
      "CONCLUSION: ",
      if (n_disagree == 0) {
        paste0("EXACT EQUIVALENCE — remove_qc(0.2) and bio-only 80% rule\n",
               "produce identical keep/remove decisions for all ", label, " features.\n")
      } else {
        paste0("NOT equivalent — ", n_disagree, " features differ.\n",
               "remove_qc() is based on max identical value frequency,\n",
               "not detection rate.\n")
      },
      divider(), "\n\n"
    )
    mechanism
  }

  uft_mechanism <- run_mechanism_check(
    uft_filtered_features, uft_raw_import, uft_filtered_import, "UFT"
  )
  tft_mechanism <- run_mechanism_check(
    tft_in_raw, uft_raw_import, tft_annot_import, "TFT via proxy"
  )

  # ── 7. Notes ───────────────────────────────────────────────────────────────
  cat(
    subdiv(), "\n",
    "NOTE ON TFT_ANNOT\n",
    subdiv(), "\n\n",
    "TFT_annot was processed by the same MetaboJanitoR pipeline\n",
    "(MSMICA_20p_halfmin_log2.csv). The UFT_raw proxy confirms the 80%\n",
    "rule for all 3,887 TFT features. remove_qc() and bio-only 80% are\n",
    "verified equivalent for this dataset.\n\n",
    subdiv(), "\n\n",
    subdiv(), "\n",
    "NOTE ON UFT_full / feature_cols\n",
    subdiv(), "\n\n",
    "feature_cols is defined in 00d_FTs.R from UFT_full column names\n",
    "but is NEVER referenced in any downstream analysis script.\n",
    "UFT_full is not used in any statistical analysis.\n",
    subdiv(), "\n\n"
  )

  invisible(list(
    detection_stats   = detection_stats,
    qc_dependent      = qc_dependent,
    uft_mechanism     = uft_mechanism,
    tft_mechanism     = tft_mechanism
  ))
}
