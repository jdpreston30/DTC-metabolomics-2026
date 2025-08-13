#* 3: Combined Clustering analysis
  #+ 3.0: Setup
    conflicts_prefer(ggplot2::margin)
  #+ 3.1: Combined PCA analysi
    pca_variant_results <- make_PCA(
      data = UFT_meta"/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Manuscripts and Projects/Active Projects/Thyroid Metabolomics/DTC-metabolomics-2025/Raw_Data/Final/tumor_pathology.xlsx")
      #_Create heatmap and get results including clade analysis
      heatmap_results <- make_heatmap(
        data = UFT_metaboanalyst_log2,
        variant_colors = variant_colors,
        top_features = 500,
        variant_levels = c("FTC", "FV-PTC", "PTC"),
        n_clades = 2,
        output_filename = "Combined_heatmap.svg"
      )
  #+ 3.3: Combined PCA analysis with hierarchical clusters
    #- 3.3.1: Create modified dataset with clade assignments
      #_Create data with clade assignments instead of variants
      UFT_with_clades <- UFT_metaboanalyst_log2 %>%
        dplyr::left_join(
          heatmap_results$clade_df %>% dplyr::select(Patient_ID, Clade),
          by = "Patient_ID"
        ) %>%
        dplyr::mutate(
          Variant = factor(paste0("Clade_", Clade), levels = c("Clade_1", "Clade_2"))
        ) %>%
        dplyr::select(-Clade) %>%
        dplyr::filter(!is.na(Variant))
    #- 3.3.2: Run PCA with clade colorss
      pca_clades_results <- make_PCA(
        data = UFT_with_clades,
        output_filename = "Combined_pca_clades.svg",
        ellipse_colors = clade_colors
      )
  #+ 3.4: Create heatmaps
    #- 3.4.1: Prepare data for heatmap
      heatmap_data_T <- UFT_metaboanalyst_log2 %>%
        left_join(tumor_pathology %>% select(Patient_ID, T), by = "Patient_ID") %>%
        select(Patient_ID, Variant, T, everything()) %>%
        rename(T_stage = T) %>%
        mutate(T_stage = case_when(
          is.na(T_stage) ~ NA_character_,
          TRUE ~ paste("T", T_stage, sep = "")
        ))
    #- 3.4.2: Create heatmaps with different feature selections
        variance_full <- make_heatmap(
          heatmap_data_T,
          variant_colors = variant_colors,
          feature_selector = "variance",
          top_features = FALSE,
          annotate_t_stage = TRUE,
          T_stage_colors = T_stage_colors,
          output_filename = NULL,
          output_png_filename = "variance_full.png",
          png_width = 2000, png_height = 2000, png_res = 220
        )
        variance_2000 <- make_heatmap(
          heatmap_data_T,
          variant_colors = variant_colors,
          feature_selector = "variance",
          top_features = 2000,
          annotate_t_stage = TRUE,
          T_stage_colors = T_stage_colors,
          output_filename = NULL,
          output_png_filename = "variance_2000.png",
          png_width = 2000, png_height = 2000, png_res = 220
        )
        variance_500 <- make_heatmap(
          heatmap_data_T,
          variant_colors = variant_colors,
          feature_selector = "variance",
          top_features = 500,
          annotate_t_stage = TRUE,
          T_stage_colors = T_stage_colors,
          output_filename = NULL,
          output_png_filename = "variance_500.png",
          png_width = 2000, png_height = 2000, png_res = 220
        )
  #+ PERMANOVA
    variance_study <- UFT_metaboanalyst_log2 %>%
      left_join(tumor_pathology, by = "Patient_ID") %>%
      rename(T_stage = T) %>%
      mutate(T_stage = case_when(
        is.na(T_stage) ~ NA_character_,
        TRUE ~ paste("T", T_stage, sep = "")
      )) %>%
      select(Patient_ID, T_stage, Sex, Age, MFC, LD, LVI, ETE, everything(), -c(Cluster, Patient, Variant_Name, N, M, Grade))
# INPUT: variance_study tibble (samples x [metadata + 20k features]), already log2-transformed & imputed

# --- 1) Define metadata and build feature matrix
meta_vars <- c("T_stage", "Sex", "Age", "MFC", "LD", "LVI", "ETE", "Variant")

# Clean up metadata types (helps avoid silent coercion issues)
meta <- variance_study |>
  dplyr::select(dplyr::all_of(meta_vars)) |>
  dplyr::mutate(
    T_stage = as.factor(T_stage),
    Sex     = as.factor(Sex),
    MFC     = as.factor(MFC), # ensure factor (fixes rows like "8 MFC 0.0132 0.887")
    LD      = as.numeric(LD),
    LVI     = as.factor(LVI),
    ETE     = as.factor(ETE),
    Variant = as.factor(Variant),
    Age     = as.numeric(Age)
  )

# Build numeric feature matrix (drop Patient_ID and all metadata columns)
features_matrix <- variance_study |>
  dplyr::select(-Patient_ID, -dplyr::all_of(meta_vars)) |>
  as.data.frame()

# Keep only samples with complete metadata, and align features
keep <- stats::complete.cases(meta)
meta_cc <- meta[keep, , drop = FALSE]
features_cc <- features_matrix[keep, , drop = FALSE]

# Drop any all-NA or zero-variance features (paranoid-safety)
all_na <- vapply(features_cc, function(z) all(is.na(z)), logical(1))
if (any(all_na)) features_cc <- features_cc[, !all_na, drop = FALSE]

zv <- vapply(features_cc, function(z) stats::var(z, na.rm = TRUE), numeric(1)) == 0
if (any(zv)) features_cc <- features_cc[, !zv, drop = FALSE]

# --- 2) Select top 20% by variance (on log2 scale as-is)
feature_vars <- vapply(features_cc, function(z) stats::var(z, na.rm = TRUE), numeric(1))
n_keep <- ceiling(0.20 * ncol(features_cc))
top_idx <- order(feature_vars, decreasing = TRUE)[1:n_keep]
features_top <- as.matrix(features_cc[, top_idx, drop = FALSE])

# --- 3) Distance matrix (Euclidean on log2; skip scaling since already log2; add scale() if desired)
dist_mat <- vegan::vegdist(features_top, method = "euclidean")

# --- 4) Run PERMANOVA for EACH metadata variable separately (R2, p)
get_permanova_stats <- function(varname) {
  f <- stats::as.formula(paste("dist_mat ~", varname))
  res <- vegan::adonis2(f, data = meta_cc, permutations = 999)
  tibble::tibble(
    Variable = varname,
    R2       = as.numeric(res$R2[1]),
    p_value  = as.numeric(res$`Pr(>F)`[1])
  )
}

results_top20 <- dplyr::bind_rows(lapply(meta_vars, get_permanova_stats)) |>
  dplyr::arrange(dplyr::desc(R2)) %>%
  arrange(p_value)




  #+ 3.5: Run t-tests for pairwise variant comparisons
    #- 3.5.1: FV-PTC vs PTC comparison
      #_Extract variant assignments (excluding FTC)
      variant_assignments_FVPTC_PTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "FTC") %>%
        dplyr::select(Patient_ID, Variant)
      FVPTC_vs_PTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FTC"),
        group_assignments = variant_assignments_FVPTC_PTC,
        group_column = "Variant",
        output_filename = "FVPTC_vs_PTC.csv",
        group1_value = "PTC",
        group2_value = "FV-PTC"
      )
    #- 3.5.2: FTC vs PTC comparison
      #_Extract variant assignments (excluding FV-PTC)
      variant_assignments_FTC_PTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "FV-PTC") %>%
        dplyr::select(Patient_ID, Variant)
      FTC_vs_PTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "FV-PTC"),
        group_assignments = variant_assignments_FTC_PTC,
        group_column = "Variant",
        output_filename = "FTC_vs_PTC.csv",
        group1_value = "PTC",
        group2_value = "FTC"
      )
    #- 3.5.3: FTC vs FV-PTC comparison
      #_Extract variant assignments (excluding PTC)
      variant_assignments_FTC_FVPTC <- UFT_metaboanalyst_log2 %>%
        dplyr::filter(Variant != "PTC") %>%
        dplyr::select(Patient_ID, Variant)
      FTC_vs_FVPTC <- mummichog_ttests(
        data = UFT_metaboanalyst_log2 %>% dplyr::filter(Variant != "PTC"),
        group_assignments = variant_assignments_FTC_FVPTC,
        group_column = "Variant",
        output_filename = "FTC_vs_FVPTC.csv",
        group1_value = "FV-PTC",
        group2_value = "FTC"
      )
  #+ 3.6: Pathway Enrichment Plot
    #- 3.6.1: Import results from mummichog
      FVPTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_PTC") %>%
        mutate(Comparisons = "FVPTC_PTC")
      FTC_PTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FTC_PTC") %>%
        mutate(Comparisons = "FTC_PTC")
      FVPTC_FTC <- read_xlsx("Outputs/Mummichog Outputs/variant_mummichog.xlsx", sheet = "FVPTC_FTC") %>%
        mutate(Comparisons = "FVPTC_FTC")
    #- 3.6.2: Bind rows then filter to important variables
    variant_enrichment <- bind_rows(FVPTC_PTC, FTC_PTC, FVPTC_FTC) %>%
      mutate(enrichment_factor = Hits.sig / Expected) %>%
      select(Comparisons, pathway_name, p_fisher, enrichment_factor) %>%
      filter(p_fisher < 0.05) %>%
      mutate(
        # Label comparisons nicely
        Comparisons = dplyr::case_when(
          Comparisons == "FVPTC_PTC" ~ "PTC vs. FV-PTC",
          Comparisons == "FTC_PTC" ~ "FTC vs. PTC",
          Comparisons == "FVPTC_FTC" ~ "FTC vs. FV-PTC",
          TRUE ~ Comparisons
        ),
        Comparisons = factor(Comparisons, levels = c("FTC vs. FV-PTC", "FTC vs. PTC", "PTC vs. FV-PTC")),
        pathway_name = forcats::fct_reorder(pathway_name, enrichment_factor, .fun = max)
      ) %>%
      tidyr::complete(pathway_name, Comparisons) %>%
      mutate(enrichment_factor = pmin(enrichment_factor, 7)) %>%
      mutate(
        # First handle Beta/β
        pathway_name = stringr::str_replace_all(
          pathway_name,
          stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE),
          "β-Alanine"
        ),
        # Capitalization fixes
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bmetabolism\\b", TRUE), "Metabolism"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bleucine\\b", TRUE), "Leucine"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bisoleucine\\b", TRUE), "Isoleucine"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bdegradation\\b", TRUE), "Degradation"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bmannose\\b", TRUE), "Mannose"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bacid\\b", TRUE), "Acid"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\boxidation\\b", TRUE), "Oxidation"),
        pathway_name = stringr::str_replace_all(pathway_name, stringr::regex("\\bperoxisome\\b", TRUE), "Peroxisome")) %>%
      mutate(
        pathway_name = stringr::str_replace_all(
          pathway_name,
          stringr::regex("\\bretinol\\b", ignore_case = TRUE),
          "Retinol"
        )
      ) %>%
      mutate(
        pathway_name = factor(
          pathway_name,
          levels = filter(., Comparisons == "FTC vs. FV-PTC") %>%
            arrange(desc(enrichment_factor)) %>%
            pull(pathway_name) %>%
            unique()
        )
      )
    #- 3.6.3: Plot
      variant_enrichment_plot <- ggplot(
        variant_enrichment,
        aes(x = 0.5, y = 0.5, size = enrichment_factor, color = p_fisher)
      ) +
        # One dummy row per facet -> avoid the warning
        geom_tile(
          data = data.frame(x = 0.5, y = 0.5),
          aes(x = x, y = y),
          width = 1, height = 1,
          fill = "white", colour = "grey80", linewidth = 0.3,
          inherit.aes = FALSE
        ) +
        geom_point(
          alpha = 0.95, shape = 16, stroke = 0,
          na.rm = TRUE, show.legend = TRUE
        ) +
        facet_grid(
          rows = vars(pathway_name),
          cols = vars(Comparisons),
          switch = "y", drop = FALSE
        ) +
        coord_fixed(clip = "off") +
        scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +

        # Keep limits ascending; reverse legend order via guide
        scale_size_continuous(
          range = c(5, 20),
          limits = c(0, 7),
          breaks = c(7, 5, 3, 1), # Labels show in this order…
          name = "Enrichment factor",
          guide = guide_legend(reverse = TRUE) # …because we reverse the legend
        ) +

        # Keep p limits ascending; reverse colorbar via guide
        scale_color_gradient(
          low = "#0a2256", high = "#c3dbe9", # Dark (small p) -> light (large p)
          limits = c(0.01, 0.05),
          oob = scales::squish,
          name = NULL,
          guide = guide_colorbar(
            reverse = TRUE, # 0.01 at top, 0.05 at bottom
            barheight = unit(5, "cm"),
            barwidth = unit(0.9, "cm")
          )
        ) +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_family = "Arial") +
        theme(
          text = element_text(family = "Arial"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0, "pt"),
          panel.spacing.y = unit(0, "pt"),
          strip.placement = "outside",
          strip.text.x.top = element_text(
            angle = 0, vjust = 1,
            face = "bold", family = "Arial", size = 6
          ),
          strip.text.y.left = element_text(
            angle = 0, hjust = 1,
            face = "bold", family = "Arial", size = 14,
            margin = margin(r = 6)
          ),
          legend.title = element_text(size = 14, face = "bold", family = "Arial"),
          legend.text = element_text(size = 14, family = "Arial"),
          plot.margin = margin(t = 20, r = 40, b = 10, l = 40)
        ) + coord_cartesian(clip = "off") 
    #- 3.6.4: Export Plot
panel_size <- 0.3 # tweak this until happy

n_rows <- length(unique(variant_enrichment$pathway_name))
n_cols <- length(unique(variant_enrichment$Comparisons))

# Add space for legends, strip labels, margins
extra_width <- 3 # inches for legends on the right
extra_height <- 1 # inches for titles/margins

total_width <- n_cols * panel_size + extra_width
total_height <- n_rows * panel_size + extra_height

ggsave(
  "variant_enrichment_plot.svg",
  variant_enrichment_plot,
  width = total_width,
  height = total_height,
  units = "in"
)
