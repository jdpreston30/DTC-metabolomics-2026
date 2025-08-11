#* 3: Clustering analysis
  #+ 3.1: HILIC PCA analysis
    #- 3.1.1: Data preparation
      df <- as.data.frame(UFT_HILIC_metaboanalyst_log2)
      #_Choose class column: prefer 'Variant' if present; else use the 2nd column
      cls_col <- if ("Variant" %in% names(df)) "Variant" else names(df)[2]
      #_X = numeric features (drop first two columns), Y = class
      X <- df[, -c(1, 2), drop = FALSE]
      #_Coerce to numeric (safely)
      X[] <- lapply(X, function(v) suppressWarnings(as.numeric(v)))
      #_Simple NA guard: feature-wise median impute (optional but robust)
      if (anyNA(X)) {
        X[] <- lapply(X, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
      Y <- factor(df[[cls_col]])
    #- 3.1.2: Perform PCA
      pca <- stats::prcomp(X, center = TRUE, scale. = TRUE)
      scores <- pca$x[, 1:2, drop = FALSE]
      explained <- round((pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100)
    #- 3.1.3: Prepare plot data
      scores_df <- data.frame(
        Comp1 = scores[, 1],
        Comp2 = scores[, 2],
        Class = Y
      )
      #_Define colors for each variant
      ellipse_colors <- c("PTC" = "#DF8D0A", "FV-PTC" = "#23744E", "FTC" = "#194992")
      point_colors <- ellipse_colors
    #- 3.1.4: Create PCA plot
      hilic_pca <- ggplot2::ggplot(scores_df, ggplot2::aes(x = Comp1, y = Comp2, color = Class)) +
        ggplot2::geom_point(
          size = 3, shape = 21, stroke = 0.8,
          fill = point_colors[as.character(scores_df$Class)]
        ) +
        ggplot2::stat_ellipse(
          geom = "polygon", ggplot2::aes(fill = Class),
          alpha = 0.3, color = NA
        ) +
        ggplot2::scale_color_manual(values = point_colors, drop = FALSE) +
        ggplot2::scale_fill_manual(values = ellipse_colors, drop = FALSE) +
        ggplot2::theme_minimal(base_family = "Arial") +
        ggplot2::labs(
          x = paste0("PC1 (", explained[1], "%)"),
          y = paste0("PC2 (", explained[2], "%)")
        ) +
        ggplot2::theme(
          axis.title = ggplot2::element_text(size = 25, face = "bold"),
          axis.text = ggplot2::element_text(size = 22, face = "bold", color = "black"),
          legend.position = "none",
          panel.grid.major = ggplot2::element_line(color = "gray80", linewidth = 0.8, linetype = "solid"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 3.2),
          panel.background = ggplot2::element_blank()
        )
    #- 3.1.5: Export PCA plot
      ggplot2::ggsave(
        filename = "Outputs/HILIC_pca.svg",
        plot = hilic_pca,
        device = "svg",
        width = 8,
        height = 8,
        units = "in",
        dpi = 600
      )
  #+ 3.2: K-means clustering analysis
    #- 3.2.1: Prepare data for k-means
      #_Load UFT data and remove variant column
      kmeans_df <- as.data.frame(UFT_metaboanalyst_log2)
      patient_ids <- kmeans_df$Patient_ID
      #_Remove Patient_ID and Variant columns for clustering
      X_kmeans <- kmeans_df[, !names(kmeans_df) %in% c("Patient_ID", "Variant"), drop = FALSE]
      #_Ensure all features are numeric
      X_kmeans[] <- lapply(X_kmeans, function(v) suppressWarnings(as.numeric(v)))
      #_Handle any remaining NAs with median imputation
      if (anyNA(X_kmeans)) {
        X_kmeans[] <- lapply(X_kmeans, function(v) {
          v[is.na(v)] <- stats::median(v, na.rm = TRUE)
          v
        })
      }
    #- 3.2.2: Perform k-means clustering (k=2)
      set.seed(2025)  # For reproducibility
      kmeans_result <- stats::kmeans(X_kmeans, centers = 2, nstart = 25)
    #- 3.2.3: Extract cluster assignments
      cluster_assignments <- data.frame(
        Patient_ID = patient_ids,
        Cluster = paste0("Cluster_", kmeans_result$cluster)
      ) %>%
        arrange(Cluster)
      #_Display cluster assignments
      print("K-means clustering results (k=2):")
      print(cluster_assignments)
      #_Summary by cluster
      cluster_summary <- table(cluster_assignments$Cluster)
      print("Patients per cluster:")
      print(cluster_summary)
  #+ 3.3: T-tests between clusters
    #- 3.3.1: Prepare data for t-tests
      #_Add cluster assignments to feature data
      ttest_data <- kmeans_df %>%
        dplyr::mutate(Cluster = kmeans_result$cluster) %>%
        dplyr::select(-Patient_ID, -Variant)
    #- 3.3.2: Perform t-tests for each feature
      #_Get feature names (exclude Cluster column)
      feature_names <- names(ttest_data)[names(ttest_data) != "Cluster"]
      #_Initialize results list
      ttest_results <- list()
      #_Perform t-test for each feature
      for (feature in feature_names) {
        cluster1_values <- ttest_data[ttest_data$Cluster == 1, feature]
        cluster2_values <- ttest_data[ttest_data$Cluster == 2, feature]
        #_Perform t-test
        test_result <- t.test(cluster1_values, cluster2_values)
        ttest_results[[feature]] <- test_result$p.value
      }
    #- 3.3.3: Parse feature names and create results tibble
      results_tibble <- tibble::tibble(
        Feature = names(ttest_results),
        p.value = unlist(ttest_results)
      ) %>%
        dplyr::mutate(
          #_Extract mode and convert: HILIC -> pos, C18 -> neg
          mode = dplyr::case_when(
            stringr::str_starts(Feature, "HILIC") ~ "positive",
            stringr::str_starts(Feature, "C18") ~ "negative",
            TRUE ~ NA_character_
          ),
          #_Extract m.z (first number after underscore)
          m.z = stringr::str_extract(Feature, "(?<=_)[0-9.]+"),
          #_Extract rt (second number - after second underscore)
          r.t = stringr::str_extract(Feature, "_[0-9.]+_([0-9.]+)") %>%
            stringr::str_extract("[0-9.]+$")
        ) %>%
        dplyr::select(m.z, p.value, mode, r.t) %>%
        dplyr::mutate(
          m.z = as.numeric(m.z),
          r.t = as.numeric(r.t)
        )
    #- 3.3.4: Export results
      write_csv(results_tibble, "Outputs/cluster_ttest_results.csv")
      print("T-test results exported to Outputs/cluster_ttest_results.csv")
      print(paste("Total features tested:", nrow(results_tibble)))
      print(head(results_tibble, 10))
    #- 3.3.4: Export results
      readr::write_csv(results_tibble, "Outputs/cluster_ttest_results.csv")
  #+ 3.4: Heatmap with HILIC Untargeted
# 1) Split groups and matrix
dat <- UFT_HILIC_metaboanalyst_log2 %>% select(-Patient_ID)
group <- factor(dat$Variant, levels = c("PTC", "FV-PTC", "FTC")) # lock desired order
X <- as.matrix(dat %>% select(-Variant)) # rows = samples, cols = features
sample_ids <- make.names(UFT_HILIC_metaboanalyst_log2$Patient_ID, unique = TRUE)

# 2) ANOVA p-values per feature (over group)
pvals <- apply(X, 2, function(x) {
  fit <- aov(x ~ group)
  summary(fit)[[1]][["Pr(>F)"]][1]
})

# 3) Top 2500 (or all if fewer)
k <- min(2500, ncol(X))
top_idx <- order(pvals)[1:k]
X_top <- X[, top_idx, drop = FALSE]

# 4) Drop zero-variance features
nzv <- apply(X_top, 2, sd, na.rm = TRUE) > 0
X_top <- X_top[, nzv, drop = FALSE]

# 5) Names + transpose (rows=features, cols=samples)
rownames(X_top) <- sample_ids # samples
colnames(X_top) <- colnames(UFT_HILIC_metaboanalyst_log2)[-(1:2)][top_idx][nzv] # features
M <- t(X_top)

# 6) Column annotation (must match columns of M)
ann_col <- data.frame(Variant = group)
rownames(ann_col) <- sample_ids
ann_col <- ann_col[colnames(M), , drop = FALSE]
# 7) Custom colors for Variant
ann_colors <- list(
  Variant = c(
    "PTC"    = "#DF8D0A",
    "FV-PTC" = "#23744E",
    "FTC"    = "#194992"
  )
)

# 8) Heatmap
pheatmap(
  M,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(255),
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize = 10,
  na_col = "#DDDDDD"
)
