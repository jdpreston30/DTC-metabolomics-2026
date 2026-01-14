#* 7: Data Not Shown
#+ 7.1: Description of metabolite counts
#- 7.1.1: Count features in full untargeted dataset
{
  hilic_counts_full <- UFT_full |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_counts_full <- UFT_full |> 
    select(starts_with("C18")) |> 
    ncol()
  total_untargeted_features <- hilic_counts_full + c18_counts_full
}
#- 7.1.2: Count filtered features 
{
  hilic_count_filtered <- UFT_filtered |> 
    select(starts_with("HILIC")) |> 
    ncol()
  c18_count_filtered <- UFT_filtered |> 
    select(starts_with("C18")) |> 
    ncol()
  total_filtered_features <- hilic_count_filtered + c18_count_filtered
}
#- 7.1.3: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "FULL UNTARGETED DATASET FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "HILIC (positive mode):  ", format(hilic_counts_full, big.mark = ","), "\n",
    "C18 (negative mode):    ", format(c18_counts_full, big.mark = ","), "\n",
    "Total features:         ", format(total_untargeted_features, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
  cat(
    "\n", strrep("=", 60), "\n",
    "QC-FILTERED UNTARGETED DATASET FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "HILIC (positive mode):  ", format(hilic_count_filtered, big.mark = ","), "\n",
    "C18 (negative mode):    ", format(c18_count_filtered, big.mark = ","), "\n",
    "Total filtered features:", format(total_filtered_features, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
}
#+ 7.2: Volcano Plot Statistics
#- 7.2.1: Extract volcano plot results (assuming you have volcano analysis results)
{
  sig_up_fc <- volcano_data$volcano_data |>
    filter(p_value < 0.05 & log2_fc > log2(1.5)) |>
    nrow()
  sig_down_fc <- volcano_data$volcano_data  |>
    filter(p_value < 0.05 & log2_fc < -log2(1.5)) |>
    nrow()
}
#- 7.2.2: Narrative printout
{
  cat(
    "\n", strrep("=", 60), "\n",
    "VOLCANO PLOT SIGNIFICANT FEATURE COUNTS\n",
    strrep("=", 60), "\n\n",
    "Significantly upregulated features (FC > 1.5, p < 0.05):   ", format(sig_up_fc, big.mark = ","), "\n",
    "Significantly downregulated features (FC < -1.5, p < 0.05): ", format(sig_down_fc, big.mark = ","), "\n",
    strrep("=", 60), "\n\n"
  )
}
