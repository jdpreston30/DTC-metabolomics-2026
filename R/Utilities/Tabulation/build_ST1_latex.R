#' Build Supplementary Table 1 LaTeX
#'
#' Converts a tibble to a landscape, multi-page LaTeX longtable with:
#' - Arial font, size 7pt
#' - Proper LaTeX formatting for column names (italics, subscripts)
#' - Multi-page spillover with repeated headers
#' - 0.5pt borders
#' - Hierarchical structure: Group (bold+underline) and Metabolite (indented plain text)
#'
#' @param data The prepared tibble with Display_Name column containing hierarchy
#' @param output_path Path to write the .tex file
#'
#' @return Invisibly returns the LaTeX code
#'
build_ST1_latex <- function(data, group_names, output_path = NULL) {
  
  library(kableExtra)
  library(dplyr)
  
  # Define formatted column names directly with larger font size (9pt vs 7pt table)
  formatted_col_names <- c(
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{\\underline{Group}:} Metabolite Name",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{Subgroup}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{\\textit{P} Value}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{\\textit{q} Value}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{log\\textsubscript{2}(FC)*}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{\\textit{m/z}}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{RT (s)}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{Column/ESI}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{Adduct}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{KEGG ID}",
    "\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{CID}"
  )
  
  # Build longtable with kableExtra
  latex_table <- kbl(
    data,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,
    col.names = formatted_col_names,
    align = c("l", "l", rep("c", ncol(data) - 2))
  ) |>
    kable_styling(
      latex_options = c("repeat_header"),
      font_size = 7
    )
  
  # Extract the raw LaTeX
  latex_code <- as.character(latex_table)
  
  # Remove any rowcolor commands to ensure no striping
  latex_code <- gsub("\\\\rowcolor\\{[^}]+\\}", "", latex_code)
  
  # Set all border thickness to 0.5pt
  latex_code <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_code, perl = TRUE)
  latex_code <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_code, perl = TRUE)
  latex_code <- gsub("\\\\bottomrule(?!\\[)", "\\\\bottomrule[0.5pt]", latex_code, perl = TRUE)
  
  # Remove (continued) text from repeated headers
  latex_code <- gsub("\\\\multicolumn\\{[0-9]+\\}\\{@\\{\\}l\\}\\{\\\\textit\\{\\(continued\\)\\}\\}\\\\\\\\\n", "", latex_code)
  
  # Remove all addlinespace commands (causes random row height variations)
  latex_code <- gsub("\\\\addlinespace(\\[[^]]*\\])?", "", latex_code)
  
  # Fix duplicate borders in endfirsthead structure - remove extra midrule between header and endfirsthead
  latex_code <- gsub("\\\\midrule\\[0\\.5pt\\]\\s*\\n\\s*\\\\endfirsthead\\s*\\n\\s*\\\\midrule\\[0\\.5pt\\]", "\\\\midrule[0.5pt]\n\\\\endfirsthead", latex_code)
  
  # Add bottomrule before \endfoot for page breaks
  latex_code <- gsub("\\n\\\\endfoot", "\n\\\\bottomrule[0.5pt]\n\\\\endfoot", latex_code)
  
  # Get number of columns for multicolumn in footer
  ncols <- ncol(data)
  
  # Add footnotes and abbreviations to \endlastfoot
  # Create blank row for spacing (using & separators for all columns)
  blank_row <- paste(c("", rep("", ncols - 1)), collapse = " & ")
  blank_row <- paste0(blank_row, " \\\\\\\\\n")
  
  footnote_text <- paste0(
    "\\\\bottomrule[0.5pt]\n",
    blank_row,
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont * Fold-change was calculated as advanced-stage/early-stage using untransformed peak area values.} \\\\\\\\\n",
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont \\\\textsuperscript{†} Metabolite annotated as a different adduct and/or using different column/ESI.} \\\\\\\\\n",
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont \\\\textsuperscript{‡} Annotated isomers: N-Acetyl-4-O-acetylneuraminate, N-Acetyl-7-O-acetylneuraminate} \\\\\\\\\n",
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont \\\\textsuperscript{§} Annotated isomer: 18-Hydroxyoleate} \\\\\\\\\n",
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont \\\\textsuperscript{‖} Annotated isomers: 11H-14,15-EETA, 12(R)-HPETE, 12(S)-HPETE, 15(S)-HPETE, 15H-11,12-EETA, 5(S)-HPETE, 8(R)-HPETE, Hepoxilin A3, Hepoxilin B3} \\\\\\\\[0.5em]\n",
    "\\\\multicolumn{", ncols, "}{@{}l@{}}{\\\\fontsize{7pt}{8.4pt}\\\\selectfont \\\\textit{FC}, fold-change; \\\\textit{RT}, retention time; \\\\textit{ESI}, electrospray ionization; \\\\textit{KEGG ID}, Kyoto Encyclopedia of Genes and Genomes ID; \\\\textit{CID}, PubChem Compound ID} \\\\\\\\\n",
    "\\\\endlastfoot"
  )
  
  # Replace the default endlastfoot with our custom footer
  latex_code <- gsub("\\\\bottomrule\\[0\\.5pt\\]\\n\\\\endlastfoot", footnote_text, latex_code)
  
  # Split into lines for processing
  lines <- strsplit(latex_code, "\n")[[1]]
  
  # Process data rows and headers
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Skip structural lines
    if (grepl("^\\\\(toprule|midrule|bottomrule|endfirsthead|endhead|endfoot|endlastfoot|addlinespace|multicolumn|rule|begin|end|caption)", line)) next
    
    # Skip header rows for data processing
    if (grepl("textbf\\{\\\\underline\\{Group\\}\\}", line)) next
    
    # Process data rows (have & separator)
    if (grepl(" & ", line)) {
      parts <- strsplit(line, " & ", fixed = TRUE)[[1]]
      if (length(parts) >= 1) {
        first_col <- parts[1]
        rest <- if(length(parts) > 1) paste(parts[-1], collapse = " & ") else ""
        
        # Trim and check if it's a known GROUP name
        first_col_trimmed <- trimws(first_col)
        
        if (first_col_trimmed %in% group_names || first_col_trimmed == "Redox Homeostasis (Continued)") {
          # GROUP: bold, underline, and larger font (9pt to match headers)
          first_col <- paste0("\\fontsize{9pt}{10.8pt}\\selectfont\\textbf{\\underline{", first_col_trimmed, "}}")
        } else if (nchar(first_col_trimmed) > 0 && !grepl("^\\\\", first_col) && first_col != "") {
          # Metabolite: add indent (0.3cm), superscript footnote symbols
          # Convert † ‡ § ‖ to \textsuperscript{}
          metabolite_text <- first_col_trimmed
          metabolite_text <- gsub("†", "\\\\textsuperscript{†}", metabolite_text)
          metabolite_text <- gsub("‡", "\\\\textsuperscript{‡}", metabolite_text)
          metabolite_text <- gsub("§", "\\\\textsuperscript{§}", metabolite_text)
          metabolite_text <- gsub("‖", "\\\\textsuperscript{‖}", metabolite_text)
          first_col <- paste0("\\hspace*{0.3cm}", metabolite_text)
        }
        # Empty/spacer rows pass through unchanged
        
        # Rebuild line
        if (length(parts) > 1) {
          lines[i] <- paste(first_col, rest, sep = " & ")
        } else {
          lines[i] <- first_col
        }
      }
    }
  }
  
  # Rejoin lines
  latex_code <- paste(lines, collapse = "\n")
  
  # Write to file if path provided
  if (!is.null(output_path)) {
    writeLines(latex_code, output_path)
    message("LaTeX table written to: ", output_path)
  }
  
  invisible(latex_code)
}
