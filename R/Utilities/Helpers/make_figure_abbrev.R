#' Create Figure Abbreviation Legend Document
#'
#' Generates a Word document with abbreviation legends for each figure.
#' Abbreviations are italicized, alphabetized within each figure, and formatted
#' as "abbrev, longform; abbrev, longform; ..."
#'
#' @param abbreviation_tibble A tibble with columns: abbreviation, longform, fig1, fig2, fig3, fig4, etc.
#'   Figure columns should contain "Y" for figures where the abbreviation appears, NA otherwise.
#' @param output_path Path to save the Word document. Default: "Outputs/Figures/Final/figure_abbreviations.docx"
#'
#' @return Invisibly returns the path to the created document
#'
#' @examples
#' \dontrun{
#' abbreviation_tibble <- tibble::tibble(
#'   abbreviation = c("Ala", "Asp", "FA"),
#'   longform = c("alanine", "aspartate", "fatty acid"),
#'   fig1 = c("Y", "Y", "Y"),
#'   fig2 = c(NA, NA, "Y")
#' )
#' make_figure_abbrev(abbreviation_tibble)
#' }
#'
#' @export
make_figure_abbrev <- function(abbreviation_tibble, 
                                output_path = "Outputs/Figures/Final/figure_abbreviations.docx") {
  
  # Load required packages
  library(officer)
  library(dplyr)
  library(stringr)
  
  # Identify figure columns (any column starting with "fig")
  fig_cols <- names(abbreviation_tibble)[grepl("^fig", names(abbreviation_tibble), ignore.case = TRUE)]
  
  # Create a new Word document
  doc <- read_docx()
  
  # Process each figure
  for (fig_col in fig_cols) {
    # Extract abbreviations for this figure (where value is "Y")
    fig_abbrevs <- abbreviation_tibble |>
      filter(.data[[fig_col]] == "Y") |>
      arrange(abbreviation) |>
      select(abbreviation, longform)
    
    # Skip if no abbreviations for this figure
    if (nrow(fig_abbrevs) == 0) next
    
    # Extract figure number from column name (e.g., "fig1" -> "1")
    fig_num <- str_extract(fig_col, "\\d+")
    
    # Create the figure label with text properties
    fig_label <- paste0("Fig", fig_num, " abbrev:")
    fp_normal <- fp_text(font.family = "Arial", font.size = 11, italic = FALSE, bold = FALSE)
    
    # Add figure label as a formatted paragraph
    doc <- doc |>
      body_add_fpar(
        fpar(ftext(fig_label, prop = fp_normal))
      )
    
    # Build the abbreviation string with proper formatting
    # We'll add this as a single paragraph with mixed formatting
    abbrev_entries <- vector("list", nrow(fig_abbrevs))
    
    for (i in seq_len(nrow(fig_abbrevs))) {
      # Create each entry: "abbrev, longform"
      abbrev_entries[[i]] <- list(
        abbrev = fig_abbrevs$abbreviation[i],
        longform = fig_abbrevs$longform[i]
      )
    }
    
    # Build plain text version for console output
    plain_text_parts <- character(nrow(fig_abbrevs))
    for (i in seq_along(abbrev_entries)) {
      plain_text_parts[i] <- paste0(abbrev_entries[[i]]$abbrev, ", ", abbrev_entries[[i]]$longform)
    }
    plain_text <- paste(plain_text_parts, collapse = "; ")
    plain_text <- paste0(plain_text, ".")
    
    # Print to console
    cat("\n", fig_label, "\n  ", plain_text, "\n", sep = "")
    
    # Start a new paragraph for abbreviations
    # We'll use fpar (formatted paragraph) to mix italic and regular text
    fp_normal <- fp_text(font.family = "Arial", font.size = 11, italic = FALSE)
    fp_italic <- fp_text(font.family = "Arial", font.size = 11, italic = TRUE)
    
    # Paragraph properties for justified text
    fp_par_justified <- fp_par(text.align = "justify", line_spacing = 1)
    
    # Build the formatted text runs
    text_runs <- list()
    
    for (i in seq_along(abbrev_entries)) {
      # Add abbreviation (italic)
      text_runs[[length(text_runs) + 1]] <- ftext(abbrev_entries[[i]]$abbrev, prop = fp_italic)
      # Add comma and space
      text_runs[[length(text_runs) + 1]] <- ftext(", ", prop = fp_normal)
      # Add longform (normal)
      text_runs[[length(text_runs) + 1]] <- ftext(abbrev_entries[[i]]$longform, prop = fp_normal)
      
      # Add semicolon separator (except for last entry)
      if (i < length(abbrev_entries)) {
        text_runs[[length(text_runs) + 1]] <- ftext("; ", prop = fp_normal)
      } else {
        # Add period at the end
        text_runs[[length(text_runs) + 1]] <- ftext(".", prop = fp_normal)
      }
    }
    
    # Create formatted paragraph with all text runs and justified alignment
    formatted_para <- do.call(fpar, c(text_runs, list(fp_p = fp_par_justified)))
    
    # Add the formatted paragraph with proper spacing
    doc <- doc |>
      body_add_fpar(formatted_para) |>
      body_add_par("", style = "Normal")  # Add blank line after each figure
  }
  
  # Ensure output directory exists
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the document
  print(doc, target = output_path)
  
  message("Figure abbreviations document created: ", output_path)
  
  invisible(output_path)
}
