#' Build Supplementary Table 1 LaTeX
#'
#' Converts a tibble to a landscape, multi-page LaTeX longtable with:
#' - Arial font, size 7pt
#' - Proper LaTeX formatting for column names (italics, subscripts)
#' - Multi-page spillover with repeated headers
#' - 0.5pt borders
#' - Hierarchical structure: Group (bold) and Metabolite (indented plain text)
#'
#' @param data The prepared tibble with Display_Name column containing hierarchy
#' @param caption Table caption (appears above table)
#' @param footnote Table footnote (appears below table)
#' @param output_path Path to write the .tex file
#'
#' @return Invisibly returns the LaTeX code
#'
build_ST1_latex <- function(data, 
                           caption = NULL, 
                           footnote = NULL,
                           output_path = NULL) {
  
  library(kableExtra)
  library(dplyr)
  
  # Format column names with LaTeX syntax
  # Display_Name -> \textbf{Group}: Metabolite Name (bold Group, plain Metabolite Name)
  # P Value -> \textit{P} Value
  # q Value -> \textit{q} Value
  # log2(FC) -> log\textsubscript{2}(FC)
  # mz -> \textit{m/z}
  format_latex_colnames <- function(names) {
    names <- gsub("^Display_Name$", "\\\\textbf{Group}: Metabolite Name", names)
    names <- gsub("^P Value$", "\\\\textit{P} Value", names)
    names <- gsub("^q Value$", "\\\\textit{q} Value", names)
    names <- gsub("^log2\\(FC\\)$", "log\\\\textsubscript{2}(FC)", names)
    names <- gsub("^mz$", "\\\\textit{m/z}", names)
    names
  }

  # Apply formatting to column names
  formatted_names <- format_latex_colnames(names(data))
  
  # Determine alignment: first two columns (Display_Name, Subgroup) left, rest centered
  n_cols <- ncol(data)
  col_align <- c("l", "l", rep("c", n_cols - 2))
  
  # Build longtable with kableExtra
  latex_table <- kbl(
    data,
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    escape = FALSE,  # Allow LaTeX in data and column names
    col.names = formatted_names,
    align = col_align
  ) |>
    kable_styling(
      latex_options = c("repeat_header"),
      font_size = 7
    )  # Don't bold entire header row - we control header formatting via col.names
  
  # Extract the raw LaTeX
  latex_code <- as.character(latex_table)
  
  # Remove any rowcolor commands to ensure no striping
  latex_code <- gsub("\\\\rowcolor\\{[^}]+\\}", "", latex_code)
  
  # Set border thickness to 0.5pt (matching example)
  latex_code <- gsub("\\\\toprule(?!\\[)", "\\\\toprule[0.5pt]", latex_code, perl = TRUE)
  latex_code <- gsub("\\\\midrule(?!\\[)", "\\\\midrule[0.5pt]", latex_code, perl = TRUE)
  
  # Add bottomrule to endfoot for consistent page bottom borders
  latex_code <- gsub("\\\\endfoot", "\\\\bottomrule\n\\\\endfoot", latex_code)
  
  # Split into lines for processing
  lines <- strsplit(latex_code, "\n")[[1]]
  
  # Process data rows to apply formatting
  # Only bold GROUP header rows - metabolites with leading spaces pass through unchanged
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Skip header and structural lines
    if (grepl("^\\\\(toprule|midrule|bottomrule|endfirsthead|endhead|endfoot|addlinespace|rule|begin|end|caption)", line)) next
    
    # Process lines that contain data (have & separator)
    if (grepl(" & ", line)) {
      parts <- strsplit(line, " & ", fixed = TRUE)[[1]]
      if (length(parts) >= 1) {
        first_col <- parts[1]
        rest <- if(length(parts) > 1) paste(parts[-1], collapse = " & ") else ""
        
        # Check if this is a GROUP header (starts with letter, no leading spaces, not empty)
        # GROUP rows have no leading spaces in the data
        if (grepl("^[A-Za-z]", first_col) && !grepl("^ ", first_col) && 
            nchar(trimws(first_col)) > 0) {
          # GROUP: bold
          first_col <- paste0("\\textbf{", first_col, "}")
        }
        # Metabolite rows (with leading spaces) and spacer rows pass through unchanged
        
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
  
  # Add caption if provided
  if (!is.null(caption)) {
    latex_code <- gsub(
      "(\\\\begin\\{longtable\\}\\{[^}]+\\})",
      paste0("\\1\n\\\\caption{", caption, "}\\\\\\\\"),
      latex_code
    )
  }
  
  # Add footnote if provided  
  if (!is.null(footnote)) {
    footnote_latex <- paste0(
      "\n\\\\multicolumn{", ncol(data), "}{l}{\\\\footnotesize ", 
      footnote, 
      "} \\\\\\\\"
    )
    latex_code <- gsub(
      "\\\\end\\{longtable\\}",
      paste0(footnote_latex, "\n\\\\end{longtable}"),
      latex_code
    )
  }
  
  # Write to file if path provided
  if (!is.null(output_path)) {
    writeLines(latex_code, output_path)
    message("LaTeX table written to: ", output_path)
  }
  
  invisible(latex_code)
}
