#' Solve Math Expressions in Render Files
#'
#' Crawls through render_figures.R (or specified file) and evaluates all
#' math expressions in x=, y=, width=, height= parameters, replacing them
#' with their exact decimal values.
#'
#' @param file Path to file to process. If NULL, automatically finds any file 
#'   with "render_figure" in its name in R/Scripts/
#' @export
#' @examples
#' m()  # Auto-find and process render_figures.R
#' m("R/Scripts/14_render_supplementary_figures.R")  # Process specific file
m <- function(file = NULL) {
  
  # If no file specified, find any file with "render_figure" in its name
  if (is.null(file)) {
    render_files <- list.files("R/Scripts", pattern = "render_figure.*\\.R$", 
                               full.names = TRUE, ignore.case = TRUE)
    
    if (length(render_files) == 0) {
      stop("No render_figure*.R files found in R/Scripts/")
    } else if (length(render_files) > 1) {
      cat("Multiple render_figure files found:\n")
      for (i in seq_along(render_files)) {
        cat(sprintf("  [%d] %s\n", i, render_files[i]))
      }
      stop("Please specify which file to process")
    }
    
    file <- render_files[1]
    cat("Processing:", file, "\n")
  }
  
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  
  # Read the file
  lines <- readLines(file, warn = FALSE)
  
  # Pattern 1: x =, y =, width =, height = followed by math expressions
  # Allow spaces in expression, be more explicit
  pattern1 <- "(width|height|x|y)(\\s*=\\s*)([0-9.]+\\s*[+*/\\-]\\s*[0-9.]+(?:\\s*[+*/\\-]\\s*[0-9.]+)*)(\\s*[,)])"
  
  # Pattern 2: expressions inside c() like c(0.76, 8.75-0.25) - expression in 2nd position
  # Made more permissive to capture complex expressions
  pattern2 <- "c\\(([0-9.]+),\\s*([0-9.+*/()\\-\\s]+)\\)"
  
  # Pattern 3: expressions inside c() like c(0.76-100/800, 10.00) - expression in 1st position
  # Made more permissive to capture complex expressions
  pattern3 <- "c\\(([0-9.+*/()\\-\\s]+),\\s*([0-9.]+)\\)"
  
  # Debug: Check for any lines with division
  debug_lines <- grep("/", lines)
  if (length(debug_lines) > 0) {
    cat("Lines with division found:", debug_lines, "\n")
  }
  
  modified <- FALSE
  
  for (i in seq_along(lines)) {
    # Keep replacing until no more matches with math operators found
    repeat {
      old_line <- lines[i]
      line_changed <- FALSE
      
      # Try pattern 1: x=, y=, width=, height=
      if (grepl(pattern1, lines[i])) {
        # Get ALL matches on this line
        temp_line <- lines[i]
        match_found <- FALSE
        
        while (grepl(pattern1, temp_line)) {
          match <- regmatches(temp_line, regexec(pattern1, temp_line))[[1]]
          
          if (length(match) >= 5) {
            param <- match[2]
            equals <- match[3]
            expr <- match[4]
            comma <- match[5]
            
            # Check if this expression has operators
            if (grepl("[+*/\\-]", expr)) {
              # Try to evaluate the expression
              tryCatch({
                result <- eval(parse(text = expr))
                
                if (is.numeric(result)) {
                  result_str <- format(result, digits = 10, nsmall = 0)
                  result_str <- gsub("\\s+", "", result_str)
                  
                  # Build replacement
                  replacement <- paste0(param, equals, result_str, comma)
                  
                  # Replace just this match in actual line
                  lines[i] <- sub(paste0(param, equals, expr, comma), replacement, lines[i], fixed = TRUE)
                  
                  cat(sprintf("Line %d: %s = %s  →  %s\n", i, param, expr, result_str))
                  modified <- TRUE
                  line_changed <- TRUE
                  match_found <- TRUE
                  break  # Process one at a time, will loop back
                }
              }, error = function(e) {})
            }
          }
          
          # Remove this match from temp_line to check for more
          temp_line <- sub(pattern1, "", temp_line)
        }
        
        if (match_found) next  # Continue repeat loop
      }
      
      # Try pattern 2: c() expressions with expression in 2nd position
      if (!line_changed && grepl(pattern2, lines[i])) {
        match <- regmatches(lines[i], regexec(pattern2, lines[i]))[[1]]
        
        if (length(match) >= 3) {
          first_val <- match[2]
          expr <- match[3]
          
          # Check if this expression has operators between numbers (not just a leading minus)
          if (grepl("[0-9]\\s*[+*/\\-]\\s*[0-9]", expr)) {
            # Try to evaluate the expression
            tryCatch({
              result <- eval(parse(text = expr))
              
              if (is.numeric(result)) {
                result_str <- format(result, digits = 10, nsmall = 0)
                result_str <- gsub("\\s+", "", result_str)
                
                # Only replace if result is different from original
                if (result_str != gsub("\\s+", "", expr)) {
                  # Build replacement
                  old_c <- paste0("c(", first_val, ", ", expr, ")")
                  new_c <- paste0("c(", first_val, ", ", result_str, ")")
                  
                  # Replace just this match in actual line
                  lines[i] <- sub(old_c, new_c, lines[i], fixed = TRUE)
                  
                  cat(sprintf("Line %d: c(..., %s)  →  c(..., %s)\n", i, expr, result_str))
                  modified <- TRUE
                  line_changed <- TRUE
                }
              }
            }, error = function(e) {})
          }
        }
      }
      
      # Try pattern 3: c() expressions with expression in 1st position
      if (!line_changed && grepl(pattern3, lines[i])) {
        match <- regmatches(lines[i], regexec(pattern3, lines[i]))[[1]]
        
        if (length(match) >= 3) {
          expr <- match[2]
          second_val <- match[3]
          
          # Check if this expression has operators between numbers (not just a leading minus)
          if (grepl("[0-9]\\s*[+*/\\-]\\s*[0-9]", expr)) {
            # Try to evaluate the expression
            tryCatch({
              result <- eval(parse(text = expr))
              
              if (is.numeric(result)) {
                result_str <- format(result, digits = 10, nsmall = 0)
                result_str <- gsub("\\s+", "", result_str)
                
                # Only replace if result is different from original
                if (result_str != gsub("\\s+", "", expr)) {
                  # Build replacement
                  old_c <- paste0("c(", expr, ", ", second_val, ")")
                  new_c <- paste0("c(", result_str, ", ", second_val, ")")
                  
                  # Replace just this match in actual line
                  lines[i] <- sub(old_c, new_c, lines[i], fixed = TRUE)
                  
                  cat(sprintf("Line %d: c(%s, ...)  →  c(%s, ...)\n", i, expr, result_str))
                  modified <- TRUE
                  line_changed <- TRUE
                }
              }
            }, error = function(e) {})
          }
        }
      }
      
      # If line didn't change, break out of repeat loop
      if (!line_changed) break
    }
  }
  
  if (modified) {
    writeLines(lines, file)
    cat(sprintf("\n✓ Updated %s\n", file))
  } else {
    cat("No math expressions found to evaluate.\n")
  }
  
  invisible(modified)
}
