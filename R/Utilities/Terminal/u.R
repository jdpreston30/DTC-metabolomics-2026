#' Reload All Utility Functions and Config
#'
#' Quick helper to re-source all utility functions from R/Utilities/ and reload config
#' Excludes Setup/ directory (one-time setup scripts, not utilities)
#' Useful during development when making changes to utility functions or config
#'
#' @export
#' @examples
#' u() # Update utilities and config
u <- function() {
  # Load config
  source("R/Utilities/Helpers/load_dynamic_config.R")
  config <<- load_dynamic_config(computer = "auto", config_path = "all_run/config_dynamic.yaml")
  
  # Load all utilities (excluding Setup/ directory)
  utils_path <- "R/Utilities/"
  all_files <- list.files(utils_path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE)
  
  # Exclude only Setup/ directory (one-time setup scripts)
  utility_files <- all_files[!grepl("Setup/", all_files)]
  
  purrr::walk(utility_files, source)
  
  cat("✅ Config and utilities reloaded\n")
  invisible(NULL)
}

#' Update, Run Script, and Render Figures
#'
#' Convenience wrapper that updates utilities/config, runs a numbered script, then renders figures
#'
#' @param n Script number to run (e.g., 1 for 01_*.R, 5 for 05_*.R)
#' @export
#' @examples
#' uf(1) # Update, run 01_clustering.R, then render figures
#' uf(5) # Update, run 05_render_figures.R, then render figures
uf <- function(n) {
  u()
  r(n)
  f()
  invisible(NULL)
}

#' Run Render Figures Script
#'
#' Quick helper to run assign_plots.R and render_figures.R scripts
#' Dynamically finds these scripts regardless of their prefix numbers
#' Useful during development for quickly re-rendering figures
#'
#' @export
#' @examples
#' f()
f <- function() {
  scripts_dir <- "R/Scripts/"
  
  # First run assign_plots (find any file with assign_plots in the name)
  assign_pattern <- "assign_plots\\.R$"
  assign_files <- list.files(scripts_dir, pattern = assign_pattern, full.names = TRUE, ignore.case = TRUE)
  
  if (length(assign_files) == 0) {
    cat("⚠️  Script not found: assign_plots.R\n")
    return(invisible(NULL))
  }
  
  if (length(assign_files) > 1) {
    cat("⚠️  Multiple assign_plots scripts found:\n")
    cat("   ", paste(basename(assign_files), collapse = ", "), "\n")
    return(invisible(NULL))
  }
  
  assign_script <- assign_files[1]
  cat("▶️  Running:", basename(assign_script), "\n")
  source(assign_script)

  # Then run render_figures (find any file with render_figures in the name)
  render_pattern <- "render_figures\\.R$"
  render_files <- list.files(scripts_dir, pattern = render_pattern, full.names = TRUE, ignore.case = TRUE)
  
  if (length(render_files) == 0) {
    cat("⚠️  Script not found: render_figures.R\n")
    return(invisible(NULL))
  }
  
  if (length(render_files) > 1) {
    cat("⚠️  Multiple render_figures scripts found:\n")
    cat("   ", paste(basename(render_files), collapse = ", "), "\n")
    return(invisible(NULL))
  }
  
  render_script <- render_files[1]
  cat("▶️  Running:", basename(render_script), "\n")
  source(render_script)
  cat("✅ Figures rendered\n")

  invisible(NULL)
}

#' Run Numbered Script
#'
#' Quick helper to run any numbered script by its number
#' Automatically pads single digit numbers with leading zero
#'
#' @param n Script number (e.g., 1 for 01_FTs.R, 5 for 05_targeted_volcano_diverge.R)
#' @export
#' @examples
#' r(1) # Runs 01_FTs.R
#' r(5) # Runs 05_targeted_volcano_diverge.R
#' r(9) # Runs 06_render_figures.R
r <- function(n) {
  # Pad single digit numbers with leading zero
  script_num <- sprintf("%02d", n)

  # Find script file that matches the pattern
  scripts_dir <- "R/Scripts/"
  pattern <- paste0("^", script_num, ".*\\.R$")

  matching_files <- list.files(scripts_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)

  if (length(matching_files) == 0) {
    cat("⚠️  No script found matching number:", n, "\n")
    return(invisible(NULL))
  }

  if (length(matching_files) > 1) {
    cat("⚠️  Multiple scripts found matching number:", n, "\n")
    cat("   ", paste(basename(matching_files), collapse = ", "), "\n")
    return(invisible(NULL))
  }

  script_path <- matching_files[1]
  cat("▶️  Running:", basename(script_path), "\n")
  source(script_path)
  cat("✅ Script completed\n")

  invisible(NULL)
}

#' @rdname u
#' @export
U <- u