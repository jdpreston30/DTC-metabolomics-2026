#' Check System Dependencies
#'
#' Checks for required system-level dependencies and provides installation
#' instructions if any are missing. This function is non-blocking and will
#' only warn users about missing dependencies.
#'
#' @details
#' This function checks for the following system dependencies:
#' - **Ghostscript**: Required for PDF to EPS conversion (pdf2ps command)
#' - **ImageMagick**: Required by the magick package for image manipulation
#'
#' @return Invisibly returns a character vector of missing dependency names.
#'   If all dependencies are found, returns an empty character vector.
#'
#' @examples
#' \dontrun{
#' # Check all system dependencies
#' check_system_dependencies()
#' }
#'
#' @author Joshua D. Preston
#' @export
check_system_dependencies <- function() {
  cat("\n🔍 Checking system dependencies...\n")
  
  # Define all system dependencies with installation instructions
  system_deps <- list(
    list(
      name = "Ghostscript",
      command = "pdf2ps",
      purpose = "PDF to EPS conversion",
      install_mac = "brew install ghostscript",
      install_linux = "sudo apt-get install ghostscript  # or: sudo yum install ghostscript",
      install_windows = "Download from https://www.ghostscript.com/download/gsdnld.html",
      optional = FALSE
    ),
    list(
      name = "ImageMagick",
      command = "magick",
      purpose = "Image processing via magick package",
      install_mac = "brew install imagemagick",
      install_linux = "sudo apt-get install imagemagick  # or: sudo yum install ImageMagick",
      install_windows = "Download from https://imagemagick.org/script/download.php",
      optional = FALSE,
      check_function = function() {
        # ImageMagick can be available via R's magick package without system install
        if (requireNamespace("magick", quietly = TRUE)) {
          tryCatch({
            magick::magick_config()
            TRUE
          }, error = function(e) FALSE)
        } else {
          FALSE
        }
      }
    ),
    list(
      name = "TinyTeX/LaTeX",
      command = "xelatex",
      purpose = "PDF generation for supplementary tables",
      install_mac = "Run in R: tinytex::install_tinytex()",
      install_linux = "Run in R: tinytex::install_tinytex()",
      install_windows = "Run in R: tinytex::install_tinytex()",
      optional = FALSE,
      check_function = function() {
        # Check if tinytex is installed and functional
        if (requireNamespace("tinytex", quietly = TRUE)) {
          tinytex::tinytex_root() != ""
        } else {
          # Fallback: check if xelatex command exists
          result <- tryCatch({
            system2("xelatex", args = "--version", stdout = FALSE, stderr = FALSE)
          }, error = function(e) 1)
          result == 0
        }
      }
    )
  )
  
  missing_deps <- character(0)
  
  # Check each dependency
  for (dep in system_deps) {
    # Use custom check function if available, otherwise check command
    is_available <- if (!is.null(dep$check_function)) {
      dep$check_function()
    } else {
      # Check if command exists
      result <- tryCatch({
        system2("command", args = c("-v", dep$command), stdout = FALSE, stderr = FALSE)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)
      result == TRUE || result == 0
    }
    
    if (is_available) {
      cat("  ✅", dep$name, "found\n")
    } else {
      missing_deps <- c(missing_deps, dep$name)
      opt_str <- if (dep$optional) " (optional)" else " (required)"
      cat("  ⚠️ ", dep$name, opt_str, " not found\n", sep = "")
      cat("     Purpose:", dep$purpose, "\n")
      
      # Show OS-appropriate installation instructions
      os <- Sys.info()["sysname"]
      if (os == "Darwin") {
        cat("     Install: ", dep$install_mac, "\n", sep = "")
      } else if (os == "Linux") {
        cat("     Install: ", dep$install_linux, "\n", sep = "")
      } else if (os == "Windows") {
        cat("     Install: ", dep$install_windows, "\n", sep = "")
      }
      cat("\n")
    }
  }
  
  # Summary
  if (length(missing_deps) == 0) {
    cat("\n✅ All system dependencies found!\n\n")
  } else {
    required_missing <- system_deps[sapply(system_deps, function(d) {
      d$name %in% missing_deps && !d$optional
    })]
    
    if (length(required_missing) > 0) {
      cat("\n⚠️  Missing", length(required_missing), "required system dependencies.\n")
      cat("   Some features may not work until these are installed.\n\n")
    } else {
      cat("\n✅ All required dependencies found!\n")
      cat("   (Optional dependencies missing but not critical)\n\n")
    }
  }
  
  invisible(missing_deps)
}
