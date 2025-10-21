#' Clean Pathway Names for Network Visualization
#'
#' This function applies specialized text transformations to pathway names specifically
#' for network visualization contexts where shorter, more compact names are preferred.
#' Includes extensive abbreviations and formatting optimized for node labels.
#'
#' @param pathway_names Character vector of pathway names to clean and abbreviate
#'
#' @return Character vector of cleaned pathway names with network-appropriate abbreviations
#'
#' @details
#' This function differs from the standard clean_pathway_names() by applying more
#' aggressive abbreviations suitable for network visualizations:
#' - Full amino acid names -> 3-letter codes (Alanine -> Ala)
#' - Long pathway names -> abbreviated forms (Polyunsaturated fatty acid -> PUFA)
#' - Complex names -> compact versions with line break considerations
#' - Ampersands and "and" -> forward slashes for compactness
#' - Specialized abbreviations for vitamins, fatty acids, and metabolic processes
#'
#' @examples
#' \dontrun{
#'   pathway_names <- c(
#'     "Valine, leucine and isoleucine degradation",
#'     "Polyunsaturated fatty acid biosynthesis",
#'     "Drug metabolism - cytochrome P450"
#'   )
#'   clean_pathway_names_for_network(pathway_names)
#'   # Returns: c("BCAA Degradation", 
#'   #            "PUFA Biosynthesis",
#'   #            "Drug Metabolism (Cyp450)")
#' }
#'
#' @seealso \code{\link{clean_pathway_names}} for standard pathway name cleaning
#' @export
clean_pathway_names_for_network <- function(pathway_names) {
    pathway_names %>%
      # Handle Greek letters and special cases first
      stringr::str_replace_all(stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE), "β-Ala") %>%
      stringr::str_replace_all(stringr::regex("\\bBile acid\\b", ignore_case = TRUE), "Bile Acid") %>%
      stringr::str_replace_all(stringr::regex("\\balpha[- ]?linolenic\\b", ignore_case = TRUE), "α-Linolenic") %>%
      stringr::str_replace_all(stringr::regex("\\bgamma[- ]?linol", ignore_case = TRUE), "γ-linol") %>%
      stringr::str_replace_all(stringr::regex("\\bbeta[- ]?oxidation\\b", ignore_case = TRUE), "β-Oxidation") %>%
      stringr::str_replace_all(stringr::regex("\\balpha[- ]?", ignore_case = TRUE), "α-") %>%
      stringr::str_replace_all(stringr::regex("\\bbeta[- ]?", ignore_case = TRUE), "β-") %>%
      stringr::str_replace_all(stringr::regex("\\bgamma[- ]?", ignore_case = TRUE), "γ-") %>%
      stringr::str_replace_all(stringr::regex("\\bdelta[- ]?", ignore_case = TRUE), "δ-") %>%
      stringr::str_replace_all(stringr::regex("\\bomega[- ]?", ignore_case = TRUE), "ω-") %>%
      # Amino acid abbreviations (excluding when part of other compounds)
      stringr::str_replace_all(stringr::regex("\\bAlanine\\b", ignore_case = TRUE), "Ala") %>%
      stringr::str_replace_all(stringr::regex("\\bArginine\\b", ignore_case = TRUE), "Arg") %>%
      stringr::str_replace_all(stringr::regex("\\bAsparagine\\b", ignore_case = TRUE), "Asn") %>%
      stringr::str_replace_all(stringr::regex("\\bAspartate\\b", ignore_case = TRUE), "Asp") %>%
      stringr::str_replace_all(stringr::regex("\\bCysteine\\b", ignore_case = TRUE), "Cys") %>%
      stringr::str_replace_all(stringr::regex("\\bGlutamate\\b", ignore_case = TRUE), "Glu") %>%
      stringr::str_replace_all(stringr::regex("\\bGlutamine\\b", ignore_case = TRUE), "Gln") %>%
      stringr::str_replace_all(stringr::regex("\\bGlycine\\b", ignore_case = TRUE), "Gly") %>%
      stringr::str_replace_all(stringr::regex("\\bHistidine\\b", ignore_case = TRUE), "His") %>%
      stringr::str_replace_all(stringr::regex("\\bIsoleucine\\b", ignore_case = TRUE), "Ile") %>%
      stringr::str_replace_all(stringr::regex("\\bLeucine\\b", ignore_case = TRUE), "Leu") %>%
      stringr::str_replace_all(stringr::regex("\\bLysine\\b", ignore_case = TRUE), "Lys") %>%
      stringr::str_replace_all(stringr::regex("\\bMethionine\\b", ignore_case = TRUE), "Met") %>%
      stringr::str_replace_all(stringr::regex("\\bPhenylalanine\\b", ignore_case = TRUE), "Phe") %>%
      stringr::str_replace_all(stringr::regex("\\bProline\\b", ignore_case = TRUE), "Pro") %>%
      stringr::str_replace_all(stringr::regex("\\bSerine\\b", ignore_case = TRUE), "Ser") %>%
      stringr::str_replace_all(stringr::regex("\\bThreonine\\b", ignore_case = TRUE), "Thr") %>%
      # stringr::str_replace_all(stringr::regex("\\bTryptophan\\b", ignore_case = TRUE), "Trp") %>%
      # stringr::str_replace_all(stringr::regex("\\bTyrosine\\b", ignore_case = TRUE), "Tyr") %>%
      stringr::str_replace_all(stringr::regex("\\bValine\\b", ignore_case = TRUE), "Val") %>%
      # Specific long pathway abbreviations
      stringr::str_replace_all(stringr::regex("\\bProstaglandin formation from arachidonate\\b", ignore_case = TRUE), "PG Formation (Arachidonate)") %>%
      stringr::str_replace_all(stringr::regex("\\bArachidonic acid metabolism\\b", ignore_case = TRUE), "Arachidonic Acid Metabolism") %>%
      stringr::str_replace_all(stringr::regex("\\bPolyunsaturated fatty acid\\b", ignore_case = TRUE), "PUFA") %>%
      stringr::str_replace_all(stringr::regex("\\bC21-steroid hormone\\b", ignore_case = TRUE), "C21-Steroid Hormone") %>%
      stringr::str_replace_all(stringr::regex("\\bPolyunsaturated FA\\b", ignore_case = TRUE), "PUFA") %>%
      stringr::str_replace_all(stringr::regex("\\bMono-unsaturated fatty acid\\b", ignore_case = TRUE), "MUFA") %>%
      stringr::str_replace_all(stringr::regex("\\bDi-unsaturated fatty acid\\b", ignore_case = TRUE), "Di-UFA") %>%
      stringr::str_replace_all(stringr::regex("\\bSaturated fatty acids?\\b", ignore_case = TRUE), "SFA") %>%
      stringr::str_replace_all(stringr::regex("\\bOmega-3 fatty acid\\b", ignore_case = TRUE), "ω-3 FA") %>%
      stringr::str_replace_all(stringr::regex("\\bPhytanic acid peroxisomal\\b", ignore_case = TRUE), "Phytanic Acid Perox.") %>%
      stringr::str_replace_all(stringr::regex("\\bPutative anti-Inflammatory metabolites formation from EPA\\b", ignore_case = TRUE), "Anti-Inflam. From EPA") %>%
      stringr::str_replace_all(stringr::regex("\\bDrug metabolism - cytochrome P450\\b", ignore_case = TRUE), "Drug Metabolism (Cyp450)") %>%
      stringr::str_replace_all(stringr::regex("\\bFructose and mannose\\b", ignore_case = TRUE), "Fructose/Mannose") %>%
      stringr::str_replace_all(stringr::regex("\\bValine, leucine and isoleucine\\b", ignore_case = TRUE), "BCAA") %>%
      stringr::str_replace_all(stringr::regex("\\bvaline, leucine & isoleucine\\b", ignore_case = TRUE), "BCAA") %>%
      stringr::str_replace_all(stringr::regex("\\bVal, Leu and Ile\\b", ignore_case = TRUE), "BCAA") %>%
      stringr::str_replace_all(stringr::regex("\\bVal, Leu & Ile\\b", ignore_case = TRUE), "BCAA") %>%
      # Urea Cycle specific replacement (must come BEFORE Met. replacement)
      stringr::str_replace_all(stringr::regex("\\bUrea cycle/amino group metabolism\\b", ignore_case = TRUE), "Urea Cycle") %>%
      # Vitamin formatting - specific replacements
      stringr::str_replace_all(stringr::regex("\\bVitamin B6 \\(pyridoxine\\) metabolism\\b", ignore_case = TRUE), "Vit. B6(Pyridoxine) Metabolism") %>%
      stringr::str_replace_all(stringr::regex("\\bVitamin A \\(retinol\\) metabolism\\b", ignore_case = TRUE), "Vitamin A Metabolism") %>%
      stringr::str_replace_all(stringr::regex("\\bVitamin B9 \\(folate\\) metabolism\\b", ignore_case = TRUE), "Vit. B9 (Folate) Metabolism") %>%
      stringr::str_replace_all(stringr::regex("\\bVitamin E metabolism\\b", ignore_case = TRUE), "Vit. E Metabolism") %>%
      # Fatty Acid abbreviation - catch all remaining variations
      stringr::str_replace_all(stringr::regex("\\bFatty acids?\\b", ignore_case = TRUE), "FA") %>%
      # Specific line breaks for long pathway names
      stringr::str_replace_all(stringr::regex("\\bFructose/Mannose Metabolism\\b", ignore_case = TRUE), "Fructose/Mannose Metabolism") %>%
      # Replace all ampersands with forward slashes
      stringr::str_replace_all(" & ", "/") %>%
      stringr::str_replace_all(" and ", "/") %>%
      # Ensure proper title case for key words
      stringr::str_replace_all(stringr::regex("\\bnitrogen\\b", ignore_case = TRUE), "Nitrogen") %>%
      stringr::str_replace_all(stringr::regex("\\bhexose\\b", ignore_case = TRUE), "Hexose") %>%
      stringr::str_replace_all(stringr::regex("\\bbutanoate\\b", ignore_case = TRUE), "Butanoate") %>%
      stringr::str_replace_all(stringr::regex("\\bubiquinone\\b", ignore_case = TRUE), "Ubiquinone") %>%
      stringr::str_replace_all(stringr::regex("\\bdrug\\b", ignore_case = TRUE), "Drug") %>%
      stringr::str_replace_all(stringr::regex("\\bperox\\.\\b", ignore_case = TRUE), "Perox.") %>%
      stringr::str_replace_all(stringr::regex("\\bde novo\\b", ignore_case = TRUE), "De Novo") %>%
      stringr::str_replace_all(stringr::regex("\\bshuttle\\b", ignore_case = TRUE), "Shuttle") %>%
      stringr::str_replace_all(stringr::regex("\\bactivation\\b", ignore_case = TRUE), "Activation") %>%
      stringr::str_replace_all(stringr::regex("\\bligation\\b", ignore_case = TRUE), "Ligation") %>%
      # Ensure proper capitalization for process terms
      stringr::str_replace_all(stringr::regex("\\bmetabolism\\b", ignore_case = TRUE), "Metabolism") %>%
      stringr::str_replace_all(stringr::regex("\\bphosphorylation\\b", ignore_case = TRUE), "Phosphorylation") %>%
      stringr::str_replace_all(stringr::regex("\\boxidation\\b", ignore_case = TRUE), "Oxidation") %>%
      stringr::str_replace_all(stringr::regex("\\bdegradation\\b", ignore_case = TRUE), "Degradation") %>%
      stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\b", ignore_case = TRUE), "Biosynthesis") %>%
      # Fix specific title case issues
      stringr::str_replace_all("anti-inflammatory", "Anti-Inflammatory") %>%
      stringr::str_replace_all("from", "From")
  }
