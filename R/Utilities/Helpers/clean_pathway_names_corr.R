#' Clean and Standardize Pathway Names for Correlation Analysis
#'
#' This function is a specialized version of clean_pathway_names() designed for
#' correlation analysis visualizations. It includes all standard transformations
#' plus additional abbreviations and custom replacements for specific pathways.
#'
#' @param pathway_names Character vector of pathway names to clean and standardize
#'
#' @return Character vector of cleaned pathway names with standardized formatting
#'
#' @details
#' The function applies all transformations from clean_pathway_names() plus:
#' - Custom pathway name replacements for better visualization
#' - Final abbreviations: "Biosynthesis" -> "Biosynth.", "Metabolism" -> "Met."
#'
#' @examples
#' \dontrun{
#'   pathway_names <- c(
#'     "Aminosugars Metabolism",
#'     "Drug Metabolism - Cyp450"
#'   )
#'   clean_pathway_names_corr(pathway_names)
#' }
#'
#' @seealso \code{\link{clean_pathway_names}} for standard pathway name cleaning
#'
#' @export
clean_pathway_names_corr <- function(pathway_names) {
  pathway_names |>
    # Specific pathway exceptions FIRST - before any other transformations
    stringr::str_replace_all(stringr::regex("\\bPhenylalanine, tyrosine and tryptophan biosynthesis\\b", ignore_case = TRUE), "Phenylalanine, Tyrosine & Tryptophan Biosynthesis") |>
    stringr::str_replace_all(stringr::regex("\\bValine, leucine and isoleucine degradation\\b", ignore_case = TRUE), "Branched-Chain Amino Acid Degradation") |>
    stringr::str_replace_all(stringr::regex("\\bUrea cycle/amino group metabolism\\b", ignore_case = TRUE), "Urea Cycle/Amino Group Metabolism") |>
    # Handle Greek letters and special cases first - EXPANDED
    stringr::str_replace_all(stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE), "β-Alanine") |>
    stringr::str_replace_all(stringr::regex("\\balpha[- ]?linolenic\\b", ignore_case = TRUE), "α-Linolenic") |>
    stringr::str_replace_all(stringr::regex("\\bgamma[- ]?linol", ignore_case = TRUE), "γ-linol") |>
    stringr::str_replace_all(stringr::regex("\\bgama[- ]?linoleic\\b", ignore_case = TRUE), "γ-Linoleic") |>
    # Additional Greek letter replacements
    stringr::str_replace_all(stringr::regex("\\bbeta[- ]?oxidation\\b", ignore_case = TRUE), "β-Oxidation") |>
    stringr::str_replace_all(stringr::regex("\\balpha[- ]?", ignore_case = TRUE), "α-") |>
    stringr::str_replace_all(stringr::regex("\\bbeta[- ]?", ignore_case = TRUE), "β-") |>
    stringr::str_replace_all(stringr::regex("\\bgamma[- ]?", ignore_case = TRUE), "γ-") |>
    stringr::str_replace_all(stringr::regex("\\bdelta[- ]?", ignore_case = TRUE), "δ-") |>
    stringr::str_replace_all(stringr::regex("\\bepsilon[- ]?", ignore_case = TRUE), "ε-") |>
    stringr::str_replace_all(stringr::regex("\\bomega[- ]?", ignore_case = TRUE), "ω-") |>
    # Specific fatty acid replacements
    stringr::str_replace_all(stringr::regex("\\bPolyunsaturated Fatty Acid\\b", ignore_case = TRUE), "PUFA") |>
    stringr::str_replace_all(stringr::regex("\\bdi-unsaturated fatty acid\\b", ignore_case = TRUE), "Di-UFA") |>
    stringr::str_replace_all(stringr::regex("\\bDi-UFAs\\b"), "Di-UFA") |>
    stringr::str_replace_all(stringr::regex("\\bmono-unsaturated fatty acid\\b", ignore_case = TRUE), "MUFA") |>
    stringr::str_replace_all(stringr::regex("\\bMono-UFAs\\b"), "MUFA") |>
    # Specific terpenoid-quinone replacement
    stringr::str_replace_all(stringr::regex("\\bterpenoid-quinone\\b", ignore_case = TRUE), "TQ") |>
    # Fix specific multi-word phrases BEFORE individual word capitalizations
    stringr::str_replace_all(stringr::regex("\\bamino sugar\\s*&\\s*nucleotide sugar\\b", ignore_case = TRUE), "Amino & Nucleotide Sugar") |>
    # Vitamin formatting - keep ONLY parenthetical content
    stringr::str_replace_all(stringr::regex("\\bVitamin A \\(Retinol\\)\\b", ignore_case = TRUE), "Retinol") |>
    stringr::str_replace_all(stringr::regex("\\bVitamin B1 \\(thiamin\\)\\b", ignore_case = TRUE), "Thiamin") |>
    stringr::str_replace_all(stringr::regex("\\bVitamin D3 \\(cholecalciferol\\)\\b", ignore_case = TRUE), "Cholecalciferol") |>
    # Amino acid abbreviations (excluding when part of beta-alanine or other compounds)
    stringr::str_replace_all(stringr::regex("\\Ala(?!\\s)", ignore_case = TRUE), "Alanine") |>
    stringr::str_replace_all(stringr::regex("\\Asp\\b", ignore_case = TRUE), "Aspartate") |>
    stringr::str_replace_all(stringr::regex("\\bHis\\b", ignore_case = TRUE), "Histidine") |>
    stringr::str_replace_all(stringr::regex("\\bacids\\b", ignore_case = TRUE), "Acid") |>
    # Fatty acid abbreviations
    stringr::str_replace_all(stringr::regex("\\bunsaturated fatty acids?\\b", ignore_case = TRUE), "UFAs") |>
    # Special abbreviations and replacements
    stringr::str_replace_all(stringr::regex("\\bcytochrome P450\\b", ignore_case = TRUE), "Cyp450") |>
    stringr::str_replace_all(stringr::regex("\\bGlycosylphosphatidylinositol \\(GPI\\)-anchor\\b", ignore_case = TRUE), "GPI-anchor") |>
    stringr::str_replace_all(stringr::regex("\\bThiamine\\b", ignore_case = TRUE), "Vitamin B1") |>
    # Add Oxford commas for 3+ item lists before replacing 'and' with '&'
    stringr::str_replace_all("Val, Leu  &  Ile Degradation", "Val, Leu, & Ile Degradation") |>
    stringr::str_replace_all("Val, Leu  &  Ile", "Val, Leu, & Ile") |>
    stringr::str_replace_all(stringr::regex("\\bVal,\\s*Leu\\s+&\\s+Ile\\b", ignore_case = TRUE), "Val, Leu, & Ile") |>
    stringr::str_replace_all(stringr::regex("\\bVal,\\s*Leu\\s*&\\s*Ile\\b", ignore_case = TRUE), "Val, Leu, & Ile") |>
    # Replace 'and' with '&'
    stringr::str_replace_all(stringr::regex("\\band\\b", ignore_case = FALSE), "&") |>
    # Universal capitalizations - specific patterns first
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\s*&\\s*metabolism\\b", ignore_case = TRUE), "Biosynthesis/Metabolism") |>
    stringr::str_replace_all(stringr::regex("\\bde\\s+novo\\s+fatty\\b", ignore_case = TRUE), "De Novo Fatty") |>
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\s+oxidation\\b", ignore_case = TRUE), "Peroxisomal Oxidation") |>
    stringr::str_replace_all(stringr::regex("\\bcarnitine\\s+shuttle\\b", ignore_case = TRUE), "Carnitine Shuttle") |>
    stringr::str_replace_all(stringr::regex("\\bmetabolism\\b", ignore_case = TRUE), "Metabolism") |>
    stringr::str_replace_all(stringr::regex("\\bdegradation\\b", ignore_case = TRUE), "Degradation") |>
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\b", ignore_case = TRUE), "Biosynthesis") |>
    stringr::str_replace_all(stringr::regex("\\bshuttle\\b", ignore_case = TRUE), "Shuttle") |>
    stringr::str_replace_all(stringr::regex("\\bnovo\\b", ignore_case = TRUE), "Novo") |>
    stringr::str_replace_all(stringr::regex("\\bfatty\\b", ignore_case = TRUE), "Fatty") |>
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\b", ignore_case = TRUE), "Peroxisomal") |>
    stringr::str_replace_all(stringr::regex("\\bformation\\b", ignore_case = TRUE), "Formation") |>
    stringr::str_replace_all(stringr::regex("\\barachidonate\\b", ignore_case = TRUE), "Arachidonate") |>
    stringr::str_replace_all(stringr::regex("\\bdihomo\\b", ignore_case = TRUE), "Dihomo") |>
    stringr::str_replace_all(stringr::regex("\\bactivation\\b", ignore_case = TRUE), "Activation") |>
    stringr::str_replace_all(stringr::regex("\\bmannose\\b", ignore_case = TRUE), "Mannose") |>
    stringr::str_replace_all(stringr::regex("\\bacid\\b", ignore_case = TRUE), "Acid") |>
    stringr::str_replace_all(stringr::regex("\\boxidation\\b", ignore_case = TRUE), "Oxidation") |>
    stringr::str_replace_all(stringr::regex("\\bperoxisome\\b", ignore_case = TRUE), "Peroxisome") |>
    stringr::str_replace_all(stringr::regex("\\bretinol\\b", ignore_case = TRUE), "Retinol") |>
    stringr::str_replace_all(stringr::regex("\\bbile\\b", ignore_case = TRUE), "Bile") |>
    stringr::str_replace_all(stringr::regex("\\bhormone\\b", ignore_case = TRUE), "Hormone") |>
    # Handle "other" before "other enzymes" to catch standalone cases
    stringr::str_replace_all(stringr::regex("\\bother\\b", ignore_case = TRUE), "Other") |>
    stringr::str_replace_all(stringr::regex("\\bOther enzymes\\b"), "Other Enzymes") |>
    stringr::str_replace_all(stringr::regex("\\bnicotinamide\\b", ignore_case = TRUE), "Nicotinamide") |>
    stringr::str_replace_all(stringr::regex("\\bnicotinate\\b", ignore_case = TRUE), "Nicotinate") |>
    stringr::str_replace_all(stringr::regex("\\bcycle\\b", ignore_case = TRUE), "Cycle") |>
    stringr::str_replace_all(stringr::regex("\\bfolate\\b", ignore_case = TRUE), "Folate") |>
    stringr::str_replace_all(stringr::regex("\\bsteroid\\b", ignore_case = TRUE), "Steroid") |>
    stringr::str_replace_all(stringr::regex("\\bprimary\\b", ignore_case = TRUE), "Primary") |>
    stringr::str_replace_all(stringr::regex("\\bsphingolipid\\b", ignore_case = TRUE), "Sphingolipid") |>
    stringr::str_replace_all(stringr::regex("\\bporphyrin\\b", ignore_case = TRUE), "Porphyrin") |>
    stringr::str_replace_all(stringr::regex("\\bglucuronate\\b", ignore_case = TRUE), "Glucuronate") |>
    stringr::str_replace_all(stringr::regex("\\binterconversions\\b", ignore_case = TRUE), "Interconversions") |>
    stringr::str_replace_all(stringr::regex("\\bascorbate\\b", ignore_case = TRUE), "Ascorbate") |>
    stringr::str_replace_all(stringr::regex("\\baldarate\\b", ignore_case = TRUE), "Aldarate") |>
    stringr::str_replace_all(stringr::regex("\\bdicarboxylate\\b", ignore_case = TRUE), "Dicarboxylate") |>
    # NEW CAPITALIZATIONS - Added per user request
    stringr::str_replace_all(stringr::regex("\\bligation\\b", ignore_case = TRUE), "Ligation") |>
    stringr::str_replace_all(stringr::regex("\\bnucleotide\\b", ignore_case = TRUE), "Nucleotide") |>
    stringr::str_replace_all(stringr::regex("\\bsugar\\b", ignore_case = TRUE), "Sugar") |>
    stringr::str_replace_all(stringr::regex("\\bphosphorylation\\b", ignore_case = TRUE), "Phosphorylation") |>
    stringr::str_replace_all(stringr::regex("\\bsucrose\\b", ignore_case = TRUE), "Sucrose") |>
    stringr::str_replace_all(stringr::regex("\\bproline\\b", ignore_case = TRUE), "Proline") |>
    stringr::str_replace_all(stringr::regex("\\basparagine\\b", ignore_case = TRUE), "Asparagine") |>
    stringr::str_replace_all(stringr::regex("\\belongation\\b", ignore_case = TRUE), "Elongation") |>
    # FINAL EXPLICIT FIXES - After all other transformations
    stringr::str_replace_all("Val, Leu  &  Ile Degradation", "Val, Leu, & Ile Degradation") |>
    stringr::str_replace_all("Val, Leu & Ile Degradation", "Val, Leu, & Ile Degradation") |>
    stringr::str_replace_all("Amino Sugar & Nucleotide Sugar", "Amino & Nucleotide Sugar") |>
    stringr::str_replace_all("Amino Sugar  &  Nucleotide Sugar", "Amino & Nucleotide Sugar") |>
    # CUSTOM REPLACEMENTS FOR CORRELATION ANALYSIS
    stringr::str_replace_all("^Aminosugars Metabolism$", "Aminosugar Metabolism") |>
    stringr::str_replace_all("^Drug Metabolism - Cyp450$", "CYP450 \\(Drug\\) Metabolism") |>
    stringr::str_replace_all("^Keratan sulfate Degradation$", "Keratan sulfate Degradation") |>
    stringr::str_replace_all("^Vitamin A \\(Retinol\\) Metabolism$", "Vitamin A Metabolism") |>
    stringr::str_replace_all("^Methionine & cysteine Metabolism$", "Methionine & Cysteine Metabolism") |>
    stringr::str_replace_all("^Squalene & cholesterol Biosynthesis$", "Squalene & Cholesterol Biosynthesis") |>
    stringr::str_replace_all("^Branched-Chain Amino Acid Degradation$", "BCAA Degradation") |>
    stringr::str_replace_all("^Glycosphingolipid Biosynthesis - globoseries$", "Glycosphingolipid Biosynthesis \\(Globoseries\\)") |>
    stringr::str_replace_all("^Glycosphingolipid Biosynthesis - lactoseries$", "Glycosphingolipid Biosynthesis \\(Lactoseries\\)") |>
    stringr::str_replace_all("^Glycosphingolipid Biosynthesis - ganglioseries$", "Glycosphingolipid Biosynthesis \\(Ganglioseries\\)") |>
    stringr::str_replace_all("^Glycine, serine, alanine & threonine Metabolism$", "Gly/Ser/Ala/Thr Metabolism") |>
    stringr::str_replace_all("^Glycosphingolipid Biosynthesis - neolactoseries$", "Glycosphingolipid Biosynthesis \\(Neolactoseries\\)") |>
    stringr::str_replace_all("^Putative anti-Inflammatory metabolites Formation from EPA$", "EPA Anti-Inflammatory")
}
