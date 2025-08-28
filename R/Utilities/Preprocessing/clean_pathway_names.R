clean_pathway_names <- function(pathway_names) {
  pathway_names %>%
    # Handle Greek letters and special cases first
    stringr::str_replace_all(stringr::regex("\\bBeta[- ]Alanine\\b", ignore_case = TRUE), "β-Alanine") %>%
    stringr::str_replace_all(stringr::regex("\\balpha[- ]?linolenic\\b", ignore_case = TRUE), "α-Linolenic") %>%
    stringr::str_replace_all(stringr::regex("\\bgamma[- ]?linol", ignore_case = TRUE), "γ-linol") %>%
    stringr::str_replace_all(stringr::regex("\\bgama[- ]?linoleic\\b", ignore_case = TRUE), "γ-Linoleic") %>%
    # Vitamin formatting - keep ONLY parenthetical content
    stringr::str_replace_all(stringr::regex("\\bVitamin A \\(Retinol\\)\\b", ignore_case = TRUE), "Retinol") %>%
    stringr::str_replace_all(stringr::regex("\\bVitamin B1 \\(thiamin\\)\\b", ignore_case = TRUE), "Thiamin") %>%
    stringr::str_replace_all(stringr::regex("\\bVitamin D3 \\(cholecalciferol\\)\\b", ignore_case = TRUE), "Cholecalciferol") %>%
    # Amino acid abbreviations (excluding when part of beta-alanine or other compounds)
    stringr::str_replace_all(stringr::regex("\\bAlanine(?!\\s)", ignore_case = TRUE), "Ala") %>%
    stringr::str_replace_all(stringr::regex("(?<!β-)\\bAlanine\\b", ignore_case = TRUE), "Ala") %>%
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
    stringr::str_replace_all(stringr::regex("\\bTryptophan\\b", ignore_case = TRUE), "Trp") %>%
    stringr::str_replace_all(stringr::regex("\\bTyrosine\\b", ignore_case = TRUE), "Tyr") %>%
    stringr::str_replace_all(stringr::regex("\\bValine\\b", ignore_case = TRUE), "Val") %>%
    stringr::str_replace_all(stringr::regex("\\bTaurine\\b", ignore_case = TRUE), "Tau") %>%
    stringr::str_replace_all(stringr::regex("\\bHypotaurine\\b", ignore_case = TRUE), "Hypotau") %>%
    # Fatty acid abbreviations
    stringr::str_replace_all(stringr::regex("\\bunsaturated fatty acids?\\b", ignore_case = TRUE), "UFAs") %>%
    # Prostaglandin abbreviation
    stringr::str_replace_all(stringr::regex("\\bProstaglandin\\b", ignore_case = TRUE), "PG") %>%
    # Special abbreviations and replacements
    stringr::str_replace_all(stringr::regex("\\bcytochrome P450\\b", ignore_case = TRUE), "Cyp450") %>%
    stringr::str_replace_all(stringr::regex("\\bGlycosylphosphatidylinositol \\(GPI\\)-anchor\\b", ignore_case = TRUE), "GPI-anchor") %>%
    stringr::str_replace_all(stringr::regex("\\bThiamine\\b", ignore_case = TRUE), "Vitamin B1") %>%
    # Replace 'and' with '&'
    stringr::str_replace_all(stringr::regex("\\band\\b", ignore_case = FALSE), " & ") %>%
    # Universal capitalizations - specific patterns first
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\s*&\\s*metabolism\\b", ignore_case = TRUE), "Biosynthesis/Metabolism") %>%
    stringr::str_replace_all(stringr::regex("\\bde\\s+novo\\s+fatty\\b", ignore_case = TRUE), "De Novo Fatty") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\s+oxidation\\b", ignore_case = TRUE), "Peroxisomal Oxidation") %>%
    stringr::str_replace_all(stringr::regex("\\bcarnitine\\s+shuttle\\b", ignore_case = TRUE), "Carnitine Shuttle") %>%
    stringr::str_replace_all(stringr::regex("\\bmetabolism\\b", ignore_case = TRUE), "Metabolism") %>%
    stringr::str_replace_all(stringr::regex("\\bdegradation\\b", ignore_case = TRUE), "Degradation") %>%
    stringr::str_replace_all(stringr::regex("\\bbiosynthesis\\b", ignore_case = TRUE), "Biosynthesis") %>%
    stringr::str_replace_all(stringr::regex("\\bshuttle\\b", ignore_case = TRUE), "Shuttle") %>%
    stringr::str_replace_all(stringr::regex("\\bnovo\\b", ignore_case = TRUE), "Novo") %>%
    stringr::str_replace_all(stringr::regex("\\bfatty\\b", ignore_case = TRUE), "Fatty") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisomal\\b", ignore_case = TRUE), "Peroxisomal") %>%
    stringr::str_replace_all(stringr::regex("\\bformation\\b", ignore_case = TRUE), "Formation") %>%
    stringr::str_replace_all(stringr::regex("\\barachidonate\\b", ignore_case = TRUE), "Arachidonate") %>%
    stringr::str_replace_all(stringr::regex("\\bdihomo\\b", ignore_case = TRUE), "Dihomo") %>%
    stringr::str_replace_all(stringr::regex("\\bactivation\\b", ignore_case = TRUE), "Activation") %>%
    stringr::str_replace_all(stringr::regex("\\bmannose\\b", ignore_case = TRUE), "Mannose") %>%
    stringr::str_replace_all(stringr::regex("\\bacid\\b", ignore_case = TRUE), "Acid") %>%
    stringr::str_replace_all(stringr::regex("\\boxidation\\b", ignore_case = TRUE), "Oxidation") %>%
    stringr::str_replace_all(stringr::regex("\\bperoxisome\\b", ignore_case = TRUE), "Peroxisome") %>%
    stringr::str_replace_all(stringr::regex("\\bretinol\\b", ignore_case = TRUE), "Retinol") %>%
    stringr::str_replace_all(stringr::regex("\\bbile\\b", ignore_case = TRUE), "Bile") %>%
    stringr::str_replace_all(stringr::regex("\\bhormone\\b", ignore_case = TRUE), "Hormone") %>%
    stringr::str_replace_all(stringr::regex("\\bother enzymes\\b", ignore_case = TRUE), "Other Enzymes") %>%
    stringr::str_replace_all(stringr::regex("\\bnicotinamide\\b", ignore_case = TRUE), "Nicotinamide") %>%
    stringr::str_replace_all(stringr::regex("\\bnicotinate\\b", ignore_case = TRUE), "Nicotinate") %>%
    stringr::str_replace_all(stringr::regex("\\bcycle\\b", ignore_case = TRUE), "Cycle") %>%
    stringr::str_replace_all(stringr::regex("\\bfolate\\b", ignore_case = TRUE), "Folate") %>%
    stringr::str_replace_all(stringr::regex("\\bsteroid\\b", ignore_case = TRUE), "Steroid") %>%
    stringr::str_replace_all(stringr::regex("\\bprimary\\b", ignore_case = TRUE), "Primary") %>%
    stringr::str_replace_all(stringr::regex("\\bsphingolipid\\b", ignore_case = TRUE), "Sphingolipid") %>%
    stringr::str_replace_all(stringr::regex("\\bporphyrin\\b", ignore_case = TRUE), "Porphyrin") %>%
    stringr::str_replace_all(stringr::regex("\\bglucuronate\\b", ignore_case = TRUE), "Glucuronate") %>%
    stringr::str_replace_all(stringr::regex("\\binterconversions\\b", ignore_case = TRUE), "Interconversions") %>%
    stringr::str_replace_all(stringr::regex("\\bascorbate\\b", ignore_case = TRUE), "Ascorbate") %>%
    stringr::str_replace_all(stringr::regex("\\baldarate\\b", ignore_case = TRUE), "Aldarate") %>%
    stringr::str_replace_all(stringr::regex("\\bdicarboxylate\\b", ignore_case = TRUE), "Dicarboxylate")
}
