# Mummichog Analysis Parameters

**Analysis Date:** 2026-01-13

**Database:** hsa_kegg

**MetaboAnalystR 'Set' Function Outputs:**
- SetPeakFormat: mpr
- SetPeakEnrichMethod: mummichog (v2)
- SetMummichogPval: 0.05 (user-specified threshold)

**Instrument Parameters (UpdateInstrumentParameters):**
- instrumentOpt: 5
- msModeOpt: mixed
- force_primary_ion: yes
- rt_frac: 0.02

**Analysis Parameters:**
- Peak filtering method: User-specified p-value threshold
- Peak filtering threshold (rounded): 0.05
- Peak filtering threshold (precise): 0.05
- Peaks analyzed: 552 out of 43851
- Pathways analyzed: 69
- Significant pathways (p < 0.05): 5
- Pathway p-values range: 0.009305 to 0.99999
- Pathway FDR: Not calculated (using raw p-values)
- Pathway enrichment FDR threshold: 0.05 (fixed)
- Minimum pathway size: 3
- Background permutations: 100

**Input Data:**
- Number of features: 43851
- Output directory: Outputs/mummichog/correlations_KEGG/LTD4

