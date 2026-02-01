# Mummichog Analysis Parameters

**Analysis Date:** 2026-02-01

**Database:** hsa_mfn

**MetaboAnalystR 'Set' Function Outputs:**
- SetPeakFormat: mpr
- SetPeakEnrichMethod: mummichog (v2)
- SetMummichogPvalFromPercent: 0.1 (top 10% of peaks)

**Instrument Parameters (UpdateInstrumentParameters):**
- instrumentOpt: 5
- msModeOpt: mixed
- force_primary_ion: yes
- rt_frac: 0.02

**Analysis Parameters:**
- Peak filtering method: Top 10% of peaks (dynamic)
- Peak filtering threshold (rounded): 0.25
- Peak filtering threshold (precise): 0.411735054725732
- Peaks analyzed: 864 out of 17270
- Pathways analyzed: 79
- Significant pathways (p < 0.05): 8
- Pathway p-values range: 0.007522 to 0.98954
- Pathway FDR: Not calculated (using raw p-values)
- Pathway enrichment FDR threshold: 0.05 (fixed)
- Minimum pathway size: 3
- Background permutations: 100

**Input Data:**
- Number of features: 17270
- Output directory: Outputs/mummichog/ttests_MFN

