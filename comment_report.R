#! 00a_environment_setup.R
#* 0a: Environment Setup
#+ 0a.1: Load conflicted and set ALL preferences BEFORE loading other packages
#+ 0a.2: Load all packages from DESCRIPTION file
#+ 0a.3: Load GitHub Packages explicitly for renv detection
#+ 0a.4: Check system dependencies

#! 00b_setup.R
#* 0b: Configuration Setup for DTC Metabolomics Analysis
#+ 0b.1: Load utility functions first (needed for dynamic config)
#+ 0b.2: Load dynamic project configuration 
#+ 0b.3: Set up global paths from config 
#+ 0b.4: Set up R environment preferences 
#+ 0b.5: Set random seed for reproducibility
#+ 0b.6: Import All Clinical and Path Data
#- 0b.6.1: Tumor IDs
#- 0b.6.2: Tumor Path
#- 0b.6.3: Load MetaboJanitoR processed CSV files  
#- 0b.6.4: Read in QC that goes with annotated TFT; deduplicate with lower p-value

#! 00c_clinical_metadata.R
#* 0c: Clinical metadata processing
#+ 0c.1: Cleanup tumor pathology data; compute T stage; compute overall stage
#- 0c.1.1: Create full version
#- 0c.1.2: Create cleanup version

#! 00d_FTs.R
#* 0d: Cleanup
#+ 0d.1: Targeted Processing
#- 0d.1.1: Remove QC from TFTs; tag with IDs; join with clinical metadata
#+ 0d.2: Untargeted Processing
#- 0d.2.1: Remove QC from UFTs (full); tag with IDs; join with clinical metadata
#- 0d.2.2: Remove QC from UFTs (filtered); tag with IDs; join with clinical metadata; apply remove_qc)
#+ 0d.3: Define any possible feature columns in vector

#! 01_clustering.R
#* 1: PCA and PLS-DA Analysis (Computation Only)
#+ 1.1: Run Volcano Analysis on UFT data (Computation Only)
#+ 1.2: Create Plots
#- 1.2.6: Volcano

#! 02_pathway_enrichment.R
#* 2: Pathway Enrichment Analysis
#+ 2.1: Run Mummichog Analysis
#- 2.1.1: Run Mummichog Ttest function
#- 2.1.2: Run Mummichog (MFN Only)
#- 2.1.3: Run Mummichog (KEGG)
#+ 2.3: Create Pathway Enrichment Plots (MFN)
#- 2.3.1: Define JSON file paths once
#- 2.3.2: Import tibbles for inspection
#- 2.3.3: Make MFN inspection tibble for visualization
#- 2.3.4: Inspect combined for visualization
#- 2.3.5: Print min and max of MFN and combined p values and enrichment for visualization scaling
#- 2.3.6: Make MFN only plot
#+ 2.4: Create Pathway Enrichment Plots (KEGG)
#- 2.4.1: Make KEGG inspection tibble for visualization
#- 2.4.2: Print min and max of KEGG p values and enrichment for visualization scaling
#- 2.4.3: Make KEGG only plot
#+ 2.5: Run Biological Network Analysis (MFN)
#+ 2.6: Plot Biological Networks (MFN)

#! 03_annotated.R
#* 3: Annotated Plots
#+ 3.1: T-tests for all metabolomic features against PGD status
#- 3.1.1: Transform feature table
#- 3.1.2: Run targeted t-tests
#- 3.1.3: Add metadata to annotated results
#- 3.1.4: Export as Excel for QC
#+ 3.2: Create Diverging Bar Plots
#- 3.2.1: Prepare data for diverging bars
#- 3.2.2: Create Plot
#+ 3.3: Create individual feature plots
#- 3.3.1: Subset to features for scatter plots
#- 3.3.2: Select only scatter plot features
#- 3.3.3: Create relevant scatter plots
#+ 3.4: Make Correlations Matrix (Full Version)
#- 3.4.1: Subset to selected features for matrix
#- 3.4.2: Create data for correlation matrix
#- 3.4.3: Create correlation matrix plot
#- 3.4.4: Inspect tibble of results
#+ 3.5: Make Correlations Matrix (Curated Version)
#- 3.5.1: Subset to selected features for matrix
#- 3.5.2: Create data for correlation matrix
#- 3.5.3: Create correlation matrix plot

#! 05_assign_plots.R
#* 5: Assign and Render Plots
#+ 5.1: Figure 1
#+ 5.2: Figure 2
#+ 5.3: Figure 3

#! 06_render_figures.R
#* 6: Render Figures
#+ 6.1: Figure 1
#+ 6.2: Figure 2
#+ 6.3: Figure 3
#+ 6.4: Figure 4
#+ 6.5: Print all figures

#! 07_data_not_shown.R
#* 7: Abstract Data Generation
#+ 7.1: Variant Type Distribution
#- 7.1.0: Join FT with path data
#- 7.1.1: Compute variant counts and percentages
#- 7.1.2: Create variant distribution sentence
#- 7.1.3: Early vs Advanced stage variant counts
#+ 7.2: Stage Distribution  
#- 7.2.1: Compute stage summary
#- 7.2.2: Individual stage counts
#- 7.2.3: Early vs Advanced stage grouping
#- 7.2.4: Create stage distribution sentence with early vs advanced grouping
#+ 7.3: Metabolite Feature Counts
#- 7.3.2: Count QC-filtered features by chromatography method  
#+ 7.4: Volcano Plot Statistics
#- 7.4.1: Extract volcano plot results (assuming you have volcano analysis results)
#+ 7.5: Stage Binning Summary (Early vs Advanced)
#- 7.5.1: Compute early vs advanced stage distribution
#+ 7.6: Display Manuscript Sentences
#+ 7.7: Summary Statistics Table
#- 7.7.1: Create summary statistics for abstract/manuscript

#! 08_tables.R
#* 8: Tables
#+ 8.1: Table 1
#- 8.1.0: Create vector of included samples
#- 8.1.1: Filter to used samples; order rows; clean data
#- 8.1.2: Build Table 1

