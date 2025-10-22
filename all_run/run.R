{
  source("R/Utilities/Helpers/load_dynamic_config.R")
  config <- load_dynamic_config(computer = "auto", config_path = "all_run/config_dynamic.yaml")
  source("R/Scripts/00a_environment_setup.R")
  source("R/Scripts/00b_setup.R")
  source("R/Scripts/00c_clinical_metadata.R")
  source("R/Scripts/00d_FTs.R")
  source("R/Scripts/01_Clustering.R")
  source("R/Scripts/02_pathway_enrichment.R")
  source("R/Scripts/03_annotated.R")
  source("R/Scripts/04_assign_plots.R")
  source("R/Scripts/05_visualization.R")
  source("R/Scripts/06_abstract_data.R")
}
