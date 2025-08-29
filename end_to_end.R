#* End to End Script Running and Compilation of All Files
#+ Get Full Tree Structure of Repo
  list_tree(getwd())
#+ Create a Compiled File of all Utilities and Scripts
  combine_R_files("R/Utilities", "Compilation/Combined_Utilities_All.R")
  combine_R_files("R/Scripts", "Compilation/Combined_Scripts_All.R")
  #!!1 rerun
#+ Run End-to-End Pipeline
  #- Set path to raw data
    raw_path <- "/Users/JoshsMacbook2015/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Manuscripts and Projects/Active Projects/Thyroid Metabolomics/Raw Data"
  #- Load progress if needed
    obj_names <- load_checkpoints(raw_path)

    # Inspect what got loaded
    print(obj_names)
  #- Run pipeline for specified files (or entire repo)
    #! Input the number(s) you want it to run to (00-06 = full)
    run_scripts("00-06", raw_path)

