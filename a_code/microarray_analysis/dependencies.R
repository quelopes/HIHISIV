dependencies <- function(code_way) {

  # ================================
  # ---------- PACKAGES ------------
  # ================================
  print("*** STARTING LIBRARIES AND FUNCTIONS ***")

  library(affy)
  library(impute)
  library(dplyr)
  
  library(limma)

  # # ----- PLOTS -----
  library(ggplot2)
  library(RColorBrewer)

  # ================================
  # --------- FUNCTIONS ------------
  # ================================
  source(paste0(code_way, "module_processing.R"))
  source(paste0(code_way, "raw_activity.R"))
  source(paste0(code_way, "e-set_activity.R"))
  source(paste0(code_way, "limma_activity.R"))
}


