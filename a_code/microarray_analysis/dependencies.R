dependencies <- function(code_way) {

  # ================================
  # ---------- PACKAGES ------------
  # ================================
  print("*** STARTING LIBRARIES AND FUNCTIONS ***")

  library(affy)
  library(impute)
  library(dplyr)
  library(limma)
  library(GOstats) 
  library(org.Hs.eg.db) 
  library(org.Mmu.eg.db) 
  library(ggplot2)
  library(RColorBrewer)
  
  library(biomaRt)

  # ================================
  # --------- FUNCTIONS ------------
  # ================================
  source(paste0(code_way, "raw_activity.R"))
  source(paste0(code_way, "e-set_activity.R"))
  source(paste0(code_way, "limma_activity.R"))
  source(paste0(code_way, "go_term_analysis.R"))
}


