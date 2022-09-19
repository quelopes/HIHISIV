# ======================
# === Configurations ===
# ======================
project <- "2932572286"
experiment_name <- "2932572286_acute_chronic"
org_annot <- "org.Hs.eg.db"

# --- ways ---
raw_data <- "../microarray/"
way_result <- "../microarray_results/"
way_code <- "../microarray_analysis/"


# --- create a dir with the results ---
if (file.exists(experiment_name)) {
  setwd(file.path(way_result, experiment_name))
  print("the directory already exists")
} else {
  dir.create(file.path(way_result, experiment_name))
  setwd(file.path(way_result, experiment_name))
}

# --- normalization ---
# norm_method <- "mas5"
norm_method <- "jac"

# --- e-set ---
pheno_data <- read.table(paste0(raw_data, project, "/", experiment_name, ".csv"), head = T, sep = ",")

# ======================
# === GERAL SETTING ===
# ======================
# --- limma ---
fdr <- 0.05 # threshold
log_fc <- 1 # threshold

source(paste0(way_code, "dependencies.R"))
dependencies(way_code)
module_processing(
  data_raw = raw_data, result_way = way_result, name_experiment = experiment_name,
  fdr = fdr, log_fc = log_fc, method_norm = norm_method, pheno_data = pheno_data,
  org_annot = org_annot
)

