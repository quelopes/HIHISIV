module_processing <- function(data_raw, result_way, pheno_data, method_norm, name_experiment, org_annot, fdr, log_fc) {
  
  # === NORMALIZED ===
  mat_normalized <- raw_activity(data_raw, result_way, pheno_data, method_norm, name_experiment)

  # === E-SET ===
  e_set <- e_set_activity(mat_normalized, pheno_data, org_annot)

  # === LIMMA ===
  result_limma <- limma_activity(e_set, fdr, log_fc, result_way, name_experiment)
}