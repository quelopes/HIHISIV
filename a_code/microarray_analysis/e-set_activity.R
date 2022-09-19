e_set_activity <- function(mat, pheno, annot) {
  print("*** ACTIVITY e-SET ***")

  rownames(pheno) <- pheno$SAMPLE_NAME
  pheno_data <- new("AnnotatedDataFrame", data = pheno)
  samples <- as.character(pheno$SAMPLE_NAME)
  mat_sel <- mat[, samples]

  colnames(mat_sel) <- pheno$SAMPLE_NAME
  e_set <- new("ExpressionSet", exprs = mat_sel, phenoData = pheno_data, annotation = annot)

  e_set
}
