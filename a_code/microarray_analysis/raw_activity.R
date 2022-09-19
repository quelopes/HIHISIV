raw_activity <- function(way_raw, way_out, pheno, method, name_exp) {
  print("*** ACTIVITY: NORMALIZATION ***")

  setwd(paste0(way_raw, project, "/CEL"))
  list_dataset <- pheno_data$SAMPLE_NAME

  # raw = ReadAffy()
  # === NORMALIZATION ===
  # --- if mas5 ---
  if (method == "mas5") {
    raw <- ReadAffy(filenames = list_dataset)
    eset_mas5 <- mas5(raw) # few minutes...
    expr_set_nologs <- exprs(eset_mas5)
    expr_set <- log(expr_set_nologs, 2)
    # summary(exprSet)
  }

  # --- if rma ---
  if (method == "rma") {
    raw <- ReadAffy(filenames = list_dataset)
    eset_rma <- rma(raw) # few minutes...
    expr_set_nologs <- exprs(eset_rma)
    expr_set <- log(expr_set_nologs, 2)
    # summary(exprSet)
  }

  if (method == "jac") {
    #   setwd(way_raw)
    mData <- list.files()
    # Merge dataset
    for (file in mData) {
      # if the merged dataset doesn't exist, create it
      if (!exists("dataset")) {
        dataset <- read.table(file, header = TRUE, sep = "\t", fill = TRUE, quote = "", check.names = FALSE)
        dataset <- dataset[, c(1, 2)]
        colnames(dataset) <- c("probe", file)
        name1 <- file
      }

      # if the merged dataset does exist, append to it
      if (exists("dataset")) {
        temp_dataset <- read.table(file, header = TRUE, sep = "\t")
        temp_dataset <- temp_dataset[, c(1, 2)]
        colnames(temp_dataset) <- c("probe", file)
        dataset <- full_join(dataset, temp_dataset, by = "probe")
        # dataset = merge(dataset, temp_dataset, by.x="probe",by.y="probe", all=TRUE)
        rm(temp_dataset)
      }
    }
    dataset <- dataset[, -2]
    colnames(dataset)[2] <- name1
    row.dataset <- as.character(dataset[, 1])
    dataset <- dataset[, -1]
    rownames(dataset) <- row.dataset

    # 	Exclude lines with NA
    percent <- 0.95 # 70% lines with NA excluded
    missing.percent.row <- apply(dataset[, 1:ncol(dataset)], 1, function(x) {
      (length(na.omit(x)) / length(x))
    })
    ind <- which(missing.percent.row <= percent)
    if (length(ind) > 0) {
      dataset <- dataset[-ind, ]
    }

    data.set <- as.matrix(log2(dataset))
    data.set <- impute.knn(data.set, k = 10, rowmax = 0.6, colmax = 1)

    dataSet <- normalizeBetweenArrays(data.set$data, method = "cyclicloess")
    #   method = "quantile"
    # 		summary(dataSet)

    expr_set <- dataSet
    #   setwd(set)
  }

  # --- writing data (table) ---
  write.table(expr_set, file = paste0(way_out, "/", name_exp, "/", name_exp, "-MatrixNormalized.csv"), quote = F, sep = ",")
  expr_set
}