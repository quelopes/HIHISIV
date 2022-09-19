limma_activity <- function(e_set, p_value, fold_change, way_out, name_exp) {

    print("*** ACTIVITY DIFFERENTIAL EXPRESSION ***")

  # ---------------------
  # *** MATRIX DESIGN ***
  # ---------------------
  targets <- pData(e_set)

  f <- targets$SETS
  f <- factor(f)
  design <- model.matrix(~ -1 + f)

  # --- colunm names ---
  names_time <- make.names(unique(f))

  colnames(design) <- names_time
  fit <- lmFit(e_set, design)

  mat_exp_probe <- fit$coefficients # adjusted value

  diff <- NULL
  for (i in 1:(length(names_time))) {
    for (j in 2:(length(names_time))) {
      if (j > i) {
        diff <- c(diff, paste(names_time[j], names_time[i], sep = "-"))
      }
    }
    rm(j)
  }
  rm(i)
  cont_matrix <- makeContrasts(contrasts = diff, levels = design)

  # -------------------------------
  # *** FITTING LINEAR FUNCTION ***
  # -------------------------------
  fit2 <- contrasts.fit(fit, cont_matrix)
  fit2 <- eBayes(fit2)

  results <- decideTests(fit2)

  summary(results)
  result_fit <- summary(results)

  lines_expr_set <- nrow(exprs(e_set))

  # -------------------
  # *** DIFFERENCES ***
  # -------------------
  for (k in 1:ncol(cont_matrix)) {
    pdf(paste0(way_out, "/", name_exp, "/", name_exp, "-VolcanoPlot-", diff[k], ".pdf"))

    diff_val <- topTable(fit2, coef = diff[k], sort.by = "none", number = lines_expr_set, adjust.method = "BH") # , sort.by="none", adjust.method="BH", p.value=1, lfc=0)

    diff_val$odds <- exp(diff_val[, "B"])
    diff_val$odds.Prob <- (diff_val$odds) / (1 + diff_val$odds)

    # Selecting values (threshold = pvalue and logFC)
    ind <- with(diff_val, which(adj.P.Val <= p_value & abs(logFC) >= fold_change))
    diff_val_de <- diff_val[ind, ]

    if (is.character(diff_val[1, 1]) == TRUE) {
      probes <- as.character(diff_val_de$ID)
      probes2 <- as.character(diff_val$ID)
    }
    if (is.character(diff_val[1, 1]) == FALSE) {
      probes <- rownames(diff_val_de)
      probes2 <- rownames(diff_val)
    }

    # --- Volcano plot ---
    volcanoplot(fit2,
      coef = diff[k], highlight = 50,
      main = paste("Volcano", diff[k], sep = ": ")
    )
    abline(h = -log10(0.05), lty = 2)
    abline(v = -1.5, lty = 2, col = "green")
    abline(v = 1.5, lty = 2, col = "red")

    if (k == 1) {
      probes_unique <- probes
      probes_num <- length(probes)
      probes_unique_all <- probes2
    }
    if (k > 1) {
      probes_unique <- unique(c(probes_unique, probes))
      probes_num <- c(probes_num, length(probes))

      probes_unique_all <- unique(c(probes_unique_all, probes2))
    }
    dev.off()
  }

  # ---
  p_num <- as.data.frame(probes_num)
  p_numNames <- diff

  p_num <- t(p_num)

  # *** MERGE PROBES COM FC AND PVALUE SELECTED ***
  vec <- c(1:length(probes_unique))
  vec_all <- c(1:length(probes_unique_all))

  if (length(probes_unique) != 0) {
    list_de_all <- data.frame(ID = probes_unique, vec)
    list_all <- data.frame(ID = probes_unique_all, vec_all)

    for (m in 1:ncol(cont_matrix)) {
      mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
      if (is.character(diff_val[1, 1]) == TRUE) {
        list_de_all <- merge(list.de.all, mat[, c(1, 2, 6)], by = "ID", all.x = TRUE)
        mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
        list_all <- merge(list_all, mat[, c(1, 2, 6)], by = "ID", all.x = TRUE)
      }
      if (is.character(diff_val[1, 1]) == FALSE) {
        list_de_all <- merge(list_de_all, mat[, c(1, 5)], by.x = "ID", by.y = "row.names", all.x = TRUE)
        mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
        list_all <- merge(list_all, mat[, c(1, 5)], by.x = "ID", by.y = "row.names", all.x = TRUE)
      }
    }

    names <- as.character(list_de_all[, 1])
    list_de_all <- list_de_all[, -c(1, 2)]
    rownames(list_de_all) <- names
    selecionados <- rownames(list_de_all)
  }


  if (length(probes_unique) == 0) {
    print("No probes/genes Differentially Expressed!")
    list_all <- data.frame(ID = probes_unique_all, vec_all)
    for (m in 1:ncol(cont_matrix)) {
      mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
      if (is.character(diff_val[1, 1]) == TRUE) {
        mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
        list_all <- merge(list_all, mat[, c(1, 2, 6)], by = "ID", all.x = TRUE)
      }
      if (is.character(diff_val[1, 1]) == FALSE) {
        mat <- topTable(fit2, coef = diff[m], number = lines_expr_set)
        list_all <- merge(list_all, mat[, c(1, 5)], by.x = "ID", by.y = "row.names", all.x = TRUE)
      }
    }
    selecionados <- NULL
  }

  names2 <- as.character(list_all[, 1])
  list.all <- list_all[, -2]

  write.table(list_all,
    file = paste0(way_out, "/", name_exp, "/", name_exp, "_diffExpression.csv"),
    sep = ",", col.names = T, row.names = F
  )
  list_all
}